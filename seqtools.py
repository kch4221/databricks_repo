#!/usr/local/bin python
from presto.Sequence import reverseComplement, alignAssembly, AssemblyRecord
from anarci import run_anarci
import pandas as pd
from Bio.Seq import Seq
from Bio.Data import CodonTable
from presto.Defaults import *
def nuc_same_prob(nc1:str, nc2:str, p1, p2) -> float:
    # Calculate P( X1==X2 | Xo1, Xo2, Q1, Q2 ) 
    if 'N' not in (nc1, nc2):
        if nc1 == nc2:
            prob = p1*p2 + (1-p1)*(1-p2)/3
        else:
            prob = p2*(1-p1)/3 + p1*(1-p2)/3 + (1-p1)*(1-p2)/2        
    else:
        prob = 0.25           
    return prob

from pvalue import precomp2_05, precomp2_01, precomp2_001, precomp2_0001

def validateP(OesScore, pval=0.01, min_overlap=10, q=0.25)->bool:
    min_overlap = min(99, min_overlap)
    q = min(q, 0.49)
    if pval >= 1:
        return True
    elif pval >= 0.05:
        table_ptr = precomp2_05[min_overlap]
    elif pval >= 0.01:
        table_ptr = precomp2_01[min_overlap]
    elif pval >= 0.001:
        table_ptr = precomp2_001[min_overlap]
    else:
        table_ptr = precomp2_0001[min_overlap]
    cutoff = table_ptr[(int)(q * 100)]
    if (OesScore > cutoff): 
        return True
    return False


def nuc_prob_per_pos(nc1:str, p1:float, nc2:str=None, p2:float=None) -> dict:
    # Calcucate probability of 4 kinds of nucleotide base in each position
    if nc2 is None:
        if nc1 != 'N':
            err_prob = (1.0 - p1) / 3.0
            result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob}
            result[nc1] = p1
        else:
            err_prob = 0.25
            result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob}
    else:
        if 'N' not in (nc1, nc2):
            p_same = nuc_same_prob(nc1, nc2, p1, p2)
            if nc1 == nc2:
                err_prob = (1 - p1)*(1 - p2) / 3.0 / p_same
                result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob}
                result[nc1] = p1*p2 / p_same
            else:
                err_prob = (1 - p1)*(1 - p2) / 4 / p_same
                result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob}
                result[nc1] = p1*(1 - p2) / p_same / 3
                result[nc2] = p2*(1 - p1) / p_same / 3
        else:
            if nc1 == nc2:
                err_prob = 0.25
                result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob}
            else:
                if nc1=='N':
                    nc1 = nc2
                    p1 = p2
                err_prob = (1.0 - p1) / 3.0
                result = {'A':err_prob, 'T':err_prob, 'C':err_prob, 'G':err_prob} 
                result[nc1] = p1
    return result


def sequence_nuc_prob(chain1, prob1, chain2=None, prob2=None)->list:
    # Calcucate probability of 4 kinds of nucleotide base for a sequence
    if chain2 is None or not chain2:
        result = [nuc_prob_per_pos(nc1, p1) for nc1, p1 in zip(chain1, prob1)]
    else:
        result = [nuc_prob_per_pos(nc1, p1, nc2, p2) for nc1, p1, nc2, p2 in zip(chain1, prob1, chain2, prob2)]
    return result


def PearAssessmentScore(seq1, seq2, prob1, prob2, ignoreN=True):
    assert len(seq1) == len(seq2)

    result = []
    match_count = 0

    for nc1, nc2, p1, p2 in zip(list(seq1), list(seq2), prob1, prob2):
        if 'N' in (nc1, nc2):
            result.append(-0.75)
            if not ignoreN: match_count += 1
        elif nc1 == nc2 :
            result.append(nuc_same_prob(nc1, nc2, p1, p2))
            match_count += 1
        else:
            result.append(nuc_same_prob(nc1, nc2, p1, p2) - 1)

    return result, 1.0 - match_count / len(seq1)


def ObservedExpectedAlignmentScore(seq1, seq2, prob1, prob2) -> float:
    assert len(seq1) == len(seq2)

    result = 0
    for nc1, nc2, p1, p2 in zip(list(seq1), list(seq2), prob1, prob2):
        p_same = nuc_same_prob(nc1, nc2, p1, p2)
        result = result + p_same - (1 - p_same)

    return result


from presto.Defaults import *
from presto.Sequence import getDNAScoreDict, overlapConsensus

def pearAssemble(head_seq, tail_seq, alpha=0.01, max_error=default_assembly_max_error,
                  min_len=default_assembly_min_len, max_len=default_assembly_max_len, scan_reverse=False,
                  assembly_stats=None, score_dict=getDNAScoreDict(mask_score=(1, 1), gap_score=(0, 0)), pre_compute_i=None)-> AssemblyRecord:
    stitch = AssemblyRecord()
    head_str = str(head_seq.seq)
    tail_str = str(tail_seq.seq)
    head_len = len(head_str)
    tail_len = len(tail_str)
    head_prob = [1- 10**(qi/(-10.0)) for qi in head_seq.letter_annotations['phred_quality']]
    tail_prob = [1- 10**(qi/(-10.0)) for qi in tail_seq.letter_annotations['phred_quality']]
    
    # Fail if sequences are too short
    if head_len <= min_len or tail_len <= min_len:
        return stitch

    # Determine if quality scores are present
    has_quality = hasattr(head_seq, 'letter_annotations') and \
                  hasattr(tail_seq, 'letter_annotations') and \
                  'phred_quality' in head_seq.letter_annotations and \
                  'phred_quality' in tail_seq.letter_annotations
    
    if not has_quality:
        return stitch

    # Determine if sub-sequences are allowed and define scan range
    if scan_reverse and max_len >= min(head_len, tail_len):
        scan_len = head_len + tail_len - min_len
    else:
        scan_len = min(max(head_len, tail_len), max_len)

    # Iterate and score overlap segments
    if pre_compute_i is None or not pre_compute_i:
        for i in range(min_len, scan_len + 1):
            a = max(0, head_len - i)
            b = head_len - max(0, i - tail_len)
            x = max(0, i - head_len)
            y = min(tail_len, i)
            score, error = PearAssessmentScore(head_str[a:b], tail_str[x:y], head_prob[a:b], tail_prob[x:y])
            score = sum(score)
            # Save stitch as optimal if z-score improves
            if score > stitch.zscore:
                stitch.head_pos = (a, b)
                stitch.tail_pos = (x, y)
                stitch.zscore = score
                stitch.error = error
                stitch.evalue = i
    else:
        i = pre_compute_i
        a = max(0, head_len - i)
        b = head_len - max(0, i - tail_len)
        x = max(0, i - head_len)
        y = min(tail_len, i)
        score, error = PearAssessmentScore(head_str[a:b], tail_str[x:y], head_prob[a:b], tail_prob[x:y])
        score = sum(score)
        stitch.head_pos = (a, b)
        stitch.tail_pos = (x, y)
        stitch.zscore = score
        stitch.error = error
        stitch.evalue = i

    # Build stitched sequences and assign best_dict values
    if stitch.head_pos is not None:
        # Correct quality scores and resolve conflicts
        a, b = stitch.head_pos
        x, y = stitch.tail_pos
        # stitch.pvalue = pear_p_value(ObservedExpectedAlignmentScore(head_str[a:b], tail_str[x:y], head_prob[a:b], tail_prob[x:y]), omega=10)
        
        overlap_seq = overlapConsensus(head_seq[a:b], tail_seq[x:y])

        if a > 0 and y < tail_len:
            # Tail overlaps end of head
            stitch.seq = head_seq[:a] + overlap_seq + tail_seq[y:]
        elif b < head_len and x > 0:
            # Head overlaps end of tail
            stitch.seq = tail_seq[:x] + overlap_seq + head_seq[b:]
        elif a == 0 and b == head_len:
            # Head is a subsequence of tail
            stitch.seq = tail_seq[:x] + overlap_seq + tail_seq[y:]
        elif x == 0 and y == tail_len:
            # Tail is a subsequence of head
            stitch.seq = head_seq[:a] + overlap_seq + head_seq[b:]
        else:
            print('Invalid overlap condition for %s.' % head_seq.id)
        
        stitch.seq.id = head_seq.id if head_seq.id == tail_seq.id else '+'.join([head_seq.id, tail_seq.id])
        stitch.seq.name = stitch.seq.id
        stitch.seq.description = ''

        stitch.pvalue = validateP(ObservedExpectedAlignmentScore(head_str[a:b], tail_str[x:y], head_prob[a:b], tail_prob[x:y]), min_overlap=min_len, pval=alpha)
        
    return stitch


def nucProbability(head_seq, tail_seq, stitch):
    result = None
    head_str = str(head_seq.seq)
    tail_str = str(tail_seq.seq)
    head_len = len(head_str)
    tail_len = len(tail_str)
    head_prob = [1- 10**(qi/(-10.0)) for qi in head_seq.letter_annotations['phred_quality']]
    tail_prob = [1- 10**(qi/(-10.0)) for qi in tail_seq.letter_annotations['phred_quality']]
    a, b = stitch.head_pos
    x, y = stitch.tail_pos
    over_lap_prob = sequence_nuc_prob(head_seq[a:b], head_prob[a:b], tail_seq[x:y], tail_prob[x:y])
    if a > 0 and y < tail_len:
        result = sequence_nuc_prob(head_seq[:a], head_prob[:a]) + over_lap_prob + sequence_nuc_prob(tail_seq[y:], tail_prob[y:])
    elif b < head_len and x > 0:
        result = sequence_nuc_prob(tail_seq[:x], tail_prob[:x]) + over_lap_prob + sequence_nuc_prob(head_seq[b:], head_prob[b:])
    elif a == 0 and b == head_len:
        result = sequence_nuc_prob(tail_seq[:x], tail_prob[:x]) + over_lap_prob + sequence_nuc_prob(tail_seq[y:], tail_prob[y:])
    elif x == 0 and y == tail_len:
        result = sequence_nuc_prob(head_seq[:a], head_prob[:a]) + over_lap_prob + sequence_nuc_prob(head_seq[b:], head_prob[b:])
    return result


def reverseStitch(stitch, evalueAsProb=False, evalueAsReads=False):
    a, b = stitch.head_pos
    x, y = stitch.tail_pos
    stitch.seq = reverseComplement(stitch.seq)
    total_len = len(stitch.seq)
    overlap = stitch.overlap

    if a > 0 and y - x == overlap:
        stitch.head_pos = (0, overlap)
        stitch.tail_pos = (total_len - overlap - a, total_len - a)
    elif b - a > overlap and x > 0:
        stitch.head_pos = (total_len - x - overlap, total_len - x)
        stitch.tail_pos = (0, overlap)
    elif a == 0 and b == overlap:
        stitch.tail_pos = (total_len - y, total_len - y + overlap)
    elif x == 0 and y == overlap:
        stitch.head_pos = (total_len - y, total_len - y + overlap)

    if evalueAsProb:
        stitch.evalue.columns = [{'A':'T', 'T':'A', 'C':'G', 'G':'C'}[c] for c in stitch.evalue.columns]
        stitch.evalue = stitch.evalue.sort_index(ascending=False).reset_index(drop=True)
    elif evalueAsReads:
        stitch.evalue = [reverseComplement(read) for read in stitch.evalue]

    return stitch


def aa_prob_per_position(triple:pd.DataFrame, ignore_stop_condon=False):
    triple = triple.reset_index(drop=True)
    aa_prob = {}
    forward = CodonTable.unambiguous_rna_by_name['Standard'].forward_table

    for k, v in forward.items():
        k = k.replace('U', 'T')
        prob = triple.loc[0, k[0]] * triple.loc[1, k[1]] * triple.loc[2, k[2]]
        if v not in aa_prob.keys():
            aa_prob[v] = 0
        aa_prob[v] += prob

    if ignore_stop_condon:
        total = 1.0
        for trp in CodonTable.unambiguous_rna_by_name['Standard'].stop_codons:
            trp = trp.replace('U', 'T')
            total -= triple.loc[0, trp[0]] * triple.loc[1, trp[1]] * triple.loc[2, trp[2]]
        for k in aa_prob.keys():
            aa_prob[k] /= total

    return aa_prob


def aa_sequence_prob(nuc_prob:pd.DataFrame):
    triple = nuc_prob.reset_index(drop=True)
    result = []
    for i in range(triple.shape[0] // 3):
        result.append(aa_prob_per_position(triple=triple.iloc[i*3: i*3+3]))
    return pd.DataFrame(result, dtype='float16')


def translate_functional(seq, return_Nuc=False):
    seq = seq.replace('.', '').replace('-', '')
    aa = Seq(seq).translate()
    for i in range(2):
        if aa.count('*') == 0: break
    # while seq and aa.find('*') >=0 :
        seq = seq[1:]
        aa = Seq(seq).translate()
    if return_Nuc:
        return str(aa), seq[: len(aa)*3]
    return str(aa)


def translate_numbering(seq, functional=False):
    try:
        seq = seq.replace('.', '').replace('-', '')
        for i in range(3):
            aa = str(Seq(seq[i: ]).translate())
            if functional and aa.count('*') > 0:
                continue
            anarci_result = run_anarci(aa.replace('*', '').replace('X', ''))
            if not anarci_result[1][0] is None:
                return seq[i: i+len(aa)*3], anarci_result
        return None, None
    except Exception:
        return None, None