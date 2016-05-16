'''
Merge pairs of reads
'''

import operator, math
from Bio import SeqIO, SeqRecord, Seq

def merge(forward, reverse, max_diffs, size, size_var, output, stagger=False):
    '''
    forward : string (filepath) or filehandle
        fastq
    max_diffs : int
    size : int
    size_var : int
    '''

    min_size = size - size_var
    max_size = size + size_var

    for_recs = SeqIO.parse(forward, 'fastq')
    rev_recs = SeqIO.parse(reverse, 'fastq')

    out_recs = merged_reads(for_recs, rev_recs, max_diffs, min_size, max_size, stagger=stagger)
    SeqIO.write(out_recs, output, 'fastq')

def merged_reads(for_recs, rev_recs, max_diffs, min_size, max_size, stagger=False):
    for for_rec, rev_rec in zip(for_recs, rev_recs):
        rc_rev_rec = rev_rec.reverse_complement()
        best_diffs, best_overlap = best_merge_position(for_rec, rc_rev_rec, min_size, max_size, stagger=stagger)

        if best_diffs <= max_diffs:
            yield merged_at(for_rec, rc_rev_rec, best_overlap, stagger=stagger)

def best_merge_position(for_rec, rev_rec, min_size, max_size, stagger=False):
    read_lens = len(for_rec) + len(rev_rec)
    window = range(read_lens - max_size, read_lens - min_size)

    if not stagger:
        results = [(hamming_distance(for_rec[-overlap: ], rev_rec[0: overlap]), overlap) for overlap in window]
    else:
        results = [(hamming_distance(rev_rec[-overlap: ], for_rec[0: overlap]), overlap) for overlap in window]
        
    best_diffs, best_overlap = min(results)
    return best_diffs, best_overlap

def hamming_distance(x, y):
    return sum(map(operator.ne, x, y))

def merged_at(x, y, overlap, stagger=False):
    if not stagger:
        for_part = x[0: -overlap]
        rev_part = y[overlap: ]
        mid_part = merge_consensus(x[-overlap: ], y[0: overlap])
    else:
        for_part = x[overlap: ]
        rev_part = y[0: -overlap]
        mid_part = merge_consensus(x[0: overlap], y[-overlap: ])

    new_record = for_part + mid_part + rev_part
    new_record.id = for_part.id
    new_record.description = ''
    return new_record

def q_to_prob(q):
    return 10.0 ** (-float(q) / 10.0)

def same_base_prob(q1, q2):
    p1 = q_to_prob(q1)
    p2 = q_to_prob(q2)
    new_p = (p1 * p2 / 3.0) / (1.0 - p1 - p2 + 4.0 / 3.0 * p1 * p2)
    return prob_to_q(new_p)

def diff_base_prob(q1, q2):
    assert(q1 >= q2)
    p1 = q_to_prob(q1)
    p2 = q_to_prob(q2)
    new_p = p1 * (1.0 - p2 / 3.0) / (p1 + p2 - 4.0 / 3.0 * p1 * p2)
    return prob_to_q(new_p)

def prob_to_q(p):
    return min([41, math.floor(-10.0 * math.log10(p))])

def consensus_one_base_quality(b1, b2, q1, q2):
    if b1 == b2:
        return (b1, same_base_prob(q1, q2))
    else:
        if q1 >= q2:
            return(b1, diff_base_prob(q1, q2))
        else:
            return(b2, diff_base_prob(q2, q1))

def consensus_base_quality(seq1, seq2, quals1, quals2):
    each_base = zip(seq1, seq2, quals1, quals2)
    each_merged_base = [consensus_one_base_quality(*base) for base in each_base]
    new_seq, new_qual = zip(*each_merged_base)
    return "".join(new_seq), list(new_qual)

def merge_consensus(x, y):
    new_seq, new_qual = consensus_base_quality(str(x.seq), str(y.seq), x.letter_annotations['phred_quality'], y.letter_annotations['phred_quality'])
    new_record = SeqRecord.SeqRecord(seq=Seq.Seq(new_seq), letter_annotations={'phred_quality': new_qual}, description='')
    return new_record
