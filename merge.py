'''
Merge pairs of reads
'''

import operator, math
from Bio import SeqIO, SeqRecord, Seq

def merge(forward, reverse, max_diffs, size, size_var, output):
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

    out_recs = merged_reads(for_recs, rev_recs, max_diffs, min_size, max_size)
    SeqIO.write(out_recs, output, 'fastq')

def merged_reads(for_recs, rev_recs, max_diffs, min_size, max_size):
    for for_rec, rev_rec in zip(for_recs, rev_recs):
        rc_rev_rec = rev_rec.reverse_complement()
        best_diffs, best_overlap = best_merge_position(for_rec, rc_rev_rec, min_size, max_size)

        if best_diffs <= max_diffs:
            yield merged_at(for_rec, rc_rev_rec, best_overlap)

def best_merge_position(for_rec, rev_rec, min_size, max_size):
    read_lens = len(for_rec) + len(rev_rec)
    window = range(read_lens - max_size, read_lens - min_size)
    results = [(hamming_distance(for_rec[-overlap: ], rev_rec[0: overlap]), overlap) for overlap in window]
    best_diffs, best_overlap = min(results)
    return best_diffs, best_overlap

def hamming_distance(x, y):
    return sum(map(operator.ne, x, y))

def merged_at(x, y, overlap):
    for_part = x[0: -overlap]
    rev_part = y[overlap: ]
    mid_part = merge_consensus(x[-overlap: ], y[0: overlap])
    new_record = for_part + mid_part + rev_part
    new_record.id = for_part.id
    new_record.description = ''
    return new_record

def merge_consensus(x, y):
    new_seq = ""
    new_qual = []
    for b1, b2, q1, q2 in zip(x.seq, y.seq, x.letter_annotations['phred_quality'], y.letter_annotations['phred_quality']):
        p1 = 10.0 ** (-q1 / 10.0)
        p2 = 10.0 ** (-q2 / 10.0)

        if b1 == b2:
            new_p = (p1 * p2 / 3.0) / (1.0 - p1 - p2 + 4.0 / 3.0 * p1 * p2)
        else:
            # swap if they are not "in order"
            if q1 < q2:
                b1, b2 = b2, b1
                p1, p2 = p2, p1
            
            new_p = p1 * (1.0 - p2 / 3.0) / (p1 + p2 - 4.0 / 3.0 * p1 * p2)

        new_seq += b1
        new_q = min([41, math.floor(-10.0 * math.log10(new_p))])
        new_qual.append(new_q)

    new_record = SeqRecord.SeqRecord(seq=Seq.Seq(new_seq), letter_annotations={'phred_quality': new_qual}, description='')
    return new_record