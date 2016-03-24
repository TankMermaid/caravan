#!/usr/bin/env python

'''
unit tests for merge.py
'''

import pytest
from caravan import merge
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

class TestBestMergePosition:
    def test_correct_small(self):
        for_seq = Seq("TCACAGTGCCAGCCGCCGCG", generic_dna)
        rev_seq =     Seq("AGTGCCAGCCGCCGCGTTTT", generic_dna)
        size = 24
        overlap = size - 8
        f = SeqRecord(for_seq, id="for")
        r = SeqRecord(rev_seq, id="rev")

        best_diffs, best_overlap = merge.best_merge_position(f, r, size - 1, size + 1)
        assert best_diffs == 0
        assert best_overlap == overlap

        f, r = r, f
        best_diffs, best_overlap = merge.best_merge_position(f, r, size - 1, size + 1, stagger=True)
        assert best_diffs == 0
        assert best_overlap == overlap

class TestMergedAt:
    def test_correct(self):
        for_seq = Seq("TCACAGTGCCAGCCGCCGCG", generic_dna)
        rev_seq =     Seq("AGTGCCAGCCGCCGCGTTTT", generic_dna)
        size = 24
        overlap = size - 8
        quals = {'phred_quality': [30] * 20}
        f = SeqRecord(for_seq, id="for", letter_annotations=quals)
        r = SeqRecord(rev_seq, id="rev", letter_annotations=quals)
        assert str(merge.merged_at(f, r, overlap, stagger=False).seq) == "TCACAGTGCCAGCCGCCGCGTTTT"

    def test_correct_merged(self):
        for_seq = Seq("TCACAGTGCCAGCCGCCGCG", generic_dna)
        rev_seq =     Seq("AGTGCCAGCCGCCGCGTTTT", generic_dna)
        size = 24
        overlap = size - 8
        quals = {'phred_quality': [30] * 20}
        f = SeqRecord(for_seq, id="for", letter_annotations=quals)
        r = SeqRecord(rev_seq, id="rev", letter_annotations=quals).reverse_complement()

        reads = merge.merged_reads([f], [r], max_diffs=1, min_size=size-1, max_size=size+1, stagger=False)
        assert str(next(reads).seq) == "TCACAGTGCCAGCCGCCGCGTTTT"
