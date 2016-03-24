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
        whole_seq = Seq("CACAGTGCCAGCCGCCGCG", generic_dna)
        l = len(whole_seq)
        offset = 3
        overlap = l - offset
        for_seq = whole_seq[0: overlap]
        rev_seq = whole_seq[offset: ].reverse_complement()
        f = SeqRecord(for_seq, id="for")
        r = SeqRecord(rev_seq, id="rev")
        best_diffs, best_overlap = merge.best_merge_position(f, r, overlap - 1, overlap + 1)
        assert best_overlap == overlap
