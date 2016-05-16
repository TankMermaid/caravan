#!/usr/bin/env python

'''
unit tests for merge.py
'''

import pytest
from caravan import primers
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

class TestSecondTrimmedRecord:
    def test_correct_trim(self):
        s =      "ACCTACGCTCCCTTTACACCCAGTAAATCCGGATAACGCTTGCCCCCTACGTATTACCGCGGCGGCTGGCACTATAAGA"
        trim_s = "ACCTACGCTCCCTTTACACCCAGTAAATCCGGATAACGCTTGCCCCCTACGTA"
        primer = "TTACCGCGGCKGCTGGCAC"
        rec = SeqRecord(Seq(s, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * len(s)})
        tr = primers.SecondTrimmedRecords(primer, iter([rec]), window=20, max_diffs=2)
        out = next(tr)

        assert str(out.seq) == trim_s
        with pytest.raises(StopIteration):
            next(tr)

    def test_correct_keep(self):
        s =      "ACCTACGCTCCCTTTACACCCAGTAAATCCGGATAACGCTTGCCCCCTACGTATTACCGCGGCGGCTGGCACTATAAGA"
        primer = "AAAAAAAAAAAAAAAA"
        rec = SeqRecord(Seq(s, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * len(s)})
        tr = primers.SecondTrimmedRecords(primer, iter([rec]), window=20, max_diffs=2)
        out = next(tr)

        assert out == rec
        with pytest.raises(StopIteration):
            next(tr)
