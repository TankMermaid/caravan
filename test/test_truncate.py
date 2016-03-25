#!/usr/bin/env python

'''
unit tests for truncate.py
'''

import pytest
from caravan import truncate
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

class TestLengthedEntries:
    def test_correct_long(self):
        l = 150
        out_length = 100
        rec = SeqRecord(Seq("A" * l, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * l})
        tr = truncate.lengthed_entries(iter([rec]), out_length, keep=False)
        out = next(tr)

        assert len(out) == out_length

    def test_correct_short_keep(self):
        l = 150
        out_length = 200
        rec = SeqRecord(Seq("A" * l, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * l})
        tr = truncate.lengthed_entries(iter([rec]), out_length, keep=True)
        out = next(tr)

        assert len(out) == l

    def test_correct_short_drop(self):
        l = 150
        out_length = 200
        rec = SeqRecord(Seq("A" * l, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * l})
        tr = truncate.lengthed_entries(iter([rec]), out_length, keep=False)

        with pytest.raises(StopIteration):
            out = next(tr)


class TestTailedEntries:
    def test_correct_trim(self):
        good_len = 100
        bad_len = 50
        total_len = good_len + bad_len
        rec = SeqRecord(Seq("A" * total_len, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * good_len + [2] * bad_len})
        tr = truncate.tailed_entries(iter([rec]), 2)

        out = next(tr)
        assert len(out) == good_len

    def test_correct_drop(self):
        total_len = 100
        rec = SeqRecord(Seq("A" * total_len, generic_dna), id="foo", letter_annotations={'phred_quality': [2] * total_len})
        tr = truncate.tailed_entries(iter([rec]), 2)

        with pytest.raises(StopIteration):
            out = next(tr)
