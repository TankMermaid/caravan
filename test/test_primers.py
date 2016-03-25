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
    s =      "ACCTACGCTCCCTTTACACCCAGTAAATCCGGATAACGCTTGCCCCCTACGTATTACCGCGGCGGCTGGCACTATAAGA"
    trim_s = "ACCTACGCTCCCTTTACACCCAGTAAATCCGGATAACGCTTGCCCCCTACGTA"
    rec = SeqRecord(Seq(s, generic_dna), id="foo", letter_annotations={'phred_quality': [30] * len(s)})
    tr = primers.SecondTrimmedRecords("TTACCGCGGCKGCTGGCAC", iter([rec]), window=20, max_diffs=2)
    out = next(tr)

    assert str(out.seq) == trim_s
