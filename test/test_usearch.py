#!/usr/bin/env python

'''
unit tests for usearch.py
'''

from caravan import usearch
import pytest

class TestMerge:
    def test_correct(self):
        opts = {'forward': 'for.fq', 'reverse': 'rev.fq', 'output': 'merge.fq'}
        out = usearch.Usearcher(debug=True).merge(**opts)
        assert out == ['usearch', '-fastq_mergepairs', 'for.fq', '-reverse', 'rev.fq', '-fastqout', 'merge.fq']