#!/usr/bin/env python

'''
unit tests for usearch.py
'''

from caravan import usearch
import pytest

class TestMerge:
    opts = {'forward': 'for.fq', 'reverse': 'rev.fq', 'output': 'merge.fq'}
    expect_out = ['usearch', '-fastq_mergepairs', 'for.fq', '-reverse', 'rev.fq', '-fastqout', 'merge.fq']

    def test_correct(self):    
        out = usearch.Usearcher(debug=True).merge(**self.opts)
        assert out == self.expect_out

    def test_with_sizes(self):
        self.opts.update({'size': 100, 'size_var': 5})
        self.expect_out += ['-fastq_minmergelen', '95', '-fastq_maxmergelen', '105']
        out = usearch.Usearcher(debug=True).merge(**self.opts)
        assert out == self.expect_out