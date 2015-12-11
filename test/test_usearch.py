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

class TestDenovo:
    opts = {'radius': '3.0', 'fasta': 'foo.fa', 'output': 'otu.fa', 'index': 'otu.up'}
    expect_out = ['usearch', '-cluster_otus', 'foo.fa', '-otu_radius_pct', '3.0', '-otus', 'otu.fa', '-uparseout', 'otu.up']

    def test_correct(self):
        out = usearch.Usearcher(debug=True).cluster_denovo(**self.opts)
        assert out == self.expect_out

    def test_raise(self):
        new_opts = dict(self.opts)
        new_opts['radius'] = 97.0
        with pytest.raises(RuntimeError):
            usearch.Usearcher(debug=True).cluster_denovo(**new_opts)

    def test_override(self):
        new_opts = dict(self.opts)
        new_opts['radius'] = 97.0
        new_opts['force'] = True
        out = usearch.Usearcher(debug=True).cluster_denovo(**new_opts)
        assert out == ['usearch', '-cluster_otus', 'foo.fa', '-otu_radius_pct', '97.0', '-otus', 'otu.fa', '-uparseout', 'otu.up']

class TestSearch:
    opts = {'fasta': 'foo.fa', 'db': 'db.fa', 'sid': '0.9', 'b6': 'bar.b6'}
    expect_out = ['usearch', '-usearch_global', 'foo.fa', '-db', 'db.fa', '-id', '0.9', '-strand', 'both', '-blast6out', 'bar.b6', '-output_no_hits']

    def test_correct(self):
        out = usearch.Usearcher(debug=True).search(**self.opts)
        assert out == self.expect_out
