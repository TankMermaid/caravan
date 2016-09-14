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
    opts = {'percent_identity': 97.5, 'fasta': 'foo.fa', 'output': 'otu.fa', 'index': 'otu.up'}
    expect_out = ['usearch', '-cluster_otus', 'foo.fa', '-otu_radius_pct', '2.5', '-otus', 'otu.fa', '-uparseout', 'otu.up']

    def test_correct(self):
        out = usearch.Usearcher(debug=True).cluster_denovo(**self.opts)
        assert out == self.expect_out

    def test_raise(self):
        new_opts = dict(self.opts)
        new_opts['percent_identity'] = 90
        with pytest.raises(RuntimeError):
            usearch.Usearcher(debug=True).cluster_denovo(**new_opts)

    def test_raise_negative(self):
        new_opts = dict(self.opts)
        new_opts['percent_identity'] = -1
        with pytest.raises(RuntimeError):
            usearch.Usearcher(debug=True).cluster_denovo(**new_opts)

    def test_raise_over100(self):
        new_opts = dict(self.opts)
        new_opts['percent_identity'] = 101
        with pytest.raises(RuntimeError):
            usearch.Usearcher(debug=True).cluster_denovo(**new_opts)

    def test_override(self):
        new_opts = dict(self.opts)
        new_opts['percent_identity'] = 90
        new_opts['force'] = True
        out = usearch.Usearcher(debug=True).cluster_denovo(**new_opts)
        assert out == ['usearch', '-cluster_otus', 'foo.fa', '-otu_radius_pct', '10.0', '-otus', 'otu.fa', '-uparseout', 'otu.up']

class TestSearch:
    opts = {'fasta': 'foo.fa', 'db': 'db.fa', 'percent_identity': 90.0, 'b6': 'bar.b6'}
    expect_out = ['usearch', '-usearch_global', 'foo.fa', '-db', 'db.fa', '-id', '0.9', '-strand', 'both', '-blast6out', 'bar.b6', '-output_no_hits']

    def test_correct(self):
        out = usearch.Usearcher(debug=True).search(**self.opts)
        assert out == self.expect_out
