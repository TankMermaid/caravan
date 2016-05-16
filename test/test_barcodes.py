#!/usr/bin/env python

'''
unit tests for barcodes.py
'''

import pytest
from caravan import barcodes
from Bio import SeqIO

@pytest.fixture
def mapper(tmpdir):
    barcode_map = {'AACCGGTT': 'sample1'}
    fastq = tmpdir.join("for.fq")
    fastq.write('@header#ACGT/1\nAACCGGTT\n+\naaaaaaaa\n')
    fastq_fh = fastq.open()
    fastq_records = SeqIO.parse(fastq_fh, 'fastq')
    max_diffs = 0
    mr = barcodes.MappedRecords(barcode_map, fastq_records, max_diffs)
    return mr


class TestMapper:
    def test_correct(self, mapper):
        assert next(mapper) == 'header#ACGT/1\tsample1\n'


class TestHammingDistance:
    def test_correct_zero(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGT') == 0

    def test_correct_one(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGN') == 1

    def test_fail_with_different_lengths(self):
        with pytest.raises(RuntimeError):
            barcodes.MappedRecords.hamming_distance('ACGT', 'ACG')
