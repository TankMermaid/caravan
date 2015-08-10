#!/usr/bin/env python

'''
unit tests for barcodes.py
'''

import pytest
from caravan import barcodes
from Bio import SeqIO

@pytest.fixture
def mapper(tmpdir):
    barcode_map = {'ACGT': 'sample1'}
    fastq = tmpdir.join("for.fq")
    fastq.write('@OURSEQ:lolapalooza1234#ACGT/1\nAACCGGTT\n+\naaaaaaaa\n')
    fastq_fh = fastq.open()
    fastq_records = SeqIO.parse(fastq_fh, 'fastq')
    max_diffs = 0
    mr = barcodes.MappedRecords(barcode_map, fastq_records, max_diffs)
    return mr


class TestMapper:
    def test_correct(self, mapper):
        record = next(mapper)
        assert record.id == 'sample=sample1;1'


class TestParseBarcode:
    def test_correct(self):
        assert barcodes.MappedRecords.parse_barcode('@any set of|chars#ACGTN/1') == 'ACGTN'
    
    def test_fail_with_nonnt_chars(self):
        with pytest.raises(RuntimeError):
            barcodes.MappedRecords.parse_barcode('@foo#bar/1') 

    def test_fail_with_bad_directionality(self):
        with pytest.raises(RuntimeError):
            barcodes.MappedRecords.parse_barcode('@foo#ACGTN/3')


class TestHammingDistance:
    def test_correct_zero(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGT') == 0

    def test_correct_one(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGN') == 1

    def test_fail_with_different_lengths(self):
        with pytest.raises(AssertionError):
            barcodes.MappedRecords.hamming_distance('ACGT', 'ACG')