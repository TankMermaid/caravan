#!/usr/bin/env python

'''
unit tests for barcodes.py
'''

from caravan import barcodes
import pytest

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