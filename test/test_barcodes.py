#!/usr/bin/env python
'''
unit tests for barcodes.py
'''

import pytest
from caravan import barcodes
from Bio import SeqIO


@pytest.fixture
def fastq_contents(request, tmpdir, scope='function'):
    'take fastq contents and make a temp file, parse with SeqIO, return iterator'
    fastq = tmpdir.join('for.fq')
    fastq.write(request.param)
    fastq_fh = fastq.open()
    fastq_records = SeqIO.parse(fastq_fh, 'fastq')
    return fastq_records


class TestMappedRecords:
    @pytest.mark.parametrize('barcode_map, fastq_contents, max_diffs, expected', [
        ({'AACCGGTT': 'sample1'}, '@header#ACGT/1\nAACCGGTT\n+\naaaaaaaa\n', 0, 'header#ACGT/1\tsample1\n'),
        ({'AACCGGTT': 'sample1'}, '@header#ACGT/1\nAACCGGGT\n+\naaaaaaaa\n', 1, 'header#ACGT/1\tsample1\n')
    ], indirect=['fastq_contents'])
    def test_correct(self, barcode_map, fastq_contents, max_diffs, expected):
        'barcode in fastq matches with barcode_map with specified differences'
        mapped_records = barcodes.MappedRecords(barcode_map, fastq_contents, max_diffs)
        assert next(mapped_records) == expected

    @pytest.mark.parametrize('barcode_map, fastq_contents, max_diffs', [
        ({'AACCGGTT': 'sample1'}, '@header#ACGT/1\nTTCCGGTT\n+\naaaaaaaa\n', 1),
        ({'GGAATTCC': 'sample1'}, '@header#ACGT/1\nAACCGGGT\n+\naaaaaaaa\n', 0)
    ], indirect=['fastq_contents'])
    def test_no_match(self, barcode_map, fastq_contents, max_diffs):
        'barcode in fastq does not match map with specified maxdiffs'
        with pytest.raises(StopIteration):
            mapped_records = barcodes.MappedRecords(barcode_map, fastq_contents, max_diffs)
            next(mapped_records)

    @pytest.mark.parametrize('barcode_map, fastq_contents, max_diffs', [
        ({'AACCGGTT': 'sample1'}, '@header#ACGT/1\nTTCCGGTT\n+\naaaaaaaa\n@header#ACGT/2\nTTCCGGTT\n+\naaaaaaaa\n', 1)
    ], indirect=['fastq_contents'])
    def test_bad_barcodes_cached(self, barcode_map, fastq_contents, max_diffs):
        'repeated bad barcode doesnt test distance again'
        with pytest.raises(StopIteration):
            mapped_records = barcodes.MappedRecords(barcode_map, fastq_contents, max_diffs)
            next(mapped_records)


class TestHammingDistance:
    def test_correct_zero(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGT') == 0

    def test_correct_one(self):
        assert barcodes.MappedRecords.hamming_distance('ACGT', 'ACGN') == 1

    def test_fail_with_different_lengths(self):
        with pytest.raises(RuntimeError):
            barcodes.MappedRecords.hamming_distance('ACGT', 'ACG')
