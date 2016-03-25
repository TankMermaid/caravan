#!/usr/bin/env python

'''
unit tests for intersect.py
'''

import pytest, tempfile, io
from caravan import intersect

def fake_fastq_lines(*headers):
    out = []
    for header in headers:
        out += ["@" + header, "AAA", "+", "###"]
    return "\n".join(out) + "\n"

def fake_fasta_lines(*headers):
    out = []
    for header in headers:
        out += [">" + header, "AAA"]
    return "\n".join(out) + "\n"

@pytest.fixture
def fastq1():
    fh = tempfile.TemporaryFile(mode='w+')
    fh.write(fake_fastq_lines("read1", "read3", "read5"))
    fh.seek(0)
    return fh

@pytest.fixture
def fasta1():
    fh = tempfile.TemporaryFile(mode='w+')
    fh.write(fake_fasta_lines("read1", "read2", "read3"))
    fh.seek(0)
    return fh

class TestExtractNumber:
    def test_correct(self):
        assert intersect.extract_number("read1234") == 1234

class TestPeek:
    def test_correct_fq(self, fastq1):
        res = intersect.peek(fastq1, None)
        assert isinstance(res, intersect.FastqRecords)

    def test_correct_fa(self, fasta1):
        res = intersect.peek(fasta1, None)
        assert isinstance(res, intersect.FastaRecords)

    def test_correct_fq_fa(self, fasta1, fastq1):
        ofa = io.StringIO()
        ofq = io.StringIO()
        intersect.intersect([fasta1, fastq1], [ofa, ofq])
        assert ofa.getvalue() == fake_fasta_lines("read1", "read3")
        assert ofq.getvalue() == fake_fastq_lines("read1", "read3")
