#!/usr/bin/env python

'''
unit tests for rdp
'''

from caravan import rdp
import pytest

class TestParseTriplet:
    def test_correct(self):
        assert rdp.FixrankLineage.parse_triplet(['"Bacteroidetes"', 'phylum', '1.0']) == ['phylum', 'Bacteroidetes', 1.0]

    def test_raise(self):
        with pytest.raises(Exception):
            rdp.FixrankLineage.parse_triplet(['FakePhylum', 'phylum', 'garbage'])

class TestSidEntry:
    def test_correct(self):
        assert rdp.FixrankLineage.parse_sid_entry('seq4;size=500') == 'seq4'

class TestParseEntries:
    def test_correct(self):
        entries = ['seq1;size=272905;', '', 'Root', 'rootrank', '1.0', 'Bacteria', 'domain', '1.0', '"Bacteroidetes"', 'phylum', '1.0', 'Flavobacteriia', 'class', '1.0', '"Flavobacteriales"', 'order', '1.0', 'Flavobacteriaceae', 'family', '1.0', 'Flavobacterium', 'genus', '1.0']
        result = 'seq1', [['domain', 'Bacteria', 1.0], ['phylum', 'Bacteroidetes', 1.0], ['class', 'Flavobacteriia', 1.0], ['order', 'Flavobacteriales', 1.0], ['family', 'Flavobacteriaceae', 1.0], ['genus', 'Flavobacterium', 1.0]]
        assert rdp.FixrankLineage.parse_entries(entries) == result

class TestTrimToConfidence:
    def test_correct(self):
        lin = rdp.FixrankLineage(['seq1', '', 'Root', 'rootrank', '1.0', 'K', 'domain', '1.0', 'P', 'phylum', '0.5', 'C', 'class', '1.0'])
        lin.trim_to_confidence(0.8)
        assert lin.triplets == [['domain', 'K', 1.0]]