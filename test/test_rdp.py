#!/usr/bin/env python

'''
unit tests for rdp
'''

from caravan import rdp
import pytest

class TestRankRepr:
    def test_correct(self):
        rank = rdp.FixrankRank('domain', 'Bacteria', '0.8')
        rep = 'FixrankRank("domain", "Bacteria", 0.8)'
        assert repr(rank) == rep

class TestRankEqual:
    def test_correct(self):
        rank1 = rdp.FixrankRank('domain', 'Bacteria', '0.8')
        rank2 = rdp.FixrankRank('domain', 'Bacteria', '0.8')
        assert rank1 == rank2

class TestParseTriplet:
    def test_correct(self):
        rank = rdp.FixrankParser.parse_triplet(['"Bacteroidetes"', 'phylum', '1.0'])
        assert rank.name == 'phylum'
        assert rank.taxon == 'Bacteroidetes'
        assert rank.confidence == 1.0

    def test_raise(self):
        with pytest.raises(Exception):
            rdp.FixrankParser.parse_triplet(['FakePhylum', 'phylum', 'garbage'])

class TestSidEntry:
    def test_correct(self):
        assert rdp.FixrankParser.parse_sid_entry('seq4;size=500') == 'seq4'

class WithEntries:
    entries = ['seq1;size=272905;', '', 'Root', 'rootrank', '1.0', 'Bacteria', 'domain', '1.0', '"Bacteroidetes"', 'phylum', '1.0', 'Flavobacteriia', 'class', '0.75', '"Flavobacteriales"', 'order', '1.0', 'Flavobacteriaceae', 'family', '1.0', 'Flavobacterium', 'genus', '1.0']
    sid = 'seq1'
    lin = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'Bacteroidetes', 1.0], ['class', 'Flavobacteriia', 0.75], ['order', 'Flavobacteriales', 1.0], ['family', 'Flavobacteriaceae', 1.0], ['genus', 'Flavobacterium', 1.0]])

class TestParseEntries(WithEntries):
    def test_correct(self):
        assert rdp.FixrankParser.parse_entries(self.entries) == (self.sid, self.lin)

class TestRanksAtConfidence(WithEntries):
    def test_correct(self):
        assert self.lin.ranks_at_confidence(0.8) == rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'Bacteroidetes', 1.0]]).ranks

class TestRanksToRank:
    def test_correct(self):
        sid, ranks = rdp.FixrankParser.parse_entries(['seq1', '', 'Root', 'rootrank', '1.0', 'K', 'domain', '1.0', 'P', 'phylum', '0.5', 'C', 'class', '1.0'])
        assert ranks.ranks_to_rank('phylum') == [rdp.FixrankRank('domain', 'K', 1.0), rdp.FixrankRank('phylum', 'P', 0.5)]