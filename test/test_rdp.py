#!/usr/bin/env python

'''
unit tests for rdp
'''

from caravan import rdp
import pytest, io, json

class WithRank:
    rank = rdp.FixrankRank('domain', 'Bacteria', '0.8')

class TestRankRepr(WithRank):
    def test_correct(self):
        rep = 'FixrankRank("domain", "Bacteria", 0.8)'
        assert repr(self.rank) == rep

class TestRankStr(WithRank):
    def test_correct(self):
        assert str(self.rank) == 'Bacteria'

class TestRankEqual:
    def test_correct(self):
        rank1 = rdp.FixrankRank('domain', 'Bacteria', '0.8')
        rank2 = rdp.FixrankRank('domain', 'Bacteria', '0.8')
        assert rank1 == rank2

class TestLineageEqual:
    def test_correct(self):
        lin1 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'P', 1.0]])
        lin2 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'P', 1.0]])
        assert lin1 == lin2

    def test_different_info(self):
        lin1 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'P', 1.0]])
        lin2 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'something_else', 1.0]])
        assert lin1 != lin2

    def test_different_lengths(self):
        lin1 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'P', 1.0]])
        lin2 = rdp.FixrankLineage([['domain', 'Bacteria', 1.0]])
        assert lin1 != lin2

class WithLineage:
    lin = rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'P', 1.0]])

class TestLineageRepr(WithLineage):
    def test_correct(self):    
        rep = 'FixrankLineage([["domain", "Bacteria", 1.0], ["phylum", "P", 1.0]])'
        assert repr(self.lin) == rep

class TestLineageStr(WithLineage):
    def test_correct(self):
        assert str(self.lin) == 'Bacteria;P'

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
        sid, lin = rdp.FixrankParser.parse_entries(['seq1', '', 'Root', 'rootrank', '1.0', 'K', 'domain', '1.0', 'P', 'phylum', '0.5', 'C', 'class', '1.0'])
        assert lin.ranks_to_rank('phylum') == [rdp.FixrankRank('domain', 'K', 1.0), rdp.FixrankRank('phylum', 'P', 0.5)]

class TestTrimToRank:
    def test_correct(self):
        lin = rdp.FixrankLineage([['domain', 'D', 1.0], ['phylum', 'P', 1.0], ['class', 'C', 0.75], ['order', 'O', 1.0], ['family', 'F', 1.0], ['genus', 'G', 1.0]])
        lin.trim_at_rank('phylum')
        assert lin == rdp.FixrankLineage([['domain', 'D', 1.0], ['phylum', 'P', 1.0]])

class TestParseLine:
    def test_correct(self):
        line = 'seq494517;size=2;\t\tRoot\trootrank\t1.0\tBacteria\tdomain\t1.0\t"Proteobacteria"\tphylum\t1.0\tAlphaproteobacteria\tclass\t0.89\tCaulobacterales\torder\t0.42\tCaulobacteraceae\tfamily\t0.41\tBrevundimonas\tgenus\t0.38\n'
        sid, lin = rdp.FixrankParser.parse_line(line)
        assert sid == 'seq494517'
        assert lin == rdp.FixrankLineage([['domain', 'Bacteria', 1.0], ['phylum', 'Proteobacteria', 1.0], ['class', 'Alphaproteobacteria', 0.89], ['order', 'Caulobacterales', 0.42], ['family', 'Caulobacteraceae', 0.41], ['genus', 'Brevundimonas', 0.38]])

class WithContent:
    content = '''seq1;size=272905;\t\tRoot\trootrank\t1.0\tBacteria\tdomain\t1.0\t"Bacteroidetes"\tphylum\t0.0\tFlavobacteriia\tclass\t1.0\t"Flavobacteriales"\torder\t1.0\tFlavobacteriaceae\tfamily\t1.0\tFlavobacterium\tgenus\t1.0
seq2;size=229776;\t\tRoot\trootrank\t1.0\tBacteria\tdomain\t1.0\t"Bacteroidetes"\tphylum\t1.0\tFlavobacteriia\tclass\t0.0\t"Flavobacteriales"\torder\t1.0\tFlavobacteriaceae\tfamily\t1.0\tFlavobacterium\tgenus\t1.0
seq3;size=212890;\t\tRoot\trootrank\t1.0\tBacteria\tdomain\t1.0\t"Bacteroidetes"\tphylum\t1.0\tFlavobacteriia\tclass\t1.0\t"Flavobacteriales"\torder\t0.0\tFlavobacteriaceae\tfamily\t1.0\tFlavobacterium\tgenus\t1.0'''

class TestParseLines(WithContent):
    def test_correct(self):
        lines = self.content.split("\n")
        expected = {'seq1': 'Bacteria', 'seq2': 'Bacteria;Bacteroidetes', 'seq3': 'Bacteria;Bacteroidetes;Flavobacteriia'}
        mapping = rdp.FixrankParser.parse_lines(lines, 0.5)
        assert expected == mapping

class TestParseFile(WithContent):
    def test_correct(self):
        output_fh = io.StringIO()
        fixrank = io.StringIO(self.content)
        rdp.FixrankParser.parse_file(fixrank, 'p', output_fh, 0.5)

        output = json.loads(output_fh.getvalue())
        expected = {'seq1': 'Bacteria', 'seq2': 'Bacteria;Bacteroidetes', 'seq3': 'Bacteria;Bacteroidetes'}
        assert output == expected