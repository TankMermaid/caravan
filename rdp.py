'''
parse fixrank files output by classifier.jar
'''

import csv, re, itertools

class FixrankLineage:
    @staticmethod
    def parse_triplet(triplet):
        '''replace quotes, cast as float, and reorder'''
        assert(len(triplet) == 3)
        taxon = re.gsub('"', '', triplet[0])
        level = triplet[1]
        confidence = float(triplet[2])

        return level, taxon, confidence

    @staticmethod
    def parse_sid_entry(entry):
        return entry.split(';')[0]

    @classmethod
    def parse_entries(cls, entries):
        '''break up entries into sid entry and triplet'''
        # first entry is sequence id
        sid_entry = cls.parse_sid_entry(entries[0])

        # second entry should be blank
        assert entries[1] == ""

        # the rest of the entries should divide into triplets
        assert (len(entries) - 2) % 3 == 0
        entry_triplets = zip(*[iter(entries[2:])] * 3)
        triplets = [cls.parse_triplet(t) for t in entry_triplets]

        # the first entry should be root, whatever that is
        assert triplets[0] == ['rootrank', 'Root', 1.0]

        return sid_entry, triplets[1:]

    def __init__(self, entries):
        self.sid, self.triples = self.parse_entries(entries)

    def trim_to_confidence(self, confidence):
        new_triplets = list(itertools.takewhile(lambda triplet: triplet[2] > confidence, self.triplets))
        self.triplets = new_triplets

class FixrankParser:
    pass
    