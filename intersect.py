'''
Intersect a mapping file with one (or two) fastq's
'''

from Bio import SeqIO
import re

def extract_number(header):
    m = re.match("read(\d+);?", header)
    return int(m.groups()[0])

class Intersect2:
    def __init__(self, mapping_lines, fastq_entries):
        self.mapping_lines = mapping_lines
        self.fastq_entries = fastq_entries

        self.sample_counts = {}

    def __iter__(self):
        return self

    def get_next_map(self):
        self.map_id, self.sample = next(self.mapping_lines).split("\t")
        self.map_num = extract_number(self.map_id)

    def get_next_entry(self):
        self.entry = next(self.fastq_entries)
        self.entry_num = extract_number(self.entry.id)

    def __next__(self):
        self.get_next_map()
        self.get_next_entry()

        while self.map_num != self.entry_num:
            if self.map_num < self.entry_num:
                self.get_next_map()
            elif self.map_num > self.entry_num:
                self.get_next_entry()

        assert(self.map_num == self.entry_num)

        if self.sample in self.sample_counts:
            self.sample_counts[self.sample] += 1
        else:
            self.sample_counts[self.sample] = 1

        self.entry.id += ';sample={};{}'.format(self.sample, self.sample_counts[self.sample])
        self.entry.description = ''
        return self.entry


class Intersect3:
    def __init__(self, mapping, for_entries, rev_entries):
        self.mapping = mapping
        self.for_entries = for_entries
        self.rev_entries = rev_entries

        self.sample_counts = {}

    def __iter__(self):
        return self

    def get_next_map(self):
        self.map_id, self.sample = next(self.mapping).rstrip().split("\t")
        self.map_num = extract_number(self.map_id)

    def get_next_for_entry(self):
        self.for_entry = next(self.for_entries)
        self.for_entry_num = extract_number(self.for_entry.id)

    def get_next_rev_entry(self):
        self.rev_entry = next(self.rev_entries)
        self.rev_entry_num = extract_number(self.rev_entry.id)

    def nums(self):
        return [self.map_num, self.for_entry_num, self.rev_entry_num]

    def nums_equal(self):
        return self.map_num == self.for_entry_num == self.rev_entry_num

    def __next__(self):
        self.get_next_map()
        self.get_next_for_entry()
        self.get_next_rev_entry()

        while not self.nums_equal():
            if min(self.nums()) == self.nums()[0]:
                self.get_next_map()
            elif min(self.nums()) == self.nums()[1]:
                self.get_next_for_entry()
            elif min(self.nums()) == self.nums()[2]:
                self.get_next_rev_entry()

        assert(self.nums_equal())

        if self.sample in self.sample_counts:
            self.sample_counts[self.sample] += 1
        else:
            self.sample_counts[self.sample] = 1

        coda = ';sample={};{}'.format(self.sample, self.sample_counts[self.sample])
        self.for_entry.id += coda
        self.for_entry.description = ''
        self.rev_entry.id += coda
        self.rev_entry.description = ''
        
        return self.for_entry, self.rev_entry


def intersect2(mapping, forward, output):
    entries = SeqIO.parse(forward, 'fastq')
    new_entries = Intersect2(mapping, entries)
    SeqIO.write(new_entries, output, 'fastq')

def intersect3(mapping, forward, reverse, forward_output, reverse_output):
    for_entries = SeqIO.parse(forward, 'fastq')
    rev_entries = SeqIO.parse(reverse, 'fastq')
    new_entries = Intersect3(mapping, for_entries, rev_entries)
    for new_for, new_rev in new_entries:
        SeqIO.write(new_for, forward_output, 'fastq')
        SeqIO.write(new_rev, reverse_output, 'fastq')
