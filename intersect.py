'''
Intersect a mapping file
'''

from Bio import SeqIO

class Intersect:
    def __init__(self, mapping_lines, fastq_entries):
        self.mapping_lines = mapping_lines
        self.fastq_entries = fastq_entries

        self.sample_counts = {}

    def __iter__(self):
        return self

    def get_next_map(self):
        self.map_id, self.sample = next(self.mapping_lines).split("\t")
        self.map_num = self.extract_number(self.map_id)

    def get_next_entry(self):
        self.entry = next(self.fastq_entries)
        self.entry_num = self.extract_number(self.entry.id)

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

    @staticmethod
    def extract_number(header):
        return int(header.lstrip('read'))

class Intersecter:
    def __init__(self, mapping, fastq, output):
        self.mapping = mapping
        self.fastq = fastq
        self.output = output

        for entry in Intersect(self.mapping, SeqIO.parse(self.fastq, 'fastq')):
            SeqIO.write(entry, self.output, 'fastq')