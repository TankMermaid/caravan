'''
Demultiplex reads by mapping the index (aka "barcode") reads to the sample
names. Sample names and barcodes are given in a from a barcode mapping
file.

Given a barcode fasta file with entries like
    >sample1
    AAAA
    >sample2
    TTTT

a read mapping to that barcode, say
    @foo
    AAAA
    +
    JJJJ

becomes a line in the output like
    foo [tab] sample1
'''

import re, operator, itertools
from Bio import SeqIO

class MappedRecords():
    def __init__(self, barcode_map, fastx_records, max_diffs):
        self.original_barcode_map = barcode_map
        self.fastx_records = fastx_records
        self.max_diffs = max_diffs

        self.bad_barcodes = {}
        self.adhoc_barcode_map = dict(barcode_map)
        self.sample_counts = {}

    def __iter__(self):
        return self

    def __next__(self):
        record = next(self.fastx_records)
        while not self.recognized_record(record):
            record = next(self.fastx_records)

        sample = self.adhoc_barcode_map[str(record.seq)]
        if sample not in self.sample_counts:
            self.sample_counts[sample] = 1
        else:
            self.sample_counts[sample] += 1

        return '{}\t{}\n'.format(record.id, sample)

    def recognized_record(self, record):
        '''does this record have a barcode in our barcode map?'''
        barcode_read = str(record.seq)

        if barcode_read in self.adhoc_barcode_map:
            return True
        elif barcode_read in self.bad_barcodes:
            return False
        else:
            # this is a new barcode. need to compare it against existing.
            for known_code in self.original_barcode_map:
                if self.hamming_distance(barcode_read, known_code) <= self.max_diffs:
                    # this is a good enough match. in our ad hoc list, point this
                    # barcode to the same sample as the matching known barcode
                    self.adhoc_barcode_map[barcode_read] = self.original_barcode_map[known_code]
                    return True

            # this code isn't similar enough to any known code. throw it away.
            self.bad_barcodes[barcode_read] = 0
            return False

    @staticmethod
    def hamming_distance(x, y):
        if len(x) != len(y):
            raise RuntimeError("tried to compare barcodes of two lengths: {} and {}".format(x, y))

        return sum(map(operator.ne, x, y))


class BarcodeMapper:
    def __init__(self, barcode_fasta, fastx, max_diffs, output, input_format, run=True):
        '''input_format is 'fasta' or 'fastq' '''

        self.barcode_fasta = SeqIO.parse(barcode_fasta, 'fasta')
        self.fastx_records = SeqIO.parse(fastx, input_format)
        self.max_diffs = max_diffs
        self.output = output
        self.input_format = input_format

        if run:
            self.run()

    def run(self):
        # parse the barcode mapping file
        self.barcode_map = {str(record.seq): record.id for record in self.barcode_fasta}

        # get a set of reads
        for pair in MappedRecords(self.barcode_map, self.fastx_records, self.max_diffs):
            self.output.write(pair)
