'''
Demultiplex reads by mapping the barcode read to the sample names from a barcode mapping
file.

Given a tab-separated barcode mapping file like
    donor1_day5   ACGT

the first read mapping to that barcode, say
    @OURSEQ:lolapalooza1234#ACGT/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @sample=donor1_day5;1
    AACCGGTT
    +
    abcdefgh
    
where the ;1 means it's the first read that mapped to donor1_day5.
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

        sample = self.adhoc_barcode_map[self.parse_barcode(record.id)]
        if sample not in self.sample_counts:
            self.sample_counts[sample] = 1
        else:
            self.sample_counts[sample] += 1

        record.id = "sample={};{}".format(sample, self.sample_counts[sample])
        record.description = ''
        return record

    def recognized_record(self, record):
        barcode_read = self.parse_barcode(record.id)

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
    def parse_barcode(record_id):
        '''
        Extract the barcode read and direction from a BioPython SeqRecord
        
        Parameters
        record_id : str
            fastq SeqRecord id
        
        returns : string
            barcode read
        '''
        
        # match, e.g. @any_set_of_chars#ACGT/1 -> ACGT
        m = re.match(".*#([ACGTN]+)/[12]$", record_id)

        if m is None:
            raise RuntimeError("fastq id did not match expected format: %s" %(record_id))

        return m.group(1)

    @staticmethod
    def hamming_distance(x, y):
        assert len(x) == len(y)
        return sum(map(operator.ne, x, y))


class BarcodeMapper:
    def __init__(self, barcode_fasta, fastx, max_diffs, output, filetype, run=True):
        '''filetype is 'fasta' or 'fastq' '''

        self.barcode_fasta = SeqIO.parse(barcode_fasta, filetype)
        self.fastx_records = SeqIO.parse(fastx, filetype)
        self.max_diffs = max_diffs
        self.output = output
        self.filetype = filetype

        if run:
            self.run()

    def run(self):
        # parse the barcode mapping file
        self.barcode_map = {str(record.seq): record.id for record in self.barcode_fasta}

        # get a set of reads
        SeqIO.write(MappedRecords(self.barcode_map, self.fastx_records, self.max_diffs), self.output, self.filetype)