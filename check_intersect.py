'''
Compare the forward and reverse fastq files. Make sure that there are pairs of forward and
reverse reads that are appropriately labeled.
'''

import itertools, subprocess
from Bio import SeqIO

class IntersectChecker():
    def __init__(self, forward, reverse):
        self.forward = forward
        self.reverse = reverse

        self.run()

    def run(self):
        self.check_for_matched_line_numbers(self.forward, self.reverse)
        self.check_for_paired_ids(self.forward, self.reverse)

    @staticmethod
    def n_lines(fn):
        '''count number of lines in fn using wc'''
        return int(subprocess.check_output(['wc', '-l', fn]).split()[0])

    @classmethod
    def check_for_matched_line_numbers(cls, f1, f2):
        '''assert that two files have the same number of lines'''
        if cls.n_lines(f1) != cls.n_lines(f2):
            raise RuntimeError("fastqs have unequal numbers of lines")

    @staticmethod
    def check_matching_fastq_ids(forward_id, reverse_id):
        '''check that reads are whatever/1 and whatever/2'''
        if not forward_id.endswith("/1"):
            raise RuntimeError("fastq id not correct direction: '%s' should end in 1" % forward_id)

        if not reverse_id.endswith("/2"):
            raise RuntimeError("fastq id not correct direction: '%s' should end in 2" % reverse_id)

        if forward_id.rstrip("/12") != reverse_id.rstrip("/12"):
            raise RuntimeError("fastq ids did not match: '%s' and '%s'" % (forward_id, reverse_id))

    @classmethod
    def check_for_paired_ids(cls, forward_fastq, reverse_fastq):
        '''check each pair of fastq records for paired read ids'''
        for forward_record, reverse_record in itertools.izip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq, 'fastq')):
            cls.check_matching_fastq_ids(forward_record.id, reverse_record.id)