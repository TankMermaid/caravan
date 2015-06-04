#!/usr/bin/env python

'''
Make a new fasta file with only the unique sequences from the original. The
new sequences are in abundance order.

If computing an index file, the sample name either is what follows "sample="
and before a semicolon or is the first part of the fasta label before a
semicolon.

Sequences that do not reach a minimum number of counts can be dropped.
'''

import sys, argparse, re, sys, json, itertools
from operator import itemgetter
from Bio import SeqIO, SeqRecord, Seq
import util

class Dereplicator():
    def __init__(self, fasta, min_counts, output, index=None, run=True):
        '''
        fasta : fasta filename or handle
            input fasta
        min_counts : int
            minimum number of counts to be included in output
        output : filehandle or filename
            output dereplicated fasta
        index : filehandle
            output index file
        '''
        
        self.fasta = fasta
        self.min_counts = min_counts
        self.output = output
        self.index = index

        if run:
            self.run()
        
    def run(self):
        self.counts, self.provenances = self.dereplicate(self.fasta)

        self.write_output(self.counts, self.output, self.provenances, self.index, self.min_counts)
        
    @staticmethod
    def dereplicate(fasta):
        '''
        Process the input fasta entries.

        Populate a dictionary
        {seq => {sample => counts}}. Use the special sample key None to
        indicate the sum of all counts for that sequence.
        '''
        
        counts = {}
        provenances = {}
        for record in SeqIO.parse(fasta, 'fasta'):
            seq = str(record.seq)
            sample = util.strip_fasta_label(record.id)

            if seq not in counts:
                counts[seq] = 1
            else:
                counts[seq] += 1

            if seq not in provenances:
                provenances[seq] = {sample: 1}
            elif sample not in provenances[seq]:
                provenances[seq][sample] = 1
            else:
                provenances[seq][sample] += 1

        return counts, provenances

    @classmethod
    def write_output(cls, counts, output, provenances=None, index=None, min_counts=0):
        '''
        Write the new fasta file in abundance order. If the provenances and
        output index file are specified, write the index file too.
        '''

        # get the seqs in abundance order
        seqs_sizes = sorted(counts.items(), key=itemgetter(1), reverse=True)

        SeqIO.write(cls.new_fasta_entries(seqs_sizes, min_counts), output, 'fasta')

        if provenances is not None and index is not None:
            # remove seqs below the count threshold and rename the provenances to be in
            # the ordered seq1, seq2, etc.
            for i, [seq, size] in enumerate(seqs_sizes):
                if size < min_counts:
                    del provenances[seq]
                else:
                    provenances["seq%d" %(i+1)] = provenances.pop(seq)

            json.dump(provenances, index)

    @classmethod
    def new_fasta_entries(cls, seqs_sizes, min_counts):
        '''yield the list in abundance order'''

        for i, [seq, size] in enumerate(seqs_sizes):
            if size >= min_counts:
                yield SeqRecord.SeqRecord(Seq.Seq(seq), id="seq%d;size=%d;" %(i+1, size), description='')