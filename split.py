#!/usr/bin/env python

'''
Split an input fasta or fastq file into k files.
'''

import itertools, os.path, sys, argparse
from math import ceil
from Bio import SeqIO
import util

class FastxSplitter():
    def __init__(self, fastx, size, filetype=None, mb=False, max_chunks=100):
        self.fastx = fastx
        self.size = self.size_parser(size)

        # guess filetype if not specified
        if filetype is None:
            self.filetype = self.infer_filetype(fastx)
        else:
            self.filetype = self.infer_filetype_from_alias(filetype)

        # how many output files will be needed?
        # size of file divided by chunk size, rounded up
        n_chunks = int(ceil(float(os.path.getsize(self.fastx)) / self.size))

        self.fns = self.output_filenames(self.fastx, n_chunks)
        util.check_for_collisions(self.fns)
        self.fhs = [open(fn, 'w') for fn in self.fns]
        self.split_entries(self.fastx, self.size, self.fhs, self.filetype)

    @staticmethod
    def infer_filetype_from_alias(x):
        '''recognize abbreviations for fasta and fastq'''
        if x in ['fasta', 'fa', 'fsa']:
            return 'fasta'
        elif x in ['fastq', 'fq']:
            return 'fastq'
        else:
            raise RuntimeError('alias {} not recognized as fasta or fastq'.format(x))

    @classmethod
    def infer_filetype(cls, fastx):
        '''infer fasta/q from the filename and test'''
        # first, make a guess based on the name
        # get the extension and scrub the dot
        extension = os.path.splitext(fastx)[-1].lstrip('.')
        filetype = cls.infer_filetype_from_alias(extension)

        # try parsing one entry
        SeqIO.parse(fastx, filetype).next()

        return filetype

    @staticmethod
    def size_parser(size_string):
        '''parse string [float][Mb/Gb] into bytes'''

        multiplier = size_string[-2:]
        if multiplier == "Gb":
            multiplier = 2**30
        elif multiplier == "Mb":
            multiplier = 2**20
        else:
            raise RuntimeError("multiplier %s not recognized" % multiplier)

        try:
            prefactor = float(size_string[:-2])
        except ValueError:
            raise RuntimeError("prefactor %s not recognized as float" % size_string[:-2])

        return prefactor * multiplier

    @staticmethod
    def output_filenames(input_filename, k):
        '''foo.fasta -> foo.fasta.0, foo.fasta.1, etc.'''
        return ['%s.%d' % (input_filename, i) for i in range(k)]

    @staticmethod
    def split_entries(inp, chunk, fhs, fastx, by_hash=False):
        '''
        Send entries in the input to filenames, cycling over each filename.
        
        inp : filehandle or filename
            input fastx file
        chunk : integer
            size of output files in bytes
        fhs : list or iterator of filehandles 
            output filehandles
        fastx : string
            'fasta' or 'fastq'
            
        returns : nothing
        '''

        fh = fhs.pop(0)
        for record in SeqIO.parse(inp, fastx):
            if fh.tell() > chunk:
                fh.close()
                fh = fhs.pop(0)

            SeqIO.write(record, fh, fastx)