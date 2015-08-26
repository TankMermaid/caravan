#!/usr/bin/env python3

import argparse, sys, os, gzip
from Bio import SeqIO

def shared_id(rid):
    return rid.split()[0]

def new_id(location, barcode, direction):
    '''(location, barcode, direction) -> location#barcode/direction'''
    return "%s#%s/%s" %(location, barcode, direction)

if __name__ == '__main__':
    # parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('forward_in', help='input forward reads fastq')
    parser.add_argument('reverse_in', help='input reverse reads fastq')
    parser.add_argument('index', help='input index reads fastq')
    parser.add_argument('forward_out', help='output forward reads fastq')
    parser.add_argument('reverse_out', help='output reverse reads fastq')
    parser.add_argument('--gzip', '-g', action='store_true', help='input files are gzipped?')
    args = parser.parse_args()

    # open up the filehandles, zipped or not
    if args.gzip:
        opener = lambda x: gzip.open(x, 'rt') # open a gzip in text format
    else:
        opener = lambda x: open(x)

    fih, rih, iih = [opener(fn) for fn in [args.forward_in, args.reverse_in, args.index]]
    
    with open(args.forward_out, 'w') as fo, open(args.reverse_out, 'w') as ro:
        fi = SeqIO.parse(fih, 'fastq')
        ri = SeqIO.parse(rih, 'fastq')
        ii = SeqIO.parse(iih, 'fastq')
        
        for fr, rr, ir in zip(fi, ri, ii):
            # check that the reads match
            if not (shared_id(fr.id) == shared_id(rr.id) == shared_id(ir.id)):
                raise RuntimeError("ids in corresponding entries did not match: %s %s %s" %(fr.id, rr.id, ir.id))
            
            # the barcode is the sequence from the index read
            barcode = str(ir.seq)
            
            # amend the ids on the forward and reverse entries
            fr.id = new_id(fr.id, barcode, '1')
            rr.id = new_id(rr.id, barcode, '2')
            fr.description = ''
            rr.description = ''
            
            # write the entries to their output files
            SeqIO.write(fr, fo, 'fastq')
            SeqIO.write(rr, ro, 'fastq')
