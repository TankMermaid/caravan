#!/usr/bin/env python
#
# author: scott olesen <swo@mit.edu>

import argparse, os, os.path, glob
from Bio import SeqIO

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("path")
    p.add_argument("forward", help="output forward fastq", type=argparse.FileType('w'))
    p.add_argument("reverse", help="output reverse fastq", type=argparse.FileType('w'))
    p.add_argument("map", help="output map fastq", type=argparse.FileType('w'))
    args = p.parse_args()

    # find the lists of files that will be used
    forwards = glob.glob("{}/*/*_1_sequence.fastq".format(args.path))
    reverses = glob.glob("{}/*/*_2_sequence.fastq".format(args.path))

    # get the sample names from that list
    sample_name = lambda p: os.path.split(os.path.split(p)[0])[1]
    samples = [sample_name(f) for f in forwards]
    assert samples == [sample_name(r) for r in reverses]

    entry_i = 1
    for sample, forward_fastq, reverse_fastq in zip(samples, forwards, reverses):
        for for_entry, rev_entry in zip(SeqIO.parse(forward_fastq, 'fastq'), SeqIO.parse(reverse_fastq, 'fastq')):
            read = "read{}".format(entry_i)
            entry_i += 1
            for_entry.id = read
            rev_entry.id = read
            for_entry.description = ""
            rev_entry.description = ""

            SeqIO.write(for_entry, args.forward, 'fastq')
            SeqIO.write(rev_entry, args.reverse, 'fastq')
            args.map.write(sample + "\n")
