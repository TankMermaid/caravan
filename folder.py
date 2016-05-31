#!/usr/bin/env python
#
# author: scott olesen <swo@mit.edu>

import argparse, os, re
from Bio import SeqIO

def matches_regexes(query, pos, neg):
    '''does query match positive regex and not match negative?'''
    return re.search(pos, query) is not None and re.search(neg, query) is None

def first_match(lst, regex, only=True):
    candidates = [elt for elt in lst if re.search(regex, elt)]

    if len(candidates) == 0:
        raise RuntimeError("could not find regex '{}' in list: {}".format(regex, lst))

    if only:
        assert len(candidates) == 1

    return candidates[0]

def fastq_paths_in_directory(target_dir, dir_pos_regex, dir_neg_regex='infosite', for_fq_regex='_1_.+\.fastq', rev_fq_regex='_2_.+\.fastq'):
    structure = {x: {'dirs': y, 'files': z} for x, y, z in os.walk(target_dir)}
    assert target_dir in structure

    top_dir_names = [d for d in structure[target_dir]['dirs'] if matches_regexes(d, dir_pos_regex, dir_neg_regex)]
    top_dir_paths = [os.path.join(target_dir, d) for d in top_dir_names]

    # check that all these top directories don't have subdirectories
    for d in top_dir_paths:
        assert structure[d]['dirs'] == []

    # find the forward and reverse reads for each 
    for_fq_names = [first_match(structure[d]['files'], for_fq_regex) for d in top_dir_paths]
    rev_fq_names = [first_match(structure[d]['files'], rev_fq_regex) for d in top_dir_paths]

    for_fq_paths = [os.path.join(d, fn) for d, fn in zip(top_dir_paths, for_fq_names)]
    rev_fq_paths = [os.path.join(d, fn) for d, fn in zip(top_dir_paths, rev_fq_names)]

    return (for_fq_paths, rev_fq_paths, top_dir_names)

def renamed_entries(fns, samples):
    for fn, sample in zip(fns, samples):
        for i, record in enumerate(SeqIO.parse(fn, 'fastq')):
            record.id = "sample={};{}".format(sample, i + 1)
            record.description = ""
            yield record

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="concatenate fastq's from multiple directories")
    p.add_argument('dir_pos_regex', help="to identify subfolders you want to keep")
    p.add_argument('for_out', type=argparse.FileType('w'))
    p.add_argument('rev_out', type=argparse.FileType('w'))
    p.add_argument('top_dirs', nargs='+', help="top directories from BMC")
    p.add_argument('--verbose', '-v', action='store_true')
    args = p.parse_args()

    for_fq_paths = []
    rev_fq_paths = []
    names = []
    for top_dir in args.top_dirs:
        new_for, new_rev, new_names = fastq_paths_in_directory(top_dir, args.dir_pos_regex)
        for_fq_paths += new_for
        rev_fq_paths += new_rev
        names += new_names

    if args.verbose:
        for f, r, n in zip(for_fq_paths, rev_fq_paths, names):
            print(f, r, n, sep='\t')

    SeqIO.write(renamed_entries(for_fq_paths, names), args.for_out, 'fastq')
    SeqIO.write(renamed_entries(rev_fq_paths, names), args.rev_out, 'fastq')
