#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, sys

def triplets_to_entry(triplets, threshold):
    out = []
    for triplet in triplets:
        name, rank, conf = triplet
        name = name.translate({ord('"'): None, ord(' '): '-'})
        if float(conf) < threshold:
            break
        else:
            out.append(name)

    return '_'.join(out)

def parse_line(line):
    fields = line.rstrip().split('\t')
    otu = fields[0].split(';')[0]

    assert fields[1] == ''

    triplets = zip(*[iter(fields[2:])] * 3)
    return (otu, triplets)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fixrank', type=argparse.FileType('r'))
    p.add_argument('--threshold', '-t', type=float, default=0.8)
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()

    print('otu', 'taxonomy', sep='\t', file=args.output)
    for line in args.fixrank:
        otu, triplets = parse_line(line)
        entry = triplets_to_entry(triplets, args.threshold)
        print(otu, entry, sep='\t', file=args.output)
