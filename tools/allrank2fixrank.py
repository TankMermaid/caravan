#!/usr/bin/env python3

import argparse, sys

fixed_ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus']

def parse_fields(fields):
    otu = fields[0]
    sign = fields[1]
    rest = fields[2:]
    assert(len(rest) % 3 == 0)
    triplets = zip(*[iter(rest)] * 3)

    return otu, sign, triplets

def convert_triplets(triplets):
    prev_taxon = "NA"
    prev_conf = 0.0
    dat = {t[1]: (t[0], t[2]) for t in triplets}

    out = []
    for rank in reversed(fixed_ranks):
        if rank in dat:
            taxon, conf = dat[rank]
            out = [taxon, rank, conf] + out
            prev_taxon = taxon
            prev_conf = conf
        else:
            out = [prev_taxon, rank, prev_conf] + out

    return out

def convert_fields(fields):
    otu, sign, triplets = parse_fields(fields)
    new_triplets = convert_triplets(triplets)
    return otu, sign, new_triplets

def convert_line(line):
    fields = line.rstrip().split("\t")
    otu, sign, new_triplets = convert_fields(fields)
    new_fields = [otu, sign] + new_triplets
    return "\t".join([str(x) for x in new_fields])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('allrank', type=argparse.FileType('r'))
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()

    for line in args.allrank:
        print(convert_line(line), file=args.output)