#!/usr/bin/env python3

import argparse, sys

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('seq_table', type=argparse.FileType('r'))
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout)
    args = p.parse_args()

    header = args.seq_table.readline().rstrip()
    samples = header.split("\t")[1:]

    for line in args.seq_table:
        fields = line.rstrip().split("\t")
        otu = fields[0]
        counts = fields[1:]

        # if this sequence has no counts, skip it
        if sum([float(x) for x in counts]) > 0:
            print(otu + ":", file=args.output)

            for sample, count in zip(samples, counts):
                if float(count) > 0:
                    print("  {}: {}".format(sample, count), file=args.output)