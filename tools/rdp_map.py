#!/usr/bin/env python

import argparse, json, re, sys

levels = ['k', 'p', 'c', 'o', 'f', 'g', 's']

def tax_map(entries, level, min_conf):
    level_idx = levels.index(level)

    tmap = {}
    for entry in entries:
        fields = entry.split(';')
        sid = fields[0]

        # look for the first + entry
        try:
            plus_index = fields.index('+')
        except ValueError as e:
            if '-' in fields:
                continue
            else:
                raise RuntimeError("i didn't find + or - in this line: {}".format(entry))

        # get pairs like ['"Bacteroidetes", "99%"']
        raw_pairs = zip(*[iter(fields[plus_index + 1: ])]*2)

        # convert to pairs like ['Bacteroidetes', 99]
        pairs = [[re.sub('"', '', name), int(re.sub('%', '', conf))] for name, conf in raw_pairs]
        assert(pairs.pop(0)[0] == 'Root')

        # construct lineage up to min confidence
        lineage_bits = []
        for i, (name, conf) in enumerate(pairs):
            if i > level_idx or conf < min_conf:
                break

            lineage_bits.append(name)

        lineage = ';'.join(lineage_bits)

        tmap[sid] = lineage

    return tmap
    

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make a mapping json using an rdp allrank file')

    p.add_argument('level', choices=levels)
    p.add_argument('allrank', type=argparse.FileType('r'))
    p.add_argument('--min_conf', '-m', type=int, default=50, help='minimum confidence 0-100')
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout)
    p.add_argument('--header', action='store_true', help='ignore first 7 lines?')

    args = p.parse_args()

    if args.header:
        for i in range(7):
            args.allrank.readline()

    tmap = tax_map(args.allrank, args.level, args.min_conf)

    json.dump(tmap, args.output)
