#!/usr/bin/env python

import argparse, json, re, sys

levels = ['k', 'p', 'c', 'o', 'f', 'g', 's']
level_idx = {level: levels.index(level) for level in levels}

def parse_entry(entry_line, min_conf):
    fields = entry_line.split(';')
    sid = fields[0]

    # look for the first + entry
    try:
        plus_index = fields.index('+')
    except ValueError as e:
        if '-' in fields:
            return sid, None
        else:
            raise RuntimeError("i didn't find + or - in this line: '{}'".format(entry_line.rstrip()))

    # get pairs like ['"Bacteroidetes", "99%"']
    raw_pairs = zip(*[iter(fields[plus_index + 1: ])]*2)

    # convert to pairs like ['Bacteroidetes', 99]
    pairs = [[re.sub('"', '', name), int(re.sub('%', '', conf))] for name, conf in raw_pairs]
    assert(pairs.pop(0)[0] == 'Root')
    trimmed_pairs = [p for p in pairs if not p[1] < min_conf]

    return sid, trimmed_pairs

def tax_maps_all(entries, min_conf):
    tmaps = {level: {} for level in levels}

    for entry in entries:
        sid, pairs = parse_entry(entry, min_conf)

        if pairs is None or pairs == []:
            continue

        for level in levels:
            idx = level_idx[level]
            out_length = idx + 1

            if len(pairs) <= out_length:
                out_pairs = pairs
            else:
                out_pairs = pairs[0: idx+1]

            lineage = ';'.join([p[0] for p in out_pairs])
            tmaps[level][sid] = lineage

    return tmaps

def output_filenames(output_base, repl):
    return {l: re.sub(repl, l, output_base) for l in levels}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make a mapping json using an rdp allrank file', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('allrank', type=argparse.FileType('r'))
    p.add_argument('--output', '-o', default='rdp_X.json', help='base output filename')
    p.add_argument('--repl', '-I', default='X', help='replacement string for output base')
    p.add_argument('--min_conf', '-m', type=int, default=50, help='minimum confidence 0-100')
    p.add_argument('--header', '-d', action='store_true', help='ignore first 7 lines?')

    args = p.parse_args()

    if args.header:
        for i in range(7):
            args.allrank.readline()

    tmaps = tax_maps_all(args.allrank, args.min_conf)

    fns = output_filenames(args.output, args.repl)
    for level in levels:
        with open(fns[level], 'w') as f:
            json.dump(tmaps[level], f)