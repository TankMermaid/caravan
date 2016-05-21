#!/usr/bin/env python

import argparse, json, re

levels = ['k', 'p', 'c', 'o', 'f', 'g', 's']

def truncate_taxonomy(tax, level):
    # remove the numbers, so Bacteria(100) becomes Bacteria etc.
    tax = re.sub('\([\d\.]+\)', '', tax)

    # remove any silly quotations
    tax = re.sub('"', '', tax)

    fields = tax.split(',')
    level_idx = min([levels.index(level), len(fields) - 1]) + 1

    return ','.join(fields[0: level_idx])

def tax_map(entries, level):
    tmap = {}
    for entry in entries:
        fields = entry.split()
        sid = fields[0]
        tax = fields[1]
        sid = re.sub(';.*$', '', sid)
        trim_tax = truncate_taxonomy(tax, level)
        tmap[sid] = trim_tax

    return tmap
    

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make a mapping json using a utax')

    p.add_argument('level', choices=levels)
    p.add_argument('utax')
    p.add_argument('output')

    args = p.parse_args()

    with open(args.utax) as tax_entries:
        tmap = tax_map(tax_entries, args.level)

    with open(args.output, 'w') as f:
        json.dump(tmap, f)