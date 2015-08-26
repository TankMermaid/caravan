#!/usr/bin/env python

import argparse, json, re

levels = ['k', 'p', 'c', 'o', 'f', 'g', 's']

def truncate_taxonomy(tax, level):
    # remove the weird prefixes
    tax = re.sub('[kpcofgs]__', '', tax)

    fields = [f for f in tax.split('; ') if f != '']
    level_idx = min([levels.index(level), len(fields) - 1]) + 1

    return ','.join(fields[0: level_idx])

def tax_map(b6_entries, tax_entries, level):
    ggid2sid = {}
    for entry in b6_entries:
        fields = entry.split()
        sid = fields[0]
        sid = re.sub(';.*$', '', sid)
        ggid = fields[1]
        ggid2sid[ggid] = sid

    sid2tax = {}
    for entry in tax_entries:
        fields = entry.rstrip().split("\t")
        ggid = fields[0]
        tax = truncate_taxonomy(fields[1], level)

        if ggid in ggid2sid:
            sid = ggid2sid.pop(ggid)
            sid2tax[sid] = tax

    return sid2tax
    

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make a mapping json using greengenes taxonomy files')

    p.add_argument('level', choices=levels)
    p.add_argument('gg_tax')
    p.add_argument('b6', help='output from van ref')
    p.add_argument('output')

    args = p.parse_args()

    with open(args.b6) as b6_entries:
        with open(args.gg_tax) as tax_entries:
            tmap = tax_map(b6_entries, tax_entries, args.level)

    with open(args.output, 'w') as f:
        json.dump(tmap, f)
