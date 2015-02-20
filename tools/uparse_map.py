#!/usr/bin/env python

import argparse, json, re

def uparse_map(entries):
    umap = {}
    for entry in entries:
        raw_sid, match_type, pid1, pid2, otu = entry.split()
        sid = re.sub(';.*$', '', raw_sid)
        if match_type in ['otu', 'match']:
            umap[sid] = otu

    return umap
    

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='make a mapping json using a uparse file')

    p.add_argument('uparse')
    p.add_argument('output')

    args = p.parse_args()

    with open(args.uparse) as uparse_entries:
        umap = uparse_map(uparse_entries)

    with open(args.output, 'w') as f:
        json.dump(umap, f)