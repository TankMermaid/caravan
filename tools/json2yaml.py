#!/usr/bin/env python3
#
# author: scott olesen <swo@mit.edu>

import argparse, json, yaml, sys

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="convert json file to yaml format", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('json', type=argparse.FileType('r'))
    p.add_argument('--output', '-o', type=argparse.FileType('w'), default=sys.stdout, help='output yaml')
    p.add_argument('--flow_style', '-s', default='False', choices=['None', 'True', 'False'], help='value to pass to yaml.dump\'s default_flow_style')
    args = p.parse_args()

    # convert string into the actual python thingy
    args.flow_style = eval(args.flow_style)

    data = json.load(args.json)
    print(yaml.dump(data, default_flow_style=args.flow_style), file=args.output, end="")
