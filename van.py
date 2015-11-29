#!/usr/bin/env python3

'''
command-line interface
'''

import argparse, sys
import convert, primers, barcodes, derep, usearch, tax, table, parse, rdp, intersect, qfilter, merge

def parse_args(args=None):
    '''
    parse command-line arguments

    returns : tuple (func, opts)
        func : function, to be called
        opts : dict, the kwargs for the function
    '''

    parser = argparse.ArgumentParser(description="caravan: a 16S pipeline", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title="commands", metavar='cmd')
    subparsers.required = True

    def subparser(name, **kwargs):
        return subparsers.add_parser(name, formatter_class=argparse.ArgumentDefaultsHelpFormatter, **kwargs)

    p = subparser('convert', help='convert raw fastqs')
    p.add_argument('fastq')
    p.add_argument('--quality', '-q', action='store_true', help='convert quality scores?')
    p.add_argument('--index', '-i', default=None, type=argparse.FileType('w'), help='extract index reads to fasta?')
    p.add_argument('--rename', '-r', action='store_true', help='rename reads sequentially?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fastq')
    p.set_defaults(func=convert.convert_fastq)

    p = subparser('trim', help='remove primer')
    p.add_argument('primer')
    p.add_argument('fastq')
    p.add_argument('--max_diffs', '-d', default=2, type=int)
    p.add_argument('--window', '-w', default=20 ,type=int)
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'))
    p.set_defaults(func=primers.PrimerRemover)

    p = subparser('demultiplex', help='assign reads to samples using index reads')
    p.add_argument('barcode_fasta')
    p.add_argument('fastx', help='input fastq (or fasta)')
    p.add_argument('--max_diffs', '-d', type=int, default=1, help='number of barcode mismatches allowed')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output tsv')
    p.add_argument('--input_format', '-t', choices=['fasta', 'fastq'], default='fastq')
    p.set_defaults(func=barcodes.BarcodeMapper)

    p = subparser('intersect2', help='intersect mapping file and fastq')
    p.add_argument('mapping', type=argparse.FileType('r'), help='tsv from demultiplex')
    p.add_argument('forward', type=argparse.FileType('r'), help='forward reads fastq')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output fastq')
    p.set_defaults(func=intersect.intersect2)

    p = subparser('intersect3', help='intersect mapping file, forward fastq, and reverse fastq')
    p.add_argument('mapping', type=argparse.FileType('r'), help='tsv from demultiplex')
    p.add_argument('forward', type=argparse.FileType('r'), help='forward input fastq')
    p.add_argument('reverse', type=argparse.FileType('r'), help='reverse input fastq')
    p.add_argument('forward_output', type=argparse.FileType('w'), help='forward output fastq')
    p.add_argument('reverse_output', type=argparse.FileType('w'), help='reverse output fastq')
    p.set_defaults(func=intersect.intersect3)

    p = subparser('merge', help='merge forward and reverse reads')
    p.add_argument('forward', help='forward fastq')
    p.add_argument('reverse', help='reverse fastq')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='merged fastq')
    p.add_argument('--max_diffs', '-d', default=2, type=int, help='maximum mismatches allowed in alignment')
    p.add_argument('--size', '-s', type=int, default=253, help='intended product size? will trash merges that don\'t fit')
    p.add_argument('--size_var', '-v', type=int, default=5, help='allowed variance in product size?')
    p.set_defaults(func=merge.merge)

    p = subparser('filter', help='remove low-quality reads')
    p.add_argument('fastq', help='input fastq')
    p.add_argument('--maxee', '-e', default=2.0, type=float, help='discard reads with > E expected errors')
    p.add_argument('--output', '-o', default=sys.stdout, help='output filtered fasta')
    p.add_argument('--output_format', '-t', choices=['fasta', 'fastq'], default='fasta')
    p.set_defaults(func=qfilter.qfilter)

    p = subparser('derep', help='find unique sequences (and write index file)')
    p.add_argument('fastx', help='input fastx')
    p.add_argument('--min_counts', '-m', type=int, default=0, help='number of times a sequence must appear to be kept')
    p.add_argument('--index', '-i', type=argparse.FileType('w'), help='output yaml index file')
    p.add_argument('--input_format', '-t', choices=['fasta', 'fastq'], default='fasta')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='dereplicated fasta')
    p.set_defaults(func=derep.Dereplicator)

    p = subparser('rdp', help='make a mapping file using an RDP fixrank')
    p.add_argument('fixrank', type=argparse.FileType('r'), help='input fixrank file from classifier.jar')
    p.add_argument('level', choices=rdp.rank_abbreviations, help='taxonomic level')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output yaml')
    p.add_argument('--min_conf', '-m', type=float, default=0.8, help='minimum confidence to assign rank')
    p.set_defaults(func=rdp.FixrankParser.parse_file)

    p = subparser('rdpall', help='make mapping files for all taxonomic levels using RDP fixrank')
    p.add_argument('fixrank', type=argparse.FileType('r'), help='input fixrank from classifier.jar')
    p.add_argument('--output_base', '-o', default='rdp_X.yml', help='output filename base')
    p.add_argument('--repl', '-I', default='X', help='pattern in output base to replace with rank-letter')
    p.add_argument('--min_conf', '-m', type=float, default=0.8, help='minimum confidence to assign rank')
    p.set_defaults(func=rdp.FixrankParser.parse_file_all_ranks)

    p = subparser('denovo', help='cluster de novo with usearch')
    p.add_argument('radius', type=float, help='0.0-100.0, recommended at most 3.0 = 97%% identity')
    p.add_argument('fasta', help='input fasta')
    p.add_argument('--output', '-o', default='out.fasta', help='representative sequences fasta')
    p.add_argument('--index', '-i', default=None, help='uparse file mapping seqs to otus?')
    p.add_argument('--rename', '-r', default='otu', help='rename otus using this prefix')
    p.set_defaults(func=usearch.Usearcher().cluster_denovo)

    p = subparser('ref', help='compare a fasta to a reference database')
    p.add_argument('fasta', help='query fasta file')
    p.add_argument('db', help='database fasta file')
    p.add_argument('sid', help='minimum identity between query and database seqs (range: 0.0 to 1.0)')
    p.add_argument('b6', help='blast6 output')
    p.set_defaults(func=usearch.Usearcher().search)

    p = subparser('tax', help='assign taxonomies based on reference alignment')
    p.add_argument('b6', help='blast6 mapping')
    p.add_argument('db', help='pickled {id => taxonomy} dict')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output tax list')
    p.set_defaults(func=tax.TaxAssigner.assign_b6_with_pickled_tax_dict)

    p = subparser('utax', help='assign taxonomies using utax (and options in utax.yml)')
    p.add_argument('fastx')
    p.add_argument('--output', '-o', default='tax.txt')
    p.set_defaults(func=usearch.Usearcher().utax)

    p = subparser('parse', help='parse a blast6 or uparse mapping file to yaml membership file')
    p.add_argument('usearch', type=argparse.FileType('r'), help='blast6 or uparse file')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output yaml')
    p.set_defaults(func=parse.Parser.usearch_to_yaml)

    p = subparser('otu_table', help='make OTU table from membership and provenances yamls')
    p.add_argument('provenances', help='yaml mapping sequence => {sample => counts}')
    p.add_argument('membership', help='yaml mapping sequence => otu')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output otu table')
    p.set_defaults(func=table.Tabler.otu_table)

    p = subparser('otu_tables', help='make multiple OTU tables')
    p.add_argument('provenances', help='yaml mapping sequence => {sample => counts}')
    p.add_argument('memberships', help='yaml files, each mapping sequence => otu', nargs='+')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output_ext', '-o', default='counts', help='output filename extension')
    p.set_defaults(func=table.Tabler.otu_tables)

    p = subparser('seq_table', help='make OTU table from provenances yaml')
    p.add_argument('provenances', help='yaml mapping sequence => {sample => counts}')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output seq table')
    p.set_defaults(func=table.Tabler.seq_table)

    args = parser.parse_args(args)
    opts = vars(args)

    # remove the "func" option from the parser results
    func = opts.pop('func')

    # replace any filenames "-" with stdin
    if "-" in opts.values():
        keys = [k for k in opts if opts[k] == "-"]
        if len(keys) > 1:
            raise RuntimeError("can only accept one '-' argument")
        else:
            opts[keys[0]] = sys.stdin

    return func, opts

if __name__ == '__main__':
    func, opts = parse_args()
    func(**opts)
