#!/usr/bin/env python3

'''
command-line interface
'''

import argparse, sys, textwrap
import convert, folder, primers, barcodes, derep, usearch, tax, table, parse, rdp, intersect, qfilter, merge, truncate

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

    p = subparser('folder', help='concatenate fastq\'s from multiple directories')
    p.add_argument('dir_pos_regex', help='regex that identifies the names of folders you want to keep')
    p.add_argument('for_out', type=argparse.FileType('w'), help='forward fastq output')
    p.add_argument('rev_out', type=argparse.FileType('w'), help='reverse fastq output')
    p.add_argument('top_dirs', nargs='+', help='top-level directories (each of which has subfolders, one per samples)')
    p.add_argument('--verbose', '-v', action='store_true', help='show directory names and paired samples')
    p.set_defaults(func=folder.folder)

    p = subparser('primer', help='remove primer')
    p.add_argument('primer')
    p.add_argument('fastq')
    p.add_argument('--max_diffs', '-d', default=2, type=int)
    p.add_argument('--window', '-w', default=20 ,type=int)
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'))
    p.set_defaults(func=primers.PrimerRemover)

    p = subparser('primer2', help='remove second primer', description=textwrap.dedent('''Looks for the primer (which
        can have degenerate bases) at the END of each input sequence. The idea is that this is the primer that would
        be seen at the tail end of the sequence. Unlike the other primer command, if the primer does not match, the
        entry is kept, NOT thrown away.'''))
    p.add_argument('primer')
    p.add_argument('fastq')
    p.add_argument('--max_diffs', '-d', default=2, type=int)
    p.add_argument('--window', '-w', default=20 ,type=int)
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'))
    p.set_defaults(func=primers.SecondPrimerRemover)

    p = subparser('demultiplex', help='assign reads to samples using index reads')
    p.add_argument('barcode_fasta')
    p.add_argument('fastx', help='input fastq (or fasta)')
    p.add_argument('--max_diffs', '-d', type=int, default=1, help='number of barcode mismatches allowed')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output tsv')
    p.add_argument('--input_format', '-t', choices=['fasta', 'fastq'], default='fastq')
    p.set_defaults(func=barcodes.BarcodeMapper)

    p = subparser('intersect', help='intersect fastx and mapping files', description=textwrap.dedent('''Looks for matching entries
        in fastq, fasta, or tsv files. It's finicky: it demands four-line fastq entries, two-line fasta entries, and one-line
        tsv entries. It can only parse files that have ids of the form 'read1', 'read2', etc., where the id would appear in the
        first line of a fastq entry like '@read1', in the first line of a fasta entry like '>read1', or in the first tab-delimited
        field of a tsv entry like 'read1'.'''))
    p.add_argument('--inputs', '-i', required=True, type=argparse.FileType('r'), nargs='+', help='input files', metavar="file")
    p.add_argument('--outputs', '-o', required=True, type=argparse.FileType('w'), nargs='+', help='output files', metavar="file")
    p.set_defaults(func=intersect.intersect)

    p = subparser('merge', help='merge forward and reverse reads')
    p.add_argument('forward', help='forward fastq')
    p.add_argument('reverse', help='reverse fastq')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='merged fastq')
    p.add_argument('--max_diffs', '-d', default=2, type=int, help='maximum mismatches allowed in alignment')
    p.add_argument('--size', '-s', type=int, default=253, help='intended product size? will trash merges that don\'t fit')
    p.add_argument('--size_var', '-v', type=int, default=5, help='allowed variance in product size?')
    p.add_argument('--stagger', '-g', action='store_true', help='staggered reads?')
    p.set_defaults(func=merge.merge)

    p = subparser('filter', help='globally filter low-quality reads')
    p.add_argument('fastq', help='input fastq')
    p.add_argument('--maxee', '-e', default=2.0, type=float, help='discard reads with > E expected errors')
    p.add_argument('--output', '-o', default=sys.stdout, help='output filtered fasta')
    p.add_argument('--output_format', '-t', choices=['fasta', 'fastq'], default='fasta')
    p.set_defaults(func=qfilter.qfilter)

    p = subparser('truncate', help='trim sequences')

    sp = p.add_subparsers()
    sp_len = sp.add_parser('length', help='truncate at a specific length', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_len.add_argument('length', type=int, help='length at which to truncate')
    sp_len.add_argument('fastx', help='input file')
    sp_len.add_argument('--from', '-f', dest='input_format', choices=['fasta', 'fastq'], default='fastq', help='input format')
    sp_len.add_argument('--keep', '-k', action='store_true', help='keep shorter sequences? otherwise discard')
    sp_len.add_argument('--to', '-t', dest='output_format', choices=['fasta', 'fastq'], default='fasta', help='output format')
    sp_len.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='truncated fastx')
    sp_len.set_defaults(func=truncate.length)

    sp_tail = sp.add_parser('tail', help='truncate a tail of low-quality nucleotides', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sp_tail.add_argument('quality', type=int, help='quality score (probably 2 = #)')
    sp_tail.add_argument('fastq', type=argparse.FileType('r'), help='input file')
    sp_tail.add_argument('--to', '-t', dest='output_format', choices=['fasta', 'fastq'], default='fasta', help='output format')
    sp_tail.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='truncated fastq')
    sp_tail.set_defaults(func=truncate.tail)

    p = subparser('derep', help='find unique sequences (and write index file)')
    p.add_argument('fastx', help='input fastx')
    p.add_argument('--min_counts', '-m', type=int, default=0, help='number of times a sequence must appear to be kept')
    p.add_argument('--index', '-i', type=argparse.FileType('w'), help='output yaml index file')
    p.add_argument('--from', '-f', dest='input_format', choices=['fasta', 'fastq'], default='fasta', help='input format')
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
    p.add_argument('--output_base', '-o', default='rdp_X.yaml', help='output filename base')
    p.add_argument('--repl', '-I', default='X', help='pattern in output base to replace with rank-letter')
    p.add_argument('--min_conf', '-m', type=float, default=0.8, help='minimum confidence to assign rank')
    p.set_defaults(func=rdp.FixrankParser.parse_file_all_ranks)

    p = subparser('denovo', help='cluster de novo with usearch')
    p.add_argument('radius', type=float, help='0.0-100.0, recommended at most 3.0 = 97%% identity')
    p.add_argument('fasta', help='input fasta')
    p.add_argument('output', help='representative sequences fasta')
    p.add_argument('--index', '-i', default=None, help='uparse file mapping seqs to otus?')
    p.add_argument('--rename', '-r', default='otu', help='rename otus using this prefix')
    p.add_argument('--force', '-f', action='store_true', help='force a radius above 3.0?')
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

    p = subparser('utax', help='assign taxonomies using utax (and options in utax.yaml)')
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
    p.add_argument('--seq_table', '-q', action='store_true', help='output a seq table?')
    p.add_argument('--seq_table_name', default='seq', help='pre-extension name for seq table (if using -q)')
    p.add_argument('--force', '-f', action='store_true', help='force overwrite of existing files?')
    p.set_defaults(func=table.Tabler.otu_tables)

    p = subparser('seq_table', help='make OTU table from provenances yaml')
    p.add_argument('provenances', help='yaml mapping sequence => {sample => counts}')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output seq table')
    p.set_defaults(func=table.Tabler.seq_table)

    args = parser.parse_args(args)
    opts = vars(args)

    if len(opts) == 0:
        raise RuntimeError("incomplete command line; try adding -h or --help")
    else:
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
