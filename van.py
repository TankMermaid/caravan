#!/usr/bin/env python

'''
command-line interface
'''

import argparse, sys
import split, check_intersect, primers, barcodes, derep, usearch, tax, table, parse, submit, rdp

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

    p = subparser('check_intersect', help='check that for and rev fastqs have matched read names')
    p.add_argument('forward', help='forward fastq')
    p.add_argument('reverse', help='reverse fastq')
    p.set_defaults(func=check_intersect.IntersectChecker)

    p = subparser('split', help='split fastx into chunks based on size')
    p.add_argument('fastx', help='file to be split')
    p.add_argument('--size', '-s', help='size for each output chunk', default='1.0Gb')
    p.set_defaults(func=split.FastxSplitter)

    p = subparser('merge', help='merge forward and reverse reads')
    p.add_argument('forward', help='forward fastq')
    p.add_argument('reverse', help='reverse fastq')
    p.add_argument('--truncqual', '-q', default=2, type=int, help='truncate the forward and reverse reads at the first Q<=q')
    p.add_argument('--output', '-o', default='merge.fq', help='merged fastq')
    p.add_argument('--size', '-s', type=int, default=None, help='intended product size? will trash merges that don\'t fit')
    p.add_argument('--size_var', '-v', type=int, default=0, help='allowed variance in product size?')
    p.set_defaults(func=usearch.Usearcher().merge)

    p = subparser('find_primers', help='find locations of forward (and reverse) primers')
    p.add_argument('primers_fasta', help='must have entry "forward" (optionally also "reverse")')
    p.add_argument('fastx', help='query fastx')
    p.add_argument('--output', '-o', default='trim.usr', help='output trim file')
    p.add_argument('--max_diffs', '-m', default=1, type=int, help='max mismatches allowed')
    p.set_defaults(func=usearch.Usearcher().primer_search)

    p = subparser('trim_primers', help='trim located primers')
    p.add_argument('trim_file', help='trim file output by find_primers')
    p.add_argument('fastq', help='trimmed fastq')
    p.add_argument('--output', '-o', default='out.fq')
    p.set_defaults(func=primers.PrimerRemover)

    p = subparser('filter', help='remove low-quality reads')
    p.add_argument('fastq', help='input fastq')
    p.add_argument('--truncqual', '-q', default=2, type=int, help='truncate the read at the first position having quality score <= N, so that all remaining Q scores are >N')
    p.add_argument('--maxee', '-e', default=2.0, type=float, help='discard reads with > E expected errors')
    p.add_argument('--output', '-o', default='out.fa', help='output filtered fasta')
    p.set_defaults(func=usearch.Usearcher().filter)

    p = subparser('demultiplex_fastq', help='assign reads to samples using index reads')
    p.add_argument('barcode_fasta')
    p.add_argument('fastx', help='input fastq')
    p.add_argument('--max_diffs', '-m', type=int, default=1, help='number of barcode mismatches allowed')
    p.add_argument('--output', '-o', default='mapped.fq', help='output fastq')
    p.set_defaults(filetype='fastq')
    p.set_defaults(func=barcodes.BarcodeMapper)

    p = subparser('demultiplex_fasta', help='assign reads to samples using index reads')
    p.add_argument('barcode_fasta')
    p.add_argument('fastx', help='input fasta')
    p.add_argument('--max_diffs', '-m', type=int, default=1, help='number of barcode mismatches allowed')
    p.add_argument('--output', '-o', default='mapped.fa', help='output fasta')
    p.set_defaults(filetype='fasta')
    p.set_defaults(func=barcodes.BarcodeMapper)

    p = subparser('both_primers', help='find and trim both primers at once')
    p.add_argument('primers_fasta')
    p.add_argument('fastx', help='query')
    p.add_argument('--output', '-o', default='trim.fa', help='output fasta')
    p.add_argument('--max_diffs', '-m', default=1, type=int, help='max mismatches allowed')
    p.set_defaults(func=usearch.Usearcher().pcr_search)

    p = subparser('derep', help='find unique sequences (and write index file)')
    p.add_argument('fasta', help='input fasta')
    p.add_argument('--min_counts', '-m', type=int, default=0, help='number of times a sequence must appear to be kept')
    p.add_argument('--index', '-i', type=argparse.FileType('w'), help='output json index file')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='dereplicated fasta')
    p.set_defaults(func=derep.Dereplicator)

    p = subparser('rdp', help='make a mapping file using an RDP fixrank')
    p.add_argument('fixrank', type=argparse.FileType('r'), help='input fixrank file from classifier.jar')
    p.add_argument('level', choices=rdp.rank_abbreviations, help='taxonomic level')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output json')
    p.add_argument('--min_conf', '-m', type=float, default=0.8, help='minimum confidence to assign rank')
    p.set_defaults(func=rdp.FixrankParser.parse_file)

    p = subparser('rdpall', help='make mapping files for all taxonomic levels using RDP fixrank')
    p.add_argument('fixrank', type=argparse.FileType('r'), help='input fixrank from classifier.jar')
    p.add_argument('--output_base', '-o', default='rdp_X.json', help='output filename base')
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

    p = subparser('utax', help='assign taxonomies using utax (and options in utax.json)')
    p.add_argument('fastx')
    p.add_argument('--output', '-o', default='tax.txt')
    p.set_defaults(func=usearch.Usearcher().utax)

    p = subparser('parse', help='parse a blast6 or uparse mapping file to json membership file')
    p.add_argument('map_fn')
    p.add_argument('json_fn')
    p.set_defaults(func=parse.Parser.map_to_json)

    p = subparser('otu_table', help='make OTU table from membership and provenances jsons')
    p.add_argument('membership', help='json mapping sequence => otu')
    p.add_argument('provenances', help='json mapping sequence => {sample => counts}')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output otu table')
    p.set_defaults(func=table.Tabler.otu_table)

    p = subparser('seq_table', help='make OTU table from provenances json')
    p.add_argument('provenances', help='json mapping sequence => {sample => counts}')
    p.add_argument('--samples', '-s', help='filename of newline separated sample names in order')
    p.add_argument('--rename', '-r', action='store_true', help='use a two-column sample list to rename them?')
    p.add_argument('--output', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output seq table')
    p.set_defaults(func=table.Tabler.seq_table)

    p = subparser('submit', help='submit jobs')
    p.add_argument('jobs', help='json jobs file')
    p.set_defaults(func=submit.Submitter.submit_jobs)

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