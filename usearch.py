'''
usearch wrapper. implements a nice parser that has more helpful help messages.
'''

import argparse, os, sys, subprocess, json, os.path
from Bio import SeqRecord, SeqIO, Seq

class Usearcher:
    def __init__(self, debug=False):
        '''debug mode doesn't run the command'''
        self.debug = debug

    def run(self, cmd):
        '''
        other methods construct commands. this method calls the command
        (unless debug is on) and returns the command arguments
        '''

        cmd = [str(x) for x in cmd]
        cmd_string = " ".join(cmd)

        if not self.debug:
            # run command and grab stdout
            try:
                self.out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                # if there is an error, look for the line just below "fatal error"
                err_lines = e.output.split("\n")
                
                fatal_error = '---Fatal error---'
                invalid_command = 'Invalid command line'
                if fatal_error in err_lines:
                    fatal_i = err_lines.index(fatal_error)
                    desc_line = err_lines[fatal_i + 1]
                elif invalid_command in err_lines:
                    err_i = err_lines.index(invalid_command)
                    desc_line = err_lines[err_i + 1]
                
                raise RuntimeError("when running with command '%s', usearch failed with message: '%s'" %(cmd_string, desc_line))

        return cmd

    def merge(self, forward, reverse, output, truncqual=None, size=None, size_var=0):
        cmd = ['usearch', '-fastq_mergepairs', forward, '-reverse', reverse, '-fastqout', output]

        if truncqual is not None:
            cmd += ['-fastq_truncqual', truncqual]

        if size is not None:
            min_len = size - size_var
            max_len = size + size_var
            cmd += ['-fastq_minmergelen', min_len, '-fastq_maxmergelen', max_len]

        return self.run(cmd)

    def filter(self, fastq, output, truncqual=None, maxee=None):
        cmd = ['usearch', '-fastq_filter', fastq, '-fastaout', output]

        if truncqual is not None:
            cmd += ['-fastq_truncqual', truncqual]

        if maxee is not None:
            cmd += ['-fastq_maxee', maxee]

        self.run(cmd)

    def cluster_denovo(self, fasta, radius, output, index=None, rename=None):
        cmd = ['usearch', '-cluster_otus', fasta, '-otu_radius_pct', radius, '-otus', output]

        if index is not None:
            cmd += ['-uparseout', index]

        if rename is not None:
            cmd += ['-relabel', rename]

        self.run(cmd)

    def cluster_smallmem(self, fasta, otuid, centroids=None, uc=None):
        if isinstance(otuid, int):
            otuid = '.%d' %(otuid)
        
        cmd = ['usearch', '-cluster_smallmem', fasta, '-id', otuid]

        if centroids is not None:
            cmd += ['-centroids', centroids]

        if uc is not None:
            cmd += ['-uc', uc]

        self.run(cmd)

    def length_sort(self, fasta, output):
        cmd = ['usearch', '-sortbylength', fasta, '-output', output]
        self.run(cmd)

    def chimera_ref(self, fasta, db, fasta_no_chimeras=None, fasta_chimeras=None, strand='plus'):
        cmd = ['usearch', '-uchime_ref', fasta, '-db', db, '-strand', strand]

        if fasta_no_chimeras is not None:
            cmd += ['-nonchimeras', fasta_no_chimeras]

        if fasta_chimeras is not None:
            cmd += ['-chimeras', fasta_chimeras]

        self.run(cmd)

    def chimera_denovo(self, fasta, fasta_no_chimeras=None, fasta_chimeras=None, strand='plus'):
        cmd = ['usearch', '-uchime_denovo', fasta, '-strand', strand]

        if fasta_no_chimeras is not None:
            cmd += ['-nonchimeras', fasta_no_chimeras]

        if fasta_chimeras is not None:
            cmd += ['-chimeras', fasta_chimeras]

        self.run(cmd)

    def search(self, fasta, db, sid, b6, not_matched=None, fasta_pairs=None, strand='both', no_hits=True):
        cmd = ['usearch', '-usearch_global', fasta, '-db', db, '-id', sid, '-strand', strand, '-blast6out', b6]

        if not_matched is not None: cmd += ['-notmatched', not_matched]
        if fasta_pairs is not None: cmd += ['-fastapairs', fasta_pairs]
        if no_hits: cmd += ['-output_no_hits']

        self.run(cmd)

    def short_search(self, fastx, db, out, max_diffs=1, sid=0.0, strand='plus'):
        '''usearch for short, well-aligned sequences. e.g., for index reads'''

        cmd = ['usearch', '-search_global', fastx, '-db', db, '-strand', strand, '-userout', out, '-id', sid,
            '-fulldp', '-maxdiffs', max_diffs, '-leftjust', '-rightjust', '-userfields', 'query+target']

        self.run(cmd)

    def primer_search(self, primers_fasta, fastx, output, max_diffs):
        '''search for oligos like primers'''

        # first check that the primers fasta has the expected content
        primers = [r for r in SeqIO.parse(primers_fasta, 'fasta')]

        if primers[0].id != 'forward':
            raise RuntimeError("primers fasta {} should have a first entry 'forward', not '{}'".format(primers_fasta, primers[0].id))

        if len(primers) == 1:
            userfields = 'query+qhi'
            strand = 'plus'
        elif len(primers) == 2:
            if primers[1].id != 'reverse':
                raise RuntimeError("primers fasta {} should have a second entry 'reverse', not '{}'".format(primers_fasta, primers[1].id))
            userfields = 'query+target+qstrand+qlo+qhi'
            strand = 'both'
        else:
            raise RuntimeError("primers fasta {} should have at most two entries".format(primers_fasta))

        cmd = ['usearch', '-search_oligodb', fastx, '-db', primers_fasta, '-userout', output,
            '-userfields', userfields, '-strand', strand, '-maxdiffs', max_diffs, '-maxhits', 2]

        self.run(cmd)

    def pcr_search(self, primers_fasta, fastx, output, max_diffs, strip=True):
        '''search for a pair of oligos'''

        cmd = ['usearch', '-search_pcr', fastx, '-db', primers_fasta, '-strand', 'both',
            '-maxdiffs', max_diffs, '-ampout', output]

        if strip:
            cmd.append('-pcr_strip_primers')

        self.run(cmd)

    def utax(self, fastx, output):
        opts_fn = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'utax.json')
        with open(opts_fn) as f:
            opts = json.load(f)

        cmd = ['usearch', '-utax', fastx, '-db', opts['db'], '-taxconfs', opts['tc'], '-tt', opts['tt'],
            '-utaxout', output]

        self.run(cmd)