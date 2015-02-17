'''
usearch wrapper. implements a nice parser that has more helpful help messages.
'''

import argparse, os, sys, subprocess, StringIO

class Usearcher:
    def __init__(self, dry=False):
        self.dry = dry

    def run(self, cmd):
        cmd = [str(x) for x in cmd]
        cmd_string = " ".join(cmd)

        if self.dry:
            print cmd_string
        else:
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

    def merge(self, forward, reverse, output, truncqual=None):
        cmd = ['usearch', '-fastq_mergepairs', forward, '-reverse', reverse, '-fastqout', output]

        if truncqual is not None:
            cmd += ['-fastq_truncqual', truncqual]

        self.run(cmd)

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