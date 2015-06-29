#!/usr/bin/env python

'''
Create OTU tables using information from
    * a membership .yml that has a hash sequence => OTU
    * a provenances .yml that has a hash sequence => {sample => counts}
'''

import re, sys, argparse, yaml, warnings
from operator import itemgetter

class Tabler:
    @staticmethod
    def table(membership, provenances, output, samples=None, rename=False, ignore_unmapped_otus=True):
        # populate the tables
        table = {}  # {sample => {otu => counts}}
        otu_abunds = {} # {otu => counts across samples}

        for seq in provenances:
            if ignore_unmapped_otus:
                if seq not in membership:
                    continue

                otu = membership[seq]
            else:
                raise RuntimeError("not ignore umapped OTUs not implemented")

            if otu not in otu_abunds:
                otu_abunds[otu] = 0

            for sample in provenances[seq]:
                c = provenances[seq][sample]

                # add counts to the table
                if sample not in table:
                    table[sample] = {otu: c}
                elif otu not in table[sample]:
                    table[sample][otu] = c
                else:
                    table[sample][otu] += c

                # add counts to abundances
                otu_abunds[otu] += c

        # get samples, if any
        if samples is None:
            samples = sorted(table.keys())
        else:
            if rename:
                # grab sample fields and extract new sample names
                old_new_samples = [line.split() for line in open(samples)]
                samples = [x[1] for x in old_new_samples]

                # check that, for every new name, the old name is actually in the list
                for old_s, new_s in old_new_samples:
                    if old_s not in table:
                        warnings.warn('sample "{}" in rename list (alias "{}") but not the provenances'.format(old_s, new_s), stacklevel=2)
                        samples.remove(new_s)

                # update the table with the new names that were found in the provenance
                for old_s, new_s in old_new_samples:
                    if new_s in samples:
                        table[new_s] = table.pop(old_s)
            else:
                samples = [line.split()[0] for line in open(samples)]
            
        # sort otu names by decreasing abundance
        sorted_otus_abunds = sorted(otu_abunds.items(), key=itemgetter(1), reverse=True)

        # first, output the header/sample line
        output.write("\t".join(['OTU_ID'] + samples) + "\n")

        for otu, counts in sorted_otus_abunds:
            output.write("\t".join([otu] + [str(table[sample].get(otu, 0)) for sample in samples]) + "\n")

    @classmethod
    def otu_table(cls, membership, provenances, output, samples=None, rename=False):
        # get the index
        with open(provenances) as f:
            provenances_dict = yaml.load(f)

        with open(membership) as f:
            membership_dict = yaml.load(f)

        cls.table(membership_dict, provenances_dict, output, samples, rename)

    @classmethod
    def seq_table(cls, provenances, output, samples=None, rename=False):
        # get the index
        with open(provenances) as f:
            provenances_dict = yaml.load(f)

        # make up membership as self => self
        membership_dict = {seq: seq for seq in provenances_dict}

        cls.table(membership_dict, provenances_dict, output, samples, rename)
