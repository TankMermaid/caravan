'''
Create OTU tables using information from
    * a membership .yml that has a hash sequence => OTU
    * a provenances .yml that has a hash sequence => {sample => counts}
'''

import re, sys, argparse, yaml, warnings, os.path
from operator import itemgetter

class Tabler:
    @staticmethod
    def table(provenances, membership, output, samples=None, rename=False, ignore_unmapped_otus=True):
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

            obj_type = type(provenances[seq])
            if obj_type is not dict:
                raise RuntimeError("malformed provenances file: sequence '{}' points to a non-dict object, type {}, with content {}".format(seq, obj_type, provenances[seq]))

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

        # sort otu names by decreasing abundance, then by name
        sorted_otus_abunds = sorted(otu_abunds.items(), key=lambda x: (-x[1], x[0]))

        # first, output the header/sample line
        output.write("\t".join(['OTU_ID'] + samples) + "\n")

        for otu, counts in sorted_otus_abunds:
            output.write("\t".join([otu] + [str(table[sample].get(otu, 0)) for sample in samples]) + "\n")

    @staticmethod
    def check_membership_format(x, fn):
        for val in x.values():
            if type(val) is not str:
                raise RuntimeError("membership file {} is not a hash of strings".format(fn))
            break

    @staticmethod
    def check_provenances_format(x, fn):
        for val in x.values():
            if type(val) is not dict:
                raise RuntimeError("provenances file {} is not a hash of hashes".format(fn))
            break

    @classmethod
    def otu_table(cls, provenances, membership, output, samples=None, rename=False):
        # get the index
        with open(provenances) as f:
            provenances_dict = yaml.load(f)

        cls.check_provenances_format(provenances_dict, provenances)

        with open(membership) as f:
            membership_dict = yaml.load(f)

        cls.table(provenances_dict, membership_dict, output, samples, rename)

    @classmethod
    def otu_tables(cls, provenances, memberships, output_ext, samples=None, rename=False, seq_table=False, seq_table_name="seq", force=False):
        # check that all the output locations are OK first
        base_output_names = [os.path.splitext(os.path.basename(m))[0] for m in memberships]
        output_names = [base + '.' + output_ext for base in base_output_names]
        existing_files = [fn for fn in output_names if os.path.exists(fn)]

        if seq_table:
            seq_table_output_name = seq_table_name + '.' + output_ext
            if os.path.exists(seq_table_output_name):
                existing_files.append(seq_table_output_name)

        if not force and len(existing_files) > 0:
            raise RuntimeError('some output files would be overwritten: {}'.format(existing_files))

        # load in the provenances first, since it will apply to all the memberships
        with open(provenances) as f:
            provenances_dict = yaml.load(f)

        cls.check_provenances_format(provenances_dict, provenances)

        # load each membership file and do its output
        for membership, output_fn in zip(memberships, output_names):
            with open(membership) as f:
                membership_dict = yaml.load(f)

            cls.check_membership_format(membership_dict, membership)

            with open(output_fn, 'w') as f:
                cls.table(provenances_dict, membership_dict, f, samples, rename)

        # make the seq table if called for
        if seq_table is not None:
            membership_dict = {seq: seq for seq in provenances_dict}
            with open(seq_table_output_name, 'w') as output:
                cls.table(provenances_dict, membership_dict, output, samples, rename)

    @classmethod
    def seq_table(cls, provenances, output, samples=None, rename=False):
        # get the index
        with open(provenances) as f:
            provenances_dict = yaml.load(f)

        # make up membership as self => self
        membership_dict = {seq: seq for seq in provenances_dict}

        cls.table(provenances_dict, membership_dict, output, samples, rename)
