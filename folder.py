# author: scott olesen <swo@mit.edu>

'''
Functionality for concatenating fastq's from multiple, already-demultiplexed directories
into a single, renamed fastq. For example, a "top-level directory" 160111Foo might have a
structure like:

    160111Foo/
    ├── 160119_AJPVA_1737T_L1_1_sequence_unmapped_barcodes.fastq
    ├── 160119_AJPVA_1737T_L1_2_sequence_unmapped_barcodes.fastq
    ├── D16-100001-1737T
    │   ├── 160111Foo_D16-100001_1_sequence_contam_report.txt
    │   ├── 160111Foo_D16-100001_1_sequence.fastq
    │   ├── 160111Foo_D16-100001_1_sequence_fastqc.html
    │   ├── 160111Foo_D16-100001_1_sequence_fastqc.zip
    │   ├── 160111Foo_D16-100001_1_sequence_tagcount.xls
    │   ├── 160111Foo_D16-100001_2_sequence.fastq
    │   ├── 160111Foo_D16-100001_2_sequence_fastqc.html
    │   ├── 160111Foo_D16-100001_2_sequence_fastqc.zip
    │   ├── 160111Foo_D16-100001_2_sequence_tagcount.xls
    │   ├── 160111Foo_D16-100001_phiX_bestmap.bam
    │   ├── 160111Foo_D16-100001_phiX_bestmap.sam
    │   └── 160111Foo_D16-100001_phiX_bestmap_stats.json
    ├── D16-100002-1737T
    │   ├── 160111Foo_D16-100002_1_sequence_contam_report.txt
    │   ├── 160111Foo_D16-100002_1_sequence.fastq
    │   ├── 160111Foo_D16-100002_1_sequence_fastqc.html
    │   ├── 160111Foo_D16-100002_1_sequence_fastqc.zip
    │   ├── 160111Foo_D16-100002_1_sequence_tagcount.xls
    │   ├── 160111Foo_D16-100002_2_sequence.fastq
    │   ├── 160111Foo_D16-100002_2_sequence_fastqc.html
    │   ├── 160111Foo_D16-100002_2_sequence_fastqc.zip
    │   ├── 160111Foo_D16-100002_2_sequence_tagcount.xls
    │   ├── 160111Foo_D16-100002_phiX_bestmap.bam
    │   ├── 160111Foo_D16-100002_phiX_bestmap.sam
    │   └── 160111Foo_D16-100002_phiX_bestmap_stats.json
    └── infosite-1737T
        ├── 160119_AJPVA_1737T_infosite.rst
            ├── images
                │   ├── 160119_AJPVA_1737T_tagcount_kmer_counts.png
                    │   └── PPPQC_1737T.jpg
                        ├── index.html
                            └── index.rst

This script can take a list of these top-level directories and a "positive directory
regex" that will identify the wanted folders. In this case, "^D16-" might be a good
regex, since it will ignore the "infosite-" directory.

The script will look for pairs of fastq's with _1_ and _2_ in each folder. It will then
concatenate all the _1_'s into the forward out, renaming each entry into the form
"@sample=[sample_name];[number of sequence in that sample". In the above situation,
the very first reads in the forward output and reverse output will have the name

    @sample=D16-100001-1737T;1

These names are compatible with caravan's expectations about a naming scheme when
dereplicating.
'''

import os, re
from Bio import SeqIO

def matches_regexes(query, pos, neg):
    '''does query match positive regex and not match negative?'''
    return re.search(pos, query) is not None and re.search(neg, query) is None

def first_match(lst, regex, only=True):
    candidates = [elt for elt in lst if re.search(regex, elt)]

    if len(candidates) == 0:
        raise RuntimeError("could not find regex '{}' in list: {}".format(regex, lst))

    if only:
        assert len(candidates) == 1

    return candidates[0]

def fastq_paths_in_directory(target_dir, dir_pos_regex, dir_neg_regex='infosite', for_fq_regex='_1_.+\.fastq', rev_fq_regex='_2_.+\.fastq'):
    structure = {x: {'dirs': y, 'files': z} for x, y, z in os.walk(target_dir)}
    assert target_dir in structure

    top_dir_names = [d for d in structure[target_dir]['dirs'] if matches_regexes(d, dir_pos_regex, dir_neg_regex)]
    top_dir_paths = [os.path.join(target_dir, d) for d in top_dir_names]

    # check that all these top directories don't have subdirectories
    for d in top_dir_paths:
        assert structure[d]['dirs'] == []

    # find the forward and reverse reads for each 
    for_fq_names = [first_match(structure[d]['files'], for_fq_regex) for d in top_dir_paths]
    rev_fq_names = [first_match(structure[d]['files'], rev_fq_regex) for d in top_dir_paths]

    for_fq_paths = [os.path.join(d, fn) for d, fn in zip(top_dir_paths, for_fq_names)]
    rev_fq_paths = [os.path.join(d, fn) for d, fn in zip(top_dir_paths, rev_fq_names)]

    return (for_fq_paths, rev_fq_paths, top_dir_names)

def renamed_entries(fns, samples):
    for fn, sample in zip(fns, samples):
        for i, record in enumerate(SeqIO.parse(fn, 'fastq')):
            record.id = "sample={};{}".format(sample, i + 1)
            record.description = ""
            yield record

def folder(dir_pos_regex, for_out, rev_out, top_dirs, verbose=False):
    for_fq_paths = []
    rev_fq_paths = []
    names = []
    for top_dir in top_dirs:
        new_for, new_rev, new_names = fastq_paths_in_directory(top_dir, dir_pos_regex)
        for_fq_paths += new_for
        rev_fq_paths += new_rev
        names += new_names

    if verbose:
        for f, r, n in zip(for_fq_paths, rev_fq_paths, names):
            print(f, r, n, sep='\t')

    SeqIO.write(renamed_entries(for_fq_paths, names), for_out, 'fastq')
    SeqIO.write(renamed_entries(rev_fq_paths, names), rev_out, 'fastq')
