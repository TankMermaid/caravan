caravan
=======

Caravan is a fork of SmileTrain that I made in order to
* rework the command-line interface
* change in the intermediate file types (into yaml)
* change the submission style (to explictly use ssub)
* move simple functions to quick perl tools

**Caveat**: This readme has fallen a little behind, especially when it comes to
the initial data format and the possibilities for merging.

## caravan requirements
* caravan is developed against Python 3.4 and Perl 5.10.1.
* Biopython

## Getting started

First, you'll need to have your data in a standard input format:

- All fastq files must be in Illumina 1.8+ format (also known as Sanger, Phred+33, or ASCII offset 33). You can convert fastq files in older Illumina format (i.e., Phred+64) using `van.py convert`.
- Intersecting multiple files requires that they have IDs like `readXXX`, where `XXX` are monotonically increasing integers. You can get your fastq files into this format using `van.py convert --rename`.
- When dereplicating, if you want provenience data (i.e., to make an OTU table), then the reads must have IDs with a `sample=XXX` field. A field is separate from other parts of the ID my semicolons. The preferred read ID format is like `read1234;sample=sample1`.

## Suggested workflow for paired-end reads
Newer versions of caravan use the new three file Illumina format: forward reads in one
fastq, reverse reads in a second, and the index reads (aka "barcode reads") in a third.

1.  Trim primers from the forward and reverse fastq's (using `van.py primer`).
2.  Demultiplex the index reads (using `van.py demultiplex`). 
3.  Intersect the forward, reverse, and mapping information (using `van.py intersect`).
4.  Merge the forward and reverse reads (using `van.py merge`).
5.  Quality filter the merged reads (using `van.py filter`).
6.  Dereplicate and provenance (using `van.py derep`).
7.  Make OTUs.

If you do not have paired-end reads, you'll probably want to use `van.py truncate` to
trim sequences by their quality or length.

## Documentation
There is sparse documentation. Sorry. But the `--help` option is usually pretty explanatory.

## Tests
Unit tests are in the `test` folder. You can run the tests from the top directory 
by running `py.test` or `make`.

Unit tests are mostly incomplete. Sorry. I'm a terrible programmer.
