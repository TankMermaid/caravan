caravan
=======

Caravan is a fork of SmileTrain that I made in order to
* rework the command-line interface
* change in the intermediate file types (into yaml)
* change the submission style (to explictly use ssub)
* move simple functions to quick perl tools

## caravan requirements
* caravan is developed against Python 3.4 and Perl 5.10.1.
* Biopython

## Suggested workflow for paired-end reads
Newer versions of caravan use the new three file Illumina format: forward reads in one
fastq, reverse reads in a second, and the index reads (aka "barcode reads") in a third.

1.  Rename the reads in the forward and reverse fastq's (using `tools/rename\_fastq.pl`).
2.  Trim primers from the forward and reverse fastq's (using `van.py trim`).
3.  Demultiplex the index reads (using `van.py demultiplex\_fastq`). 
4.  Intersect the forward, reverse, and mapping information (using `van.py intersect3`).
5.  Merge the forward and reverse reads (using `van.py merge`).
6.  Quality filter the merged reads (using `van.py filter`).
7.  Dereplicate and provenance (using `van.py derep`).
8.  Make OTUs.

There are some kludges you could use:
* If you have two-file format, you can extract the barcodes and pop them into a fasta. 
Then you can use `van.py demultiplex_fasta` to work with it.

## Documentation
There is sparse documentation.

## Tests
Unit tests are in the `test` folder. You can run the tests from the top directory 
by running `py.test` or `make`.

Unit tests are mostly incomplete. Sorry. I'm a terrible programmer.
