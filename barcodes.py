'''
Demultiplex reads by mapping the barcode read to the sample names from a barcode mapping
file.

Given a tab-separated barcode mapping file like
    donor1_day5   ACGT

the first read mapping to that barcode, say
    @OURSEQ:lolapalooza1234#ACGT/1
    AACCGGTT
    +
    abcdefgh

becomes output like
    @sample=donor1_day5;1
    AACCGGTT
    +whatever
    abcdefgh
    
where the ;1 means it's the first read that mapped to donor1_day5.
'''

import usearch_python.primer, util
import sys, argparse, string, itertools, re, csv
from Bio import SeqIO

class BarcodeMapper:
    def __init__(self, barcode_fn, fastq, max_barcode_diffs, output=sys.stdout, run=True):
        self.barcode_fn = barcode_fn
        self.fastq = fastq
        self.max_barcode_diffs = max_barcode_diffs
        self.output = output

        if run:
            self.run()

    def run(self):
        # parse the barcode mapping file
        with open(self.barcode_fn, 'r') as f:
            self.barcode_map = self.barcode_file_to_dictionary(f)

        # get a set of reads
        for record in self.renamed_fastq_records(self.fastq, self.barcode_map, self.max_barcode_diffs):
            SeqIO.write(record, self.output, 'fastq')

    @staticmethod
    def barcode_file_to_dictionary(barcode_lines):
        '''parse a barcode mapping file into a dictionary {barcode: sample}'''
        barcode_map = {}
        for line in barcode_lines:
            sample, barcode = line.split()
            barcode_map[barcode] = sample

        return barcode_map

    @staticmethod
    def best_barcode_match(known_barcodes, barcode):
        '''
        Find the best match between a known barcode a list of known barcodes

        Parameters
        known_barcodes : sequence of iterator of sequences
            list of known barcodes
        barcode : string
            the barcode read to be matched against the known barcodes

        Returns
        min_mismatches : int
            number of mismatches in the best alignment
        best_known_barcode : string
            known barcode that aligned best
        '''
        
        # get a list of pairs (n_mismatches, known_barcode)
        #n_mismatches = lambda known_barcode: util.mismatches(barcode, known_barcode, 1)[1]
        n_mismatches = lambda known_barcode: usearch_python.primer.MatchPrefix(barcode, known_barcode)

        alignments = [(n_mismatches(known_barcode), known_barcode) for known_barcode in known_barcodes]

        # find the alignment that has the minimum number of mismatches
        min_mismatches, best_known_barcode = min(alignments, key=lambda x: x[0])

        return min_mismatches, best_known_barcode

    @staticmethod
    def parse_barcode(record):
        '''
        Extract the barcode read and direction from a BioPython SeqRecord
        
        Parameters
        record : SeqRecord
            fastq record
        
        returns : tuple
            (barcode read, read direction), where direction is either '1' or '2'
        '''
        
        # match, e.g. @any_set_of_chars#ACGT/1 -> ACGT
        m = re.match(".*#([ACGTN]+)/(\d)$", record.id)

        if m is None:
            raise RuntimeError("fastq id did not match expected format: %s" %(record.id))

        # pull out the read and direction from the match
        barcode_read = m.group(1)
        read_direction = m.group(2)
        
        if read_direction not in ['1', '2']:
            raise RuntimeError('read direction not 1 or 2: %s' %(record.id))
        
        return (barcode_read, read_direction)
    
    @classmethod
    def renamed_fastq_records(cls, fastq, barcode_map, max_barcode_diffs):
        '''
        Rename the read IDs in a fastq file with the corresponding sample name. Get the barcode
        read right from the ID line, look it up in the barcode map, and pick the best match.

        Parameters
        fastq : filename or filehandle
            input
        barcode_map : dictionary
            entries are {barcode: sample_name}
        max_barcode_diffs : int
            maximum number of mismatches between a barcode read and known barcode before throwing
            out that read

        yields : SeqRecord
            fastq records
        '''

        # keep track of the computations where we align the barcode read to the known barcodes
        barcode_read_to_sample = {}
        sample_counts = {}

        for record in SeqIO.parse(fastq, 'fastq'):
            # look for the barcode from the read ID line
            barcode_read, read_direction = cls.parse_barcode(record)
            
            if barcode_read in barcode_read_to_sample:
                sample = barcode_read_to_sample[barcode_read]
                sample_counts[sample] += 1
            else:
                # try aligning to every known barcode
                n_mismatches, best_known_barcode = cls.best_barcode_match(barcode_map.keys(), barcode_read)

                # if the match was good, assign that barcode read to the sample that the best read
                # matches
                if n_mismatches > max_barcode_diffs:
                    continue
                else:
                    # get the name for this sample; record which sample we mapped this barcode
                    # read to
                    sample = barcode_map[best_known_barcode]
                    barcode_read_to_sample[barcode_read] = sample
                    sample_counts[sample] = 1

            record.id = "sample=%s;%d/%s" %(sample, sample_counts[sample], read_direction)
            
            # expunge other parts of title
            record.description = ''
            yield record


class BarcodeCounter(BarcodeMapper):
    '''
    For counting the number of each barcode that appear in the first N entries of a fastq.
    Currently it only looks for exact matches.
    '''

    def __init__(self, barcode_fn, fastq, max_entries, output=sys.stdout, run=True):
        self.barcode_fn = barcode_fn
        self.fastq = fastq
        self.max_entries = max_entries
        self.output = output

        if run:
            self.run()

    def run(self):
        with open(self.barcode_fn) as f:
            self.barcode_map = self.barcode_file_to_dictionary(f)
    
        self.counts = self.count_barcodes(self.fastq, self.barcode_map, self.max_entries)
    
        self.write_counts_report(self.counts, self.output)

    @classmethod
    def count_barcodes(cls, fastq, barcode_map, max_entries):
        '''count the number of appearances of each barcode in a fastq'''
    
        # initialize the count dictionary
        counts = {sample: 0 for sample in barcode_map.values()}
        counts['mapped'] = 0
        counts['total'] = 0
    
        for i, record in enumerate(SeqIO.parse(fastq, 'fastq')):
            if i >= max_entries:
                break
            
            barcode_read, read_direction = cls.parse_barcode(record)
            if barcode_read in barcode_map:
                sample = barcode_map[barcode_read]
                counts[sample] += 1
                counts['mapped'] += 1
            
            counts['total'] += 1
        
        return counts

    @staticmethod
    def write_counts_report(counts, f):
        '''log the counts in a filehandle f'''
    
        # sort the results by abundance
        sorted_counts = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)
        
        percent = lambda n: "%.2f%%" %(float(n) / counts['total'] * 100)
        
        w = csv.writer(f, delimiter='\t', quoting=csv.QUOTE_NONE)
        w.writerow(['sample', '# reads', '% reads'])
            
        for sample, n in sorted_counts:
            w.writerow([sample, str(n), percent(n)])