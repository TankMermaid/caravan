'''
remove primers from a fastq using the trim file

trim files have lines like...

for forward primers only:
read_id start_idx end_idx

for reverse primers too:
read_id primer_id strand start_idx end_idx

this script looks through the fastq and the trim file and either
trims the fastq sequence (or throws it out if not all the primers
hit)
'''

from Bio import SeqRecord, SeqIO, Seq

class TrimmedRecords():
    '''class for iterating through fastq and trimming'''

    def __init__(self, trims, fastq_records, reverse):
        '''
        trims : dict
            if forward, {record id => primer index}
            if reverse, {record id => [start, end index]}
        '''
        self.trims = trims
        self.fastq_records = fastq_records
        self.reverse = reverse

    def __iter__(self):
        return self

    def next(self):
        record = self.fastq_records.next()
        while not self.is_valid_record(record.id):
            record = self.fastq_records.next()

        trim_idx = self.trims[record.id]
        if self.reverse:
            return record[trim_idx[0]: trim_idx[1]]
        else:
            return record[trim_idx:]

    def is_valid_record(self, rid):
        '''does this record id have a corresponding trim entry?'''
        if self.reverse:
            return (rid in self.trims and all([idx is not None for idx in self.trims[rid]]))
        else:
            return (rid in self.trims)


class PrimerRemover():
    def __init__(self, trim_file, fastq, output, run=True):
        # determine if the trim file has reverse primers
        with open(trim_file) as f:
            n_fields = len(f.readline().split())

        if n_fields == 2:
            self.reverse = False
        elif n_fields == 5:
            self.reverse = True
        else:
            raise RuntimeError("trim file {} should have either 2 or 5 columns".format(trim_file))

        self.trim_records = open(trim_file)
        self.fastq_records = SeqIO.parse(fastq, 'fastq')
        self.output = output

        if run:
            self.run()

    def __del__(self):
        self.trim_records.close()
        self.fastq_records.close()

    def run(self):
        '''print the successfully trimmed entries'''

        # parse the trim file
        trims = self.parse_trim_file(self.trim_records, self.reverse)

        SeqIO.write(TrimmedRecords(trims, self.fastq_records, self.reverse), self.output, 'fastq')

    @staticmethod
    def parse_trim_file(trim_records, reverse):
        '''parse file produced by find_primers command'''

        trims = {}
        if reverse:
            for record in trim_records:
                rid, primer, strand, start_idx, end_idx = record.split()

                # check that primer and strand are in agreement
                # either fill in or create and entry [start index, end index]
                if primer == 'forward' and strand == '+':
                    if rid in trims:
                        trims[rid][0] = int(end_idx)
                    else:
                        trims[rid] = [int(end_idx), None]
                elif primer == 'reverse' and strand == '-':
                    if rid in trims:
                        trims[rid][1] = int(start_idx) - 1
                    else:
                        trims[rid] = [None, int(start_idx) - 1]
        else:
            # just add the start index
            for record in trim_records:
                rid, end_idx = record.split()
                trims[rid] = int(end_idx)

        return trims