'''
Intersect multiple files
'''

from Bio import SeqIO

class Peeker:
    def __init__(self, fh):
        if isinstance(fh, str):
            self.fh = open(fh)
        else:
            self.fh = fh

        self.first = True

        self.first_line = next(self.fh)

        if self.first_line.startswith('>'):
            self.filetype = 'fasta'
        elif self.first_line.startswith('@'):
            self.filetype = 'fastq'
        elif self.first_line.count("\t") == 1:
            self.filetype = 'two-column'
        else:
            raise RuntimeError("did not recognize file starting with: {}".format(self.first_line.rstrip()))

    def __iter__(self):
        return self

    def __next__(self):
        if self.first:
            self.first = False
            return self.first_line
        else:
            return next(self.fh)

    def readline(self):
        return next(self)

class TsvRecord:
    def __init__(self, line):
        self.line = line
        self.id = line.split("\t")[0]

    @staticmethod
    def tsv_records(lines):
        for line in lines:
            yield TsvRecord(line)

class RecordKeeper:
    def __init__(self, fh, output):
        if isinstance(fh, str):
            self.fh = open(fh)
        else:
            self.fh = fh

        self.output = output

        self.peeker = Peeker(self.fh)

        if self.peeker.filetype in ['fasta', 'fastq']:
            self.records = SeqIO.parse(self.peeker, self.peeker.filetype)
            self.output_record = lambda: SeqIO.write(self.record, self.output, self.peeker.filetype)
        elif self.peeker.filetype == 'two-column':
            self.records = TsvRecord.tsv_records(self.peeker)
            self.output_record = lambda: self.output.write(self.record.line + "\n")
        else:
            raise RuntimeError("don't recognize filetype: {}".format(self.peeker.filetype))

        self.skip_record()

    def skip_record(self):
        print("  obj skip")
        self.record = next(self.records)
        print("  done skip")

    def current_id(self):
        return self.record.id


def extract_number(header):
    return int(header.lstrip('read'))

def intersect(inputs, outputs):
    keepers = [RecordKeeper(inp, outp) for inp, outp in zip(inputs, outputs)]

    while True:
        print("start while")
        current_is = [extract_number(keeper.current_id()) for keeper in keepers]
        print(current_is)
        max_i = max(current_is)
        
        to_skip = [keeper for i, keeper in zip(current_is, keepers) if i < max_i]
        print("skipping", [i for i in current_is if i < max_i])

        if len(to_skip) == 0:
            # all entries are aligned; output and update
            print("outputting all")
            for keeper in keepers:
                print("  output")
                keeper.output_record()
                keeper.skip_record()
        else:
            # skip records that are the minimum ones
            print("skipping some")
            for keepr in to_skip:
                print("  skip")
                keeper.skip_record()

            print("done skipping")
