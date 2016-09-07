'''
Intersect multiple files

Note that this only works for 4-line fastq's. I ditched the Bio.SeqIO because it gets
confusing when you have multiple files open.
'''

def extract_number(header):
    m = re.match("read(\d+);?", header)
    return int(m.groups()[0])

def peek(input_fh, output_fh):
    first_line = next(input_fh)
    if first_line.startswith('@'):
        return FastqRecords(input_fh, output_fh, first_line)
    elif first_line.startswith('>'):
        return FastaRecords(input_fh, output_fh, first_line)
    elif first_line.count("\t") >= 1:
        return TsvRecords(input_fh, output_fh, first_line)
    else:
        raise RuntimeError("did not recognize file starting with: '{}'".format(first_line.rstrip()))

class Records:
    def __init__(self, input_fh, output_fh, first_line):
        self.input_fh = input_fh
        self.output_fh = output_fh
        self.first_line = first_line
        self.first = True


class FastqRecords(Records):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # read in three more lines
        self.record = [self.first_line] + [next(self.input_fh) for i in range(3)]

    def read_number(self):
        return extract_number(self.record[0].lstrip('@'))

    def increment(self):
        self.record = [next(self.input_fh) for i in range(4)]

    def output(self):
        for line in self.record:
            self.output_fh.write(line)


class FastaRecords(Records):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.record = [self.first_line, next(self.input_fh)]

    def read_number(self):
        return extract_number(self.record[0].lstrip('>'))

    def increment(self):
        self.record = [next(self.input_fh) for i in range(2)]

    def output(self):
        for line in self.record:
            self.output_fh.write(line)


class TsvRecords(Records):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.record = first_line

    def read_number(self):
        return extract_number(self.record.split("\t")[0])

    def increment(self):
        self.record = next(self.input_fh)

    def output(self):
        self.output_fh.write(self.record)


def intersect(inputs, outputs):
    if len(inputs) != len(outputs):
        raise RuntimeError("number of inputs and outputs not matched")

    keepers = [peek(inp, outp) for inp, outp in zip(inputs, outputs)]

    try:
        while True:
            current_nums = [k.read_number() for k in keepers]
            max_num = max(current_nums)

            to_skip = [i for i, num in enumerate(current_nums) if num < max_num]

            if len(to_skip) == 0:
                # all entries are aligned; output and update
                for k in keepers:
                    k.output()

                for k in keepers:
                    k.increment()
            else:
                # skip records that are the minimum ones
                for i in to_skip:
                    keepers[i].increment()

    except StopIteration:
        pass
