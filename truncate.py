'''
truncate a fastx at some position or quality
'''

from Bio import SeqRecord, SeqIO, Seq

def lengthed_entries(entries, length, keep):
    for entry in entries:
        if len(entry.seq) < length:
            if keep:
                yield entry
        else:
            yield entry[0: length]

def length(length, fastx, input_format, keep, output):
    entries = lengthed_entries(SeqIO.parse(fastx, input_formap), length, keep)
    SeqIO.write(entries, output, input_format)

def tailed_entries(entries, quality):
    for entry in entries:
        quals = entry.letter_annotations['phred_quality']
        if quals[-1] <= quality:
            i = 1
            while i <= len(quals) and quals[-i] <= quality:
                i += 1

            i -= 1 # since we got an extra increment at the end

            # if i == len(quals), then whole sequence was trash
            if i < len(quals):
                yield entry[0: -i]
        else:
            yield entry

def tail(quality, fastq, output):
    entries = tailed_entries(SeqIO.parse(fastq, 'fastq'), quality)
    SeqIO.write(entries, output, 'fastq')
