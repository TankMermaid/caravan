'''
Quality filter reads based on expected errors
'''

from Bio import SeqIO

def ee(entry):
    return sum([10 ** (-q / 10) for q in entry.letter_annotations['phred_quality']])

def qfilter(fastq, maxee, output):
    entries = SeqIO.parse(fastq, 'fastq')
    criterion = lambda entry: ee(entry) <= maxee

    filtered_entries = filter(criterion, entries)
    SeqIO.write(filtered_entries, output, 'fastq')