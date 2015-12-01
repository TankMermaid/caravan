'''
Convert a raw fastq. This means you can
  - quality: Convert from old-school "Illumina" quality scores to Sanger scores
  - rename: Rename from read names like "@MISEQ:1:1101:15783:1330#NCTCCTCC/1"
    to names like "read1", "read2", etc.
  - index: Extract reads. Like from the above, get ">read1 \n NCTCCTCC"
'''

from Bio import SeqRecord, SeqIO, Seq
from Bio.Alphabet import IUPAC
import re

quality_table = str.maketrans({chr(x): chr(x - 31) for x in range(65, 105)})

def convert_fastq(fastq, output, quality=False, index=None, rename=False):
    if quality:
        entries = SeqIO.parse(fastq, 'fastq-illumina')
    else:
        entries = SeqIO.parse(fastq, 'fastq')

    def new_entries(entries):
        for entry_i, entry in enumerate(entries):
            if rename:
                new_id = "read{}".format(entry_i + 1)
            else:
                new_id = entry.id

            if index is not None:
                # try to find the sample
                try:
                    read = re.search('#([ACGTN]+)/', entry.id).groups()[0]
                    index.write(">{}\n{}\n".format(new_id, read))
                except:
                    raise RuntimeError("failed to parse entry {}".format(entry.id))

            if rename:
                entry.id = new_id

            # output the new sample
            entry.description = ""
            yield entry

    SeqIO.write(new_entries(entries), output, 'fastq')
