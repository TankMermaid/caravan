'''
class for parsing blast6 and uparse files
'''

import util, yaml

class Parser:
    @classmethod
    def usearch_to_yaml(cls, usearch, output):
        '''convert uparse or blast6 to yaml'''

        # parse the first line differently from the rest
        first_line = next(usearch)
        n_fields = len(first_line.split())

        if n_fields == 5:
            cls.uparse_to_yaml(first_line, usearch, output)
        elif n_fields == 12:
            cls.blast6_to_yaml(first_line, usearch, output)
        else:
            raise RuntimeError("mapping file with {} columsn not recognized as blast6 or uparse".format(n_fields))

    @staticmethod
    def uparse_to_yaml(first_line, usearch, output):
        write_out = lambda seq, otu: output.write("{}: {}\n".format(seq, otu))

        # for the first line, the 5th field is the otu name
        fields = first_line.rstrip().split('\t')
        seq = util.strip_fasta_label(fields[0])
        assert(fields[1] == 'otu')
        otu = fields[4]
        write_out(seq, otu)

        for line in usearch:
            fields = line.rstrip().split('\t')
            seq = util.strip_fasta_label(fields[0])
            hit_type = fields[1]

            if hit_type == 'otu':
                otu = fields[5]
                write_out(seq, otu)
            elif hit_type == 'match':
                otu = fields[4]
                write_out(seq, otu)

    @staticmethod
    def blast6_to_yaml(first_line, usearch, output, no_hit="*", save_no_hit=True):
        def parse_line(line):
            fields = line.rstrip().split('\t')
            query = util.strip_fasta_label(fields[0])
            target = fields[1]
            return query, target

        write_out = lambda query, target: output.write("{}: {}\n".format(query, target))

        def process_line(line):
            query, target = parse_line(line)
            if target != no_hit or save_no_hit:
                write_out(query, target)

        process_line(first_line)
        for line in usearch:
            process_line(line)
