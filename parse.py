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
        fields = first_line.rstrip().split()
        n_fields = len(fields)

        if n_fields == 5:
            cls.uparse_to_yaml(first_line, usearch, output)
        elif n_fields == 12:
            cls.blast6_to_yaml(first_line, usearch, output)
        else:
            raise RuntimeError("mapping file with {} columns not recognized as blast6 or uparse: {}".format(n_fields, fields))

    @staticmethod
    def uparse_to_yaml(first_line, usearch, output):
        def process_line(line):
            fields = line.rstrip().split('\t')
            seq = util.strip_fasta_label(fields[0])
            hit_type = fields[1]

            if hit_type == 'otu' and fields[2] == '*':
                otu = fields[4]
            elif hit_type == 'otu' and fields[2] != '*':
                otu = fields[5]
            elif hit_type == 'match':
                otu = fields[4]

            if hit_type != 'chimera':
                output.write("{}: {}\n".format(seq, otu))

        process_line(first_line)

        for line in usearch:
            process_line(line)

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
