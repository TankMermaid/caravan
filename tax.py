'''
use a blast6 alignment file to get taxonomies
'''

import util

class TaxAssigner:
    @staticmethod
    def b6_line_to_query_hit(line):
        query, hit = line.rstrip().split("\t")[0: 2]
        query = util.strip_fasta_label(query)
        return query, hit

    @classmethod
    def assign_b6_tax(cls, b6, db, output, no_hit_target='*', no_hit_tax='no_hit'):
        '''
        look up each reference id in a pickled database of names

        b6 : filename
            blast6 mapping from usearch
        db : filename
            tab-separated file {id => taxonomy}, like the Greengenes 97_otu_taxonomy.txt
        output : filehandle
            membership yaml
        '''

        # read the taxonomic database file
        tax = dict([line.rstrip().split("\t") for line in db])
        tax[no_hit_target] = no_hit_tax

        # read the b6, extracting the queries (my OTU IDs) and hits (database OTU IDs)
        for line in b6:
            query, hit = cls.b6_line_to_query_hit(line)

            try:
                query_tax = tax[hit]
            except KeyError:
                raise RuntimeError("could not find taxonomy for '{}' in database".format(hit))

            # eliminate spaces
            query_tax = query_tax.translate({ord(' '): None})

            # if the query is an integer, write it unquoted
            try:
                query = int(query)
                fmt = '{}: {}'
            except ValueError:
                fmt = '{}: "{}"'

            print(fmt.format(query, query_tax), file=output)
