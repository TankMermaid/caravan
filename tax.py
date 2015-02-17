'''
use a blast6 alignment file to get taxonomies
'''

import cPickle as pickle

class TaxAssigner:
    @staticmethod
    def assign_b6_with_pickled_tax_dict(b6, db, output, no_hit_target='*', no_hit_tax='no_hit'):
        '''
        look up each reference id in a pickled database of names

        b6 : filename
            blast6 mapping from usearch
        db : filename
            pickled dictionary gg_id => tax
        output : filehandle
            output taxonomy list
        '''

        raise RuntimeError("tax assignment not implemented")

        with open(db) as f:
            taxes = pickle.load(f)

        taxes[no_hit_target] = no_hit_tax

        with open(b6) as f:
            targets = [o['target_label'] for o in b6m.B6.parse_lines(f)]

        outs = [taxes[target] for target in targets]
        output.write("\n".join(outs) + "\n")