'''
class for parsing blast6 files
'''

import util, yaml

class Parser:
    @staticmethod
    def b6_to_dict(b6, no_hit="*", save_no_hit=True):
        '''from a blast6 mapping file, make a dict {query => target}'''

        membership = {}
        with open(b6) as f:
            for line in f:
                fields = line.split()
                query = util.strip_fasta_label(fields[0])
                target = fields[1]

                if target != no_hit or save_no_hit:
                    membership[query] = target

        return membership

    @staticmethod
    def up_to_dict(up, keep_chimera=False):
        '''from a uparse mapping file, make a dict {seq => otu}'''

        membership = {}
        with open(up) as f:
            for line in f:
                fields = line.split()
                seq = util.strip_fasta_label(fields[0])
                hit_type = fields[1]
                otu = fields[4]

                if hit_type == 'chimera':
                    if keep_chimera:
                        otu = 'chimera'
                    else:
                        # ignore this entry
                        continue

                membership[seq] = otu

        return membership

    @classmethod
    def map_to_yaml(cls, map_fn, yaml_fn):
        with open(map_fn) as f:
            line = f.readline()

        # count the number of fields
        n_fields = len(line.split())

        if n_fields == 5:
            d = cls.up_to_dict(map_fn)
        elif n_fields == 12:
            d = cls.b6_to_dict(map_fn)
        else:
            raise RuntimeError("mapping file with {} columsn not recognized as blast6 or uparse".format(n_fields))

        with open(yml_fn, 'w') as f:
            yaml.dump(d, f, default_flow_style=False)
