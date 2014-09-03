from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import VariantCounterStrategy, GCRatioStrategy, NCounterStrategy

# Future: Ensure the file is valid.
def load_variant_files(g, paths_filename):
    """Load and process the paths file."""

    if paths_filename is not None:
        path_list = open(paths_filename)
        files = {}
        current_group = None
        for line in path_list:
            if line.startswith("#"):
                current_group = line[1:].strip()

                if current_group in g.groups:
                    raise Exception("[FAIL] Group %s has already been processed" % current_group)

                g.groups[current_group] = {}
                g.group_buckets[current_group] = {}
                g.group_counts[current_group] = []
                continue
            if current_group is not None:
                fields = line.split("\t")

                if fields[0] in g.files:
                    raise Exception("[FAIL] File %s has already been processed" % fields[0])

                g.files[fields[0]] = {
                    "path": fields[1].strip(),
                    "group": current_group
                }
        path_list.close()

def load_variants_from_file(g, path, group):
    """Load each variant position record from a Variant Query file into a given group."""

    f = open(path)
    for line in f:
        fields = line.strip().split("\t")

        chrno, pos = fields[0].split(":")
        chrno = int(chrno) # NOTE Explodes for allosomes
        pos = int(pos)     # NOTE Positions are 1-indexed

        # Check group exists
        if group not in g.groups:
            raise Exception()
        if chrno not in g.groups[group]:
            g.groups[group][chrno] = []

            if chrno not in g.chr_max_len:
                g.chr_max_len[chrno] = 1

        # NOTE No duplicate checking to prevent list lookups
        g.groups[group][chrno].append(pos)

        # Check whether this is the highest variant position seen on this chr
        if pos > g.chr_max_len[chrno]:
            g.chr_max_len[chrno] = pos

    f.close()

#TODO Methods may take a list of locations or may need to actually analyze
#     a proper genomic sequence
"""Execute Goldilocks search."""
g = Goldilocks(VariantCounterStrategy)
g.groups = {"ONE": {1: [1,2,5]}}
g.group_buckets = {"ONE": {}, "total": {}}
g.group_counts = {"ONE": [], "total": []}
g.chr_max_len = {1:12}
g.STRIDE = 1
g.LENGTH = 3
g.regions = g.search_regions()

candidates = g._filter("max", actual_distance=1)


#########################################
g = Goldilocks(GCRatioStrategy)

g.groups = {"ONE": {1: "_CCCGGGAGATTT"}}
g.group_buckets = {"ONE": {}, "total": {}}
g.group_counts = {"ONE": [], "total": []}
g.chr_max_len = {1:12}
g.STRIDE = 1
g.LENGTH = 3
g.regions = g.search_regions()

candidates = g._filter("max", actual_distance=1)


#########################################
g = Goldilocks(NCounterStrategy)

g.groups = {"ONE": {1: "_ANAGGGANACAN"}}
g.group_buckets = {"ONE": {}, "total": {}}
g.group_counts = {"ONE": [], "total": []}
g.chr_max_len = {1:12}
g.STRIDE = 1
g.LENGTH = 3
g.regions = g.search_regions()

candidates = g._filter("min", percentile_distance=75)

#for reg in g.regions:
#    print g.regions[reg]


#g.candidates = g.initial_filter()
#g.winners = g.enrich()

