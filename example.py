from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import VariantCounterStrategy, GCRatioStrategy, NCounterStrategy

#TODO Methods may take a list of locations or may need to actually analyze
#     a proper genomic sequence
"""Execute Goldilocks search."""
data = {"ONE": {1: [1,2,5]}}
g = Goldilocks(VariantCounterStrategy, data, is_seq=False, stride=1, length=3)

candidates = g._filter("max", actual_distance=1)


#########################################
data = {"ONE": {1: "_CCCGGGAGATTT"}}
g = Goldilocks(GCRatioStrategy, data, stride=1, length=3)

candidates = g._filter("max", actual_distance=1)


#########################################
data = {"ONE": {1: "_ANAGGGANACAN"}}
g = Goldilocks(NCounterStrategy, data, stride=1, length=3)

candidates = g._filter("min", percentile_distance=75)

#for reg in g.regions:
#    print g.regions[reg]


#g.candidates = g.initial_filter()
#g.winners = g.enrich()

