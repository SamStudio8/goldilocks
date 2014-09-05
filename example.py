from goldilocks.goldilocks import Goldilocks
from goldilocks.strategies import VariantCounterStrategy, GCRatioStrategy, NCounterStrategy

#TODO Methods may take a list of locations or may need to actually analyze
#     a proper genomic sequence
"""Execute Goldilocks search."""
data = {"ONE": {1: [1,2,5]}}
g = Goldilocks(VariantCounterStrategy, data, is_seq=False, stride=1, length=3)

candidates = g._filter("max", actual_distance=1)
print("#WND\tVAL\tCHR\tPOSITIONS (INC.)")
for region in candidates:
    print("%d\t%.2f\t%s\t%10d - %10d" % (region["id"],
                                    region["group_counts"]["total"],
                                    region["chr"],
                                    region["pos_start"],
                                    region["pos_end"],
    ))


#########################################
data = {"ONE": {1: "_CCCGGGAGATTT"}}
g = Goldilocks(GCRatioStrategy, data, stride=1, length=3)

candidates = g._filter("max", actual_distance=1)
print("#WND\tVAL\tCHR\tPOSITIONS (INC.)")
for region in candidates:
    print("%d\t%.2f\t%s\t%10d - %10d" % (region["id"],
                                    region["group_counts"]["total"],
                                    region["chr"],
                                    region["pos_start"],
                                    region["pos_end"],
    ))


#########################################
data = {"ONE": {1: "_ANAGGGANACAN"}}
g = Goldilocks(NCounterStrategy, data, stride=1, length=3)

candidates = g._filter("min", percentile_distance=75, limit=0,
        exclusions={
            "start_lte": 3,
            "end_gte": 10
        })
print("#WND\tVAL\tCHR\tPOSITIONS (INC.)")
for region in candidates:
    print("%d\t%.2f\t%s\t%10d - %10d" % (region["id"],
                                    region["group_counts"]["total"],
                                    region["chr"],
                                    region["pos_start"],
                                    region["pos_end"],
    ))

#########################################
import sys
sys.exit()
from Mountain.IO import fastareaders
f = fastareaders.FASTA_Library("/store/sanger/ref/hs37d5.fa")
data = {"F": {
            1: f.get_next().seq,
        }}

# Will assume that files follow the recommendation that sequence lines
# are no longer than 80 characters
g = Goldilocks(GCRatioStrategy, data, stride=500000, length=1000000)

# Avoid human leukocyte antigen loci on chr6
candidates = g._filter("max", percentile_distance=10, limit=10, exclusions={"chr": [6]})

print("#WND\tVAL\tCHR\tPOSITIONS (INC.)")
for region in candidates:
    print("%d\t%.2f\t%s\t%10d - %10d" % (region["id"],
                                    region["group_counts"]["total"],
                                    region["chr"],
                                    region["pos_start"],
                                    region["pos_end"],
    ))
