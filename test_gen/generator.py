import math

sequence = "AAACCCNNNTTTAAACCCGGGNNNAAACCCGGGTTTNNNCCCGGGTTT"
sequence = "GCGTNANNGANGGCTANTCTANAGCNNTTTCTNTNNGCANCANTTGNN"
lenseq = len(sequence)

def mark(seq):
    markseq = ["="] * len(seq)
    for i, base in enumerate(seq):
        if base=="N":
            markseq[i] = "*"
    return "".join(markseq)

def entry_contains(subseq):
    count = 0
    for base in subseq:
        if base=="N":
            count += 1
    return count

print
print "#       1 |%s| %d" % (mark(sequence), lenseq)
print "#          %s" % sequence

#TODO Length must be at least 3
stride = 1
length = 3

counts = []
entries = []
entry_size = 0
for i, region_s in enumerate(range(0, lenseq-length+1, stride)):
    region_e = region_s + length - 1

    count = entry_contains(sequence[region_s:region_e+1])
    counts.append(count)
    if length % 2 == 0:
        bucket_str = ("|%s%d%s|" % (" "*int(math.floor((length-2)/2)), count, " "*int(math.floor((length-2)/2)-1)))
    else:
        bucket_str = ("|%s%d%s|" % (" "*int(math.floor((length-2)/2)), count, " "*int(math.floor((length-2)/2))))
    entries.append(bucket_str)
entry_size = len(bucket_str)+stride

# Calculate how many entries can fit on each line beneath the sequence
# when taking the indentation in to consideration
allocations = []
allocated = 0
indent = 0
while allocated < len(entries):
    allocatable = ( int(math.floor((lenseq-indent+stride)/entry_size)) )
    allocations.append(allocatable)
    allocated += allocatable
    indent += 1

# Build the lines
lines = []
for i in range(0, len(allocations)):
    line_i = ""
    for j in range(0, allocations[i]):
        stride_gap = ""
        if j > 0:
            stride_gap = " "*(stride)
        try:
            to_append = stride_gap + entries[ ((j*len(allocations)) + i) ]
            line_i += to_append
        except IndexError:
            # No data left to add
            break
    lines.append(line_i)

for i, l in enumerate(lines):
    print "#%s          %s" % (" "*((i*stride)), "".join(l))

import numpy as np
#for i in reversed(np.argsort(counts)):
#    print "%d\t%d" % (i, counts[i])

buckets = {}
print "EXPECTED_REGIONS = {"
for i, count in enumerate(counts):
    if count not in buckets:
        buckets[count] = []
    buckets[count].append(i)
    print "    %d: %d," % (i, count)
print "}"

print "EXPECTED_BUCKETS = {"
for b in buckets:
    print "    %d: %s," % (b, buckets[b])
print "}"
