samples = [
        "my_sample",
        "my_other_sample"
]

seqs = {
    "one": [
        "CATCANCAT",
        "TATANTATA"
    ],
    2: [
        "NANANANANA",
        "GANGANGAN",
    ],
    "X": [
        "GATTACAGATTACAN",
        "GATTACAGATTACAN"
    ],
}


bases = {
        "A": 0,
        "G": 0,
        "T": 0,
        "C": 0,
        "N": 0,
}

"""
for i in range(0, len(seqs[0])-2, 1):
    for seq in seqs:
        triplet = seq[i:i+3]
        for base in triplet:
            if base in bases:
                bases[base] += 1

    print "%d: %s," % (i, bases)
    bases = {
            "A": 0,
            "G": 0,
            "T": 0,
            "C": 0,
            "N": 0,
    }
"""

for chrom in seqs:
    print "\"%s\": {" % chrom
    for seqi, seq in enumerate(seqs[chrom]):
        print "\t\"%s\": {" % samples[seqi]
        for i in range(0, len(seqs[chrom][0])-2, 1):
            triplet = seq[i:i+3]
            for base in triplet:
                if base in bases:
                    bases[base] += 1

            print "\t\t%d: %s," % (i, bases)
            bases = {
                    "A": 0,
                    "G": 0,
                    "T": 0,
                    "C": 0,
                    "N": 0,
            }
        print "\t},"
    print "},"
