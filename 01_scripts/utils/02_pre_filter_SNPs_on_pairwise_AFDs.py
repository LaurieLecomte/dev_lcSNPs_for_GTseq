#!/usr/bin/env python3
"""Report SNPs whith a maximum pairwise AFD value above user theshold

Usage:
    <program> input_afds min_afd output_ids

Examples input format (number of columns dependents on the number of pairwise tests):

ChromName	pos	afd1	afd2	afd3
Chr1		33	0.07	0.14	0.30
Chr2		98	0.61	0.22	0.18
"""

# Modules
import sys

# Parse user input
try:
    input_afds = sys.argv[1]
    min_afd = float(sys.argv[2])
    output_ids = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read and filter
with open(output_ids, "wt") as outfile:
    with open(input_afds) as infile:
        for line in infile:
            if line.startswith("ChromName"):
                continue

            l = line.strip().split("\t")
            infos = l[:2]
            data = [round(float(x), 2) for x in l[2:]]
            if max(data) >= min_afd:
                outfile.write("\t".join(infos) + "\n")
