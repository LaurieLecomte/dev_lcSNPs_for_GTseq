#!/usr/bin/env python3
""" From the computed pairwise AFD, keep only the lines for the SNPs that
passed the filters

Usage:
    <program> input_AFDs input_retained_SNPs output_AFDs
"""

# Modules
import sys

# Parsing user input
try:
    input_AFDs = sys.argv[1]
    input_retained_SNPs = sys.argv[2]
    output_AFDs = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read info about retained SNPs
retained = set([tuple(x.strip().split("\t")[:2]) for x in open(input_retained_SNPs).readlines()])

# Keep only retained SNPs
header_written = False

with open(input_AFDs) as infile:
    with open(output_AFDs, "wt") as outfile:
        for line in infile:
            l = line.strip().split("\t")

            if not header_written:
                outfile.write("ChromName\tPosition\t" + "\t".join(["AFD"] * (len(l) - 2)) + "\n")
                header_written = True

            if tuple(l[:2]) in retained:
                outfile.write(line)
