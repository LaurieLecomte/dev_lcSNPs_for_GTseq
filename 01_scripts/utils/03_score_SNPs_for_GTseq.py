#!/usr/bin/env python3
"""For each SNP of interest, extract information about flanking SNPs, sequence
complexity and GC content, etc.

Usage:
    <program> input_selected_snps input_all_snps input_genome window_size output_file

`input_selected_snps` file format:

NC_036838.1	60649
NC_036838.1	84727
NC_036838.1	434621
NC_036838.1	981627
NC_036838.1	1702986
NC_036838.1	1758963
NC_036838.1	1761652

`input_all_snps` format:

chromo       position  major  minor  anc  knownEM   nInd
NC_036838.1  25        G      A      G    0.010773  138
NC_036838.1  29        G      T      G    0.844871  139
NC_036838.1  53        A      C      A    0.036794  156
NC_036838.1  67        C      T      C    0.011812  174
NC_036838.1  86        G      A      G    0.161940  178
NC_036838.1  88        T      G      T    0.022172  174

`input_genome` format: fasta or gzip-compressed fasta

`window_size`: Minimum distance between two retained SNPs <int>

`output_file`: Name of output file
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[:31]

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with myopen(input_file) as f:
        sequence = []
        name = ""
        begun = False

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if begun:
                    yield Fasta(name, "".join(sequence))

                name = line[1:]
                sequence = ""
                begun = True

            else:
                sequence += line

        if name != "":
            yield Fasta(name, "".join(sequence))

# Parsing user input
try:
    input_selected_snps = sys.argv[1]
    input_all_snps = sys.argv[2]
    input_genome = sys.argv[3]
    window_size = int(sys.argv[4])
    output_file = sys.argv[5]
except:
    print(__doc__)
    sys.exit(1)

# Load genome
print("Loading genome sequences")
genome = dict()

fasta = fasta_iterator(input_genome)

for f in fasta:
    # Report only long scaffolds
    if len(f.sequence) > 100000:
        print("    Scaffold: " + f.name.split(" ")[0])

    genome[f.name.split(" ")[0]] = f.sequence.upper()

# Load all SNPs
all_snps = defaultdict(lambda: defaultdict(list))

print("Loading SNPs")
with myopen(input_all_snps) as infile:
    for line in infile:
        if line.startswith("chromo"):
            continue

        l = line.strip().split("\t")
        chrom, pos, major, minor, ancestral, EM, nind = l
        pos = int(pos)
        EM = float(EM)
        nind = int(nind)

        pos_truncated = pos // window_size

        all_snps[chrom][pos_truncated].append((chrom, pos, major, minor, ancestral, EM, nind))

# Score selected SNPs using all the SNPs
print("Scoring SNPs")
seen_chrom = set()
with open(input_selected_snps) as infile:

    with open(output_file, "wt") as outfile:
        # Header
        outfile.write("Chrom\tPosition\tNumSNPs\tSumMAFs\tComplexity\tGCcontent\tChrom2\tPosition2\tMajor\tMinor\tAncestral\tFst\tNumSamples\tSequenceMod\tSequence\n")

        for line in infile:
            if "position" in line:
                continue

            l = line.strip().split("\t")
            chrom = l[0]
            if not chrom in seen_chrom:
                print(f"    Scaffold: {chrom}")
                seen_chrom.add(chrom)

            pos = int(l[1])

            pos_truncated = pos // window_size

            flanking_snps = (all_snps[chrom][pos_truncated - 1] +
                    all_snps[chrom][pos_truncated] +
                    all_snps[chrom][pos_truncated + 1])

            center_snp = [x for x in flanking_snps if abs(x[1] - pos) == 0][0]

            flanking_snps = [x for x in flanking_snps if
                    abs(x[1] - pos) and abs(x[1] - pos) <= window_size]

            sum_mafs = sum([x[5] for x in flanking_snps])
            num_snps = len(flanking_snps)

            left = max(pos - window_size - 1, 0)
            right = min(pos + window_size, len(genome[chrom]))
            
            seq = list(genome[chrom][left: right])
            seq_original = seq[:]

            # Use compression, patterns, GC content...
            complexity = len(gzip.compress("".join(seq).encode()))
            gc_content = (seq.count("G") + seq.count("C")) / (len(seq) - seq.count("N"))

            for s in flanking_snps:
                seq[s[1] - pos + window_size] = "N"

            seq[window_size] = "[" + center_snp[2] + "/" + center_snp[3] + "]"

            outfile.write("\t".join([str(x) for x in l +
                [num_snps, round(sum_mafs, 4), complexity, round(gc_content, 4)] +
                list(center_snp) +
                ["".join(seq), seq_original]
                ]) + "\n")
