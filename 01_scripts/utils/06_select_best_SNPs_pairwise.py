#!/usr/bin/env python3
"""Choose the best SNPs for all pairs of populations

For each pair of populations
    - sort the AFDs by decreasing order
    - Pick the top ones (exclude them if they are within `min_distance` of a previously retained SNPs)
    - Go down the list until you have reached `target_sum` of AFDs <float>

Continue for following populations and report retained SNPs

Usage:
    <program> input_afds target_sum exponent min_distance output_snps

Format of input_afds file (tab-separated):

ChromName    Position  AFDScore  AFDScore  AFDScore  AFDScore  AFDScore  AFDScore  AFDScore
NC_036838.1  1743435   0.171305  0.272733  0.115156  0.056149  0.387889  0.444038  0.014845
NC_036838.1  1760446   0.182961  0.384065  0.077204  0.105757  0.306861  0.201104  0.350584
NC_036838.1  1902438   0.03789   0.048098  0.122343  0.160233  0.170441  0.010208  0.319688
NC_036838.1  1926887   0.187807  0.031834  0.06648   0.254287  0.034646  0.219641  0.386849
NC_036838.1  1942787   0.155984  0.147081  0.02978   0.185764  0.117301  0.303065  0.263861
NC_036838.1  1949218   0.200558  0.065416  0.034061  0.166497  0.099477  0.265974  0.368012
NC_036838.1  2009390   0.106727  0.03674   0.379363  0.272636  0.342623  0.069987  0.205130
NC_036838.1  2171262   0.194908  0.427275  0.248636  0.053728  0.178639  0.232367  0.203546
NC_036838.1  3253985   0.144088  0.479436  0.307792  0.163704  0.171644  0.335348  0.023561
"""

# Modules
from collections import defaultdict
import pandas as pd
import sys

# Parse user input
try:
    input_afds = sys.argv[1]
    target_sum = float(sys.argv[2])
    exponent = float(sys.argv[3])
    min_distance = int(sys.argv[4])
    output_snps = sys.argv[5]
except:
    print(__doc__)
    sys.exit(1)

# Set variables
snp_positions = defaultdict(list)
snp_positions["fake"].append(-min_distance)
retained_snps = list()
retained = 0

df = pd.read_csv(input_afds, sep="\t")

# Iterate over columns
for col in list(df.columns):
    num_retained_col = 0
    current_sum = 0.0

    if "AFD" in col:
        # Sort df by that column
        df.sort_values(col, ascending=False, inplace=True)
        df[col] = df[col] ** exponent

        values = list(df[col])
        chroms = list(df["ChromName"])
        positions = list(df["Position"])
        afds = list(zip(values, chroms, positions))

        for a in afds:
            afd, chrom, pos = a

            # Check if too close to another SNP
            is_too_close = len([x for x in snp_positions[chrom] if abs(x-pos) < min_distance])

            if not is_too_close:
                snp_positions[chrom].append(pos)
                retained_snps.append((chrom, pos))
                current_sum += afd
                num_retained_col += 1

            if current_sum > target_sum:
                print(f"  {col}:\tRetained {num_retained_col} SNPs for {col}")
                break

print(f"Retained a total of {len(retained_snps)}")

# Output results
with open(output_snps, "wt") as outfile:
    for snp in sorted(retained_snps):
        outfile.write(snp[0] + "\t" + str(snp[1]) + "\n")
