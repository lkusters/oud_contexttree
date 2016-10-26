# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 14:48:24 2016

@author: lkusters
"""

import argparse
from sequences_operate import seqdistance
from sequences_load import seqgenerator
import numpy as np

parser = argparse.ArgumentParser(
                    description="Load all data files (zipped fasta/fastq "
                                "format), use FullTree to estimate the "
                                "similarity (KL-Divergence) and store in "
                                "histogram "
                                "(txt format)"
                                )
parser.add_argument('-i', nargs='+', required=True,
                    help='input filenames (.fna.gz or .fastq.gz)')
parser.add_argument('-o', nargs=1, required=True,
                    help='output filename (storing the entropies)')
parser.add_argument('-d', nargs=1, type=int, required=True, dest='dmax',
                    help='value of tree depth')
parser.add_argument('-r', nargs=1, type=str, required=True, dest='revcomp',
                    help="automatically also include reverse complement "
                    "of each sequence? y/n")

args = parser.parse_args()
filenamesin = args.i
filenameout = args.o[0]
depth = args.dmax[0]
if args.revcomp[0] == 'y':
    revcom = True
elif args.revcomp[0] == 'n':
    revcom = False
else:
    raise ValueError("not clear wether reverse complement should be "
                     "included, choose -r y/n", args.revcomp[0])

sequences = seqgenerator(filenamesin)
distances = list()
for seq1 in sequences:
    for seq2 in sequences:
        distances.append(seqdistance(depth, seq1, seq2, revcom))

hist = np.histogram(distances, 100)  # only store the histogram of H
edges = hist[1]
values = hist[0]
# print(entropies)
f = open(filenameout, 'w')
f.write("Result from experiment estimating similarity using full trees of "
        "{0} depth with KL-Divergence, reverse complement is included? {1}\n"
        .format(depth, args.revcomp[0])
        )
f.write('depth = {0}\n'.format(depth))
f.write('bin edges = ')
for h in edges:
    f.write(str(h)+' ')
f.write('bin counts = ')
for h in values:
    f.write(str(h)+' ')
f.close()
