# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 12:15:53 2016

@author: Lieneke Kusters

Compare all input sequences by constructing the trees and computing distance,
based on KL-Divergence
Note, reverse complement is not considered here
"""

import argparse
import numpy as np
from contexttree.sequences_load import seqgenerator
from contexttree.sequences_operate import seqdistance

parser = argparse.ArgumentParser(description='Load all data files (zipped ' +
                                 'fasta/fastq format), use FullTree to get ' +
                                 'trees and then calculate distances ' +
                                 'with respect to each other, based on KL-' +
                                 'Divergence and store it (txt format)')
parser.add_argument('-i', nargs='+', required=True,
                    help='input filenames (.fna.gz or .fastq.gz)')
parser.add_argument('-o', nargs=1, required=True,
                    help='output filename .txt (storing histogram of ' +
                    'distances)')
parser.add_argument('-d', nargs=1, type=int, required=True, dest='depth',
                    help='depth of tree')

args = parser.parse_args()
filenamesin = args.i
filenameout = args.o[0]
depth = args.depth[0]

sequences1 = seqgenerator(filenamesin)
sequences2 = seqgenerator(filenamesin)

distances = list()
for seq1 in sequences1:
    for seq2 in sequences2:
        distances.append(seqdistance(depth, seq1, seq2))

hist = np.histogram(distances, 100)  # only store the histogram
edges = hist[1]
values = hist[0]

f = open(filenameout, 'w')
f.write('Result from experiment estimating distances between sequences ' +
        'using full context trees and KL-Divergence calculation.\n')
f.write('tree depth = {0}\n'.format(depth))
f.write('bin edges = ')
for h in edges:
    f.write(str(h)+' ')
f.write('\nbin counts = ')
for h in values:
    f.write(str(h)+' ')
f.close()
