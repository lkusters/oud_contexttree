# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:55:56 2016

@author: Lieneke Kusters
"""

FILENAME = 'human.22.rna.fna.gz'

from contexttree.loadsequences import seqgenerator
from contexttree.TreeCounts import TreeCounts

# first use a short sequence to verify
depth = 2
seq = 'AAAAAAAAAAAA' # 10 * A after context AA
tree = TreeCounts(depth,seq)
print(tree)
# now add another string
seq = 'CCCCCC' # 4 * C after context CC
tree.updatesymbolcounts(seq)
print(tree)
# now make a copy and add the trees together
treecopy = tree.getcopy()
tree.combine(treecopy)
print(tree)

# now load example sequences and construct corresponding models
depth = 6
sequences = seqgenerator([FILENAME])
seq = next(sequences)
tree = TreeCounts(depth,seq)
print(tree)
