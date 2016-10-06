# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:55:56 2016

@author: Lieneke Kusters
"""

import numpy as np
from contexttree.sequences_load import seqgenerator
from contexttree.sequences_operate import seqdistance
from contexttree.TreeCounts import TreeCounts
from contexttree.FullTree import FullTree

FILENAME = 'human.22.rna.fna.gz'

# 1 test functionality of TreeCounts
# first use a short sequence to verify
depth = 2
seq = 'AAAAAAAAAAAA'  # 10 * A after context AA
tree = TreeCounts(depth, seq)
print(tree)
# now add another string
seq = 'CCCCCC'  # 4 * C after context CC
tree.updatesymbolcounts(seq)
print(tree)
# now make a copy and add the trees together
treecopy = tree.getcopy()
tree.combine(treecopy)
print(tree)

# 2 test functionality of TreeCounts when loading sequences
depth = 6
sequences = seqgenerator([FILENAME])
seq = next(sequences)
tree = TreeCounts(depth, seq)
print(tree)

# 3 test functionality of FullTree when loading sequences
depth = 2
seq = 'AAAAAAAAAAAA'  # 10 * A after context AA
tree = FullTree(depth, seq)
print(tree)
print(tree.getprobs())
print('Verify: ', np.log2((10+1/2)/(10+2)))
print('Verify: ', np.log2((1/2)/(10+2)))
# now add another string
seq = 'CCCCCC'  # 4 * C after context CC
tree.updatesymbolcounts(seq)
print(tree)
print('Rself: ', tree.getrself())
print(tree.getprobs())
seq = 'CCCCCC'  # 4 * C after context CC
tree2 = FullTree(depth, seq)
print('Divergence: ', [tree.getdivergence(tree2), tree2.getdivergence(tree)])
print('Distance: ', [tree.getdistance(tree2), tree2.getdistance(tree)])

# 4 test functionality of sequence compare
depth = 2
seq1 = 'AAAAAAAAAAAA'  # 10 * A after context AA
tree1 = FullTree(depth, seq1)
seq2 = 'CCCCCCCCCCAAAA'  # 10 * A after context AA
tree2 = FullTree(depth, seq2)
print([tree1.getdistance(tree2), seqdistance(depth, seq1, seq2)])

