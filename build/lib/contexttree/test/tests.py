# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:55:56 2016

@author: Lieneke Kusters
"""

FILENAME = 'human.22.rna.fna.gz'


from contexttree.loadsequences import seqgenerator
from contexttree import TreeCounts

# load example sequences
sequences = seqgenerator(list(FILENAME))