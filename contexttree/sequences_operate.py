# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:52:50 2016

@author: Lieneke Kusters

functions that use FullTree class to operate on DNA sequences
"""

from contexttree.FullTree import FullTree
from Bio.Seq import Seq


def seqdistance(depth, seq1, seq2, rev):
    """
    Compare 2 DNA sequences, by constructing the trees and computing distance,
    based on KL-Divergence
    rev = True means that reverse complement of the sequences should be
    included in the models
    """

    tree1 = FullTree(depth, seq1)
    tree2 = FullTree(depth, seq2)
    if rev:
        tree1.updatesymbolcounts(depth, str(Seq(seq1).reverse_complement()))
        tree2.updatesymbolcounts(depth, str(Seq(seq1).reverse_complement()))

    return tree1.getdistance(tree2)
