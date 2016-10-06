# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:52:50 2016

@author: Lieneke Kusters

functions that use FullTree class to operate on DNA sequences
"""

from contexttree.FullTree import FullTree


def seqdistance(depth, seq1, seq2):
    """
    Compare 2 DNA sequences, by constructing the trees and computing distance,
    based on KL-Divergence
    """

    tree1 = FullTree(depth, seq1)
    tree2 = FullTree(depth, seq2)
    return tree1.getdistance(tree2)
