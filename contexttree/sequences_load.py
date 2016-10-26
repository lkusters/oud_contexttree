# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 11:03:06 2016

@author: Lieneke Kusters
"""

import gzip
from Bio import SeqIO


# Functions that can be used to load sequences from a file
def seqgenerator(filenames_list):
    """ This function will load sequences from a zipped fasta or fastq file
    and return a generator of sequences (without identifiers)
    """
    for filename in filenames_list:
        handle = gzip.open(filename, 'rt')
        checkextension = filename.split('.')
        if checkextension[-2] == 'fna' or checkextension[-2] == 'fa':
            for record in SeqIO.parse(handle, 'fasta'):
                yield str(record.seq)
        elif checkextension[-2] == 'fastq':
            for record in SeqIO.parse(handle, 'fastq'):
                yield str(record.seq)
        else:
            print("filename extension {0} not recognised".format(
                                                      checkextension[-2]))
            continue
