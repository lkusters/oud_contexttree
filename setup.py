# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:12:48 2016

@author: Lieneke Kusters
"""

from distutils.core import setup

setup(
    name='contexttree',
    version='0.2dev',
    description='contexttree package',
    author='Lieneke Kusters',
    packages=['contexttree', 'contexttree.test'],
    long_description=open('README.txt').read(),
    scripts=['bin/parse_seqdistances.py'],
    # modules=['TreeCounts', 'loadsequences'],
    # test_suite='contexttree.module.tests',
)
