# -*- coding: utf-8 -*-
__author__ = 'bars'

import os
from dove import examples
from shutil import copyfile


def main(args):
    data_path = examples.__path__[0]
    cwd = os.getcwd()
    examplefile = args.examplefile

    pedigree = '/examplepedigree.csv'
    filterfile = '/examplefilter.txt'

    if examplefile == 'both':
        copyfile(data_path + pedigree, cwd + pedigree)
        copyfile(data_path + filterfile, cwd + filterfile)
    if examplefile == 'pedigree':
        copyfile(data_path + pedigree, cwd + pedigree)
    if examplefile == 'filterfile':
        copyfile(data_path + filterfile, cwd + filterfile)
