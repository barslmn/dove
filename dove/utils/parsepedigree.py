# -*- coding: utf-8 -*-
__author__ = 'bars'

import pandas as pd
import sys


class ParsePedigree():

    """Docstring for ParsePedigree. """

    def __init__(self, pedigree):
        self.pedigree = pedigree

    def parsepedigree(self):
        df = pd.read_csv(self.pedigree)
        if df['index'].sum() > 1:
            sys.exit('There might be more than one index. Exiting')
        if df.isnull().values.any() == True:
            sys.exit('There are missing Values. Exiting')

        return df.to_dict(orient='records')
