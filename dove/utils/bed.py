# -*- coding: utf-8 -*-
__author__ = 'bars'

from io import StringIO
import pandas as pd
from collections import defaultdict


class Bed:
    """description"""

    def __init__(self, bed_file, mode='file'):
        self.bed_file = bed_file
        self.mode = mode

    def get_header(self):
        lines_to_skip = 0
        header = defaultdict(list)
        if self.mode == 'str':
            for line in self.bed_file.split('\n'):
                if line.startswith('track'):
                    header['track'].append(line)
                    lines_to_skip += 1
                elif line.startswith('browser'):
                    header['browser'].append(line)
                    lines_to_skip += 1
                else:
                    break
        else:
            with open(self.bed_file) as f:
                lines = f.read().splitlines()
                for line in lines:
                    if line.startswith('track'):
                        header['track'].append(line)
                        lines_to_skip += 1
                    elif line.startswith('browser'):
                        header['browser'].append(line)
                        lines_to_skip += 1
                    else:
                        break
        return lines_to_skip, header

    def from_file(self):
        lines_to_skip, header = self.get_header()
        df_bed = pd.read_csv(
            self.bed_file,
            sep='\t',
            usecols=[0, 1, 2],
            names=['CHR', 'START', 'END'],
            dtype={'START': int, 'END': int},
            skiprows=lines_to_skip
        )
        return df_bed

    def from_string(self):
        lines_to_skip, header = self.get_header()
        df_bed = pd.read_csv(
            StringIO(self.bed_file),
            sep='\t',
            usecols=[0, 1, 2],
            names=['CHR', 'START', 'END'],
            dtype={'START': int, 'END': int},
            skiprows=lines_to_skip
        )
        return df_bed
