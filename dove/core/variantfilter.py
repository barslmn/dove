# -*- coding: utf-8 -*-
__author__ = 'bars'

import pandas as pd
from dove.utils.vcf import Vcf


class VariantFilter:

    """Docstring for VariantFilter. """

    def __init__(self, df, file_type, filter_columns, drop_columns=None, bed_file=None):
        self.df = df
        self.file_type = file_type
        self.filter_columns = filter_columns
        self.drop_columns = drop_columns
        self.bed_file = bed_file

    def filter_variants(self):
        if self.bed_file is not None:
            self.df = self.filter_with_bed()

        if self.filter_columns is not None:
            list_of_filter_dicts = []
            for filter_column in self.filter_columns:
                filter_dict = {'column': filter_column[0],
                               'filter_option': filter_column[1].lower(),
                               'filter_args': filter_column[2:]}
                list_of_filter_dicts.append(filter_dict)

            for filter_dict in list_of_filter_dicts:
                if filter_dict['filter_option'] in self.in_ex_options('both'):
                    self.df = self.includes_excludes(
                        filter_dict['column'], filter_dict['filter_option'], filter_dict['filter_args'])
                if filter_dict['filter_option'] in self.eq_ne_options('both'):
                    self.df = self.equals_notequals(
                        filter_dict['column'], filter_dict['filter_option'], filter_dict['filter_args'])
                if filter_dict['filter_option'] in self.l_g_options('all'):
                    self.df = self.less_greater(
                        filter_dict['column'], filter_dict['filter_option'], filter_dict['filter_args'])

        if self.drop_columns is not None:
            self.df.drop_duplicates(self.drop_columns, inplace=True)

        return self.df

    def filter_with_bed(self):
        df_bed = pd.read_table(self.bed_file,
                               usecols=[0, 1, 2],
                               delim_whitespace=True,
                               names=['CHR', 'START', 'END'],
                               dtype={'START': int, 'END': int})

        if self.file_type == 'annotation':
            self.df['POS'] = self.df['LOC'].str.split('-').str[0]
            self.df['POS'] = self.df['POS'].apply(pd.to_numeric)

        idx = pd.IntervalIndex.from_arrays(
            df_bed['START'], df_bed['END'], closed='both')
        df_bed.set_index(idx, inplace=True)
        if self.file_type == 'annotation':
            mask = self.df.apply(lambda x: [
                                 x['POS'] in y for y in df_bed.loc[df_bed.CHR == x.CHR, ].index], axis=1)
        if self.file_type == 'vcf':
            mask = self.df.apply(lambda x: [
                                 x['POS'] in y for y in df_bed.loc[df_bed.CHR == x.CHROM, ].index], axis=1)
        mask = mask.apply(lambda x: sum(x)) > 0
        return self.df[mask]

    def equals_notequals(self, column, filter_option, filter_args):
        df_filter = pd.DataFrame(filter_args, columns=[column])
        if filter_option in self.eq_ne_options('eq'):
            return self.df[self.df[column].isin(df_filter[column])]
        if filter_option in self.eq_ne_options('ne'):
            return self.df[~self.df[column].isin(df_filter[column])]

    def less_greater(self, column, filter_option, filter_args):
        if filter_option in self.l_g_options('lt'):
            return self.df[pd.to_numeric(self.df[column], downcast='float') < float(filter_args[0])]
        if filter_option in self.l_g_options('le'):
            return self.df[pd.to_numeric(self.df[column], downcast='float') <= float(filter_args[0])]
        if filter_option in self.l_g_options('gt'):
            return self.df[pd.to_numeric(self.df[column], downcast='float') > float(filter_args[0])]
        if filter_option in self.l_g_options('ge'):
            return self.df[pd.to_numeric(self.df[column], downcast='float') >= float(filter_args[0])]

    def includes_excludes(self, column, filter_option, filter_args):
        if filter_option in self.in_ex_options('in'):
            filter_option = True
        if filter_option in self.in_ex_options('ex'):
            filter_option = False

        return self.df[self.df[column].str.contains('|'.join(filter_args), case=False, na=False) == filter_option]

    def eq_ne_options(self, option):
        eqs = ['eq', 'equal', 'equals']
        nes = ['ne', 'notequals', 'not_equals']
        if option == 'both':
            return eqs + nes
        if option == 'eq':
            return eqs
        if option == 'ne':
            return nes

    def l_g_options(self, option):
        lts = ['lt', 'lessthan', 'less_than']
        les = ['le', 'lessequal', 'less_equals',
               'lessthanorequals', 'less_than_or_equals']
        gts = ['gt', 'greaterthan', 'greater_than']
        ges = ['ge', 'greaterequals', 'greater_equals',
               'greaterthanorequals', 'greater_than_or_equals']
        if option == 'all':
            return lts + les + gts + ges
        if option == 'l':
            return lts + les
        if option == 'g':
            return gts + ges
        if option == 'lt':
            return lts
        if option == 'le':
            return les
        if option == 'gt':
            return gts
        if option == 'ge':
            return ges

    def in_ex_options(self, option):
        ins = ['in', 'include', 'includes']
        exs = ['ex', 'exclude', 'excludes']
        if option == 'both':
            return ins + exs
        if option == 'in':
            return ins
        if option == 'ex':
            return exs

    def check_input(self):
        pass


def main(args):
    input_file = args.input
    output_file = args.output
    bed_file = args.bed_file
    filter_columns = args.column
    drop_columns = args.drop

    if input_file.endswith('.csv'):
        df = pd.read_csv(input_file, low_memory=False)
        VF = VariantFilter(df, 'annotation', filter_columns,
                           drop_columns, bed_file)
        df = VF.filter_variants()
        df.to_csv(output_file, index=False)
    if input_file.endswith(('.vcf', '.vcf.gz')):
        with Vcf(input_file) as vcf:
            VF = VariantFilter(vcf.vdf, 'vcf', filter_columns,
                               drop_columns, bed_file)
            vcf.vdf = VF.filter_variants()
            vcf.to_vcf(output_file)
