# -*- coding: utf-8 -*-
__author__ = 'bars'

from collections import OrderedDict
from collections import defaultdict
from io import StringIO
import gzip
import pandas as pd
import sys


class Vcf():
    '''
    Main class for parsing vcf and its columns like sample format and variant info.
    '''

    def __init__(self, vcf_file, mode='file', compression='infer', chunksize=None):
        self.vcf_file = vcf_file
        self.mode = mode
        self.compression = compression
        # Chunksize not used right now.
        self.chunksize = chunksize
        # Generic VCF header without sample_names. For getting the sample_names.
        self.generic_header = ['CHROM',
                               'POS',
                               'ID',
                               'REF',
                               'ALT',
                               'QUAL',
                               'FILTER',
                               'INFO',
                               'FORMAT']

    def __enter__(self):
        if self.mode == 'file':
            self.get_header_n_meta_info()
            self.from_file()
            return self
        if self.mode == 'str':
            self.get_header_n_meta_info()
            self.from_string()
            return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    def __str__(self):
        '''sample_names and number of variants'''
        return '{} {}'.format(' '.join(self.sample_names), len(self.vdf))

    def __repr__(self):
        '''sample_names and number of variants'''
        return '{} {} {}'.format(self.__class__.__name__, ' '.join(self.sample_names), len(self.vdf))

    def __len__(self):
        '''Number of variants'''
        return len(self.vdf)

    @property
    def sample_names(self):
        '''sample_names is header - generic_header'''
        return [colname for colname in self.header if colname not in self.generic_header]

    @property
    def vdf(self):
        '''These functions for setting vdf attribute.'''
        return self.variant_df

    @vdf.setter
    def vdf(self, new_df):
        '''These functions for setting vdf attribute.'''
        self.variant_df = new_df

    def get_header_n_meta_info(self):
        '''Get meta_info and header from vcf file.
           This function returns a dictionary.'''
        if self.mode == 'file':
            # Open file
            if self.vcf_file.endswith('.vcf.gz'):
                with gzip.open(self.vcf_file, 'rb') as f:
                    lines = f.readlines()
                # If gzip already pass as an argument pass
                # else set it to gzip
                if self.compression != 'infer':
                    pass
                else:
                    self.compression = 'gzip'
            elif self.vcf_file.endswith('.vcf'):
                with open(self.vcf_file, 'r') as f:
                    lines = f.readlines()
            else:
                sys.stderr.write('Are you sure this is a vcf file??')
                sys.exit()
        elif self.mode == 'str':
            lines = self.vcf_file.split('\n')
        # meta_info dict
        self.meta_info = OrderedDict()
        # Loop over lines
        for line in lines:
            if self.compression == 'gzip':
                line = line.decode('utf-8')
            # When it reaches the header break
            if '#CHROM' in line:
                if len(line.split('\t')) > 1:
                    sep = '\t'
                else:
                    sep = None
                self.header = [colname for colname in line.strip('#').replace(
                    '\n', '').split(sep)]
                break
            # Remove the pound and new line.  Split by equal sign.
            k, v = line.strip('#').replace('\n', '').split('=', 1)
            if k not in self.meta_info:
                self.meta_info[k] = [v]
            else:
                self.meta_info[k].append(v)

    def from_file(self):
        '''Reads the table. vdf stands for Variant Data Frame.
        comment='#' passes meta info and header.
        header is parsed from the file with header() function
        low_memory=False reads columns as object
        compression for gzip vcf files.
        chunksize pass'''
        self.variant_df = pd.read_table(self.vcf_file,
                                        comment='#',
                                        delim_whitespace=True,
                                        names=self.header,
                                        low_memory=False,
                                        compression=self.compression)
        # chunks WIP
        # if self.chunksize is not None:
        #     for vdf_chunk in pd.read_table(self.vcf_file,
        #                                    comment='#',
        #                                    names=self.header,
        #                                    low_memory=False,
        #                                    compression=compression,
        #                                    chunksize=self.chunksize):
        #         yield vdf_chunk

    def from_string(self):
        self.variant_df = pd.read_table(StringIO(self.vcf_file),
                                        comment='#',
                                        delim_whitespace=True,
                                        names=self.header,
                                        low_memory=False)

    def get_variant_info(self, concat=False, drop=False):
        # Create key set from all rows
        all_key_set = set()
        for info_row in self.vdf['INFO']:
            # This line checks if the row is empty.
            # Fix for test vcf 'string_as_flag.vcf'
            if info_row == info_row:
                for variant_info in info_row.split(';'):
                    all_key_set.add(variant_info.split('=')[0])

        # Create info_dict that will be use to create info_df
        variant_info_dict = defaultdict(list)
        for index, info_row in zip(self.vdf.index, self.vdf['INFO']):
            variant_info_dict['index'].append(index)

            # If row is empty add NA to all keys and continue with the loop.
            if not info_row == info_row:
                for key in all_key_set:
                    variant_info_dict[key].append('NA')
                continue

            # If row is not empty collect its keys to another set.
            row_key_set = set()
            for variant_info in info_row.split(';'):
                row_key_set.add(variant_info.split('=')[0])

            # Add NA to missing keys.
            for key in all_key_set - row_key_set:
                variant_info_dict[key].append('NA')

            # Add values to dict
            for variant_info in info_row.split(';'):
                try:
                    row_key, v = variant_info.split('=')
                # If Its only a flag just add 1. E.g. ;NegativeTrainError;
                except ValueError:
                    row_key = variant_info.split('=')[0]
                    v = 1
                variant_info_dict[row_key].append(v)

        variant_info_df = pd.DataFrame.from_dict(variant_info_dict)
        variant_info_df.set_index('index', inplace=True)

        if concat == True:
            new_df = pd.concat([self.vdf, variant_info_df], axis=1)
            if drop == True:
                new_df.drop('INFO', axis=1, inplace=True)
            return new_df
        else:
            return variant_info_df

    def get_sample_format(self, sample_names=None, concat=False, drop=False):
        # If sample names is None self.sample_names is used.
        if sample_names is None:
            sample_names = self.sample_names
        else:
            sample_names = sample_names

        if len(sample_names) < 1:
            sys.exit('No sample column found.')


        samples_df_dict = {}
        for sample in sample_names:

            all_key_set = set()
            for format_row in self.vdf['FORMAT']:
                if format_row == format_row:
                    for sample_format in format_row.split(':'):
                        all_key_set.add(sample_format)

            sample_format_dict = defaultdict(list)
            for index, format_row, sample_row in zip(self.vdf.index, self.vdf['FORMAT'], self.vdf[sample]):
                sample_format_dict['index'].append(index)

                row_key_set = set()
                for sample_format in format_row.split(':'):
                    row_key_set.add(sample_format)

                for key in all_key_set - row_key_set:
                    sample_format_dict[key].append('NA')

                # Whole sample row missing. Fix for test vcf "strelka.vcf"
                if not sample_row == sample_row:
                    sample_row = 'NA'

                # Sample missing format. Fix for test vcf "example-4.1.vcf"
                while len(format_row.split(':')) > len(sample_row.split(':')):
                    sample_row += ':NA'

                for sample_format, sample_info in zip(format_row.split(':'), sample_row.split(':')):
                    sample_format_dict[sample_format].append(sample_info)

            sample_format_df = pd.DataFrame.from_dict(sample_format_dict)
            sample_format_df.set_index('index', inplace=True)
            samples_df_dict[sample] = sample_format_df

        if concat == True:
            for sample in samples_df_dict.keys():
                samples_df_dict[sample] = pd.concat(
                    [self.vdf, samples_df_dict[sample]], axis=1)

            if drop == True:
                if len(samples_df_dict) == 1:
                    samples_df_dict[sample_names[0]] = samples_df_dict[sample].drop(
                        sample_names[0], axis=1)
                else:
                    for sample in sample_names:
                        samples_df_dict[sample] = samples_df_dict[sample].drop(
                            sample_names, axis=1)

            if len(samples_df_dict) == 1:
                return samples_df_dict[list(samples_df_dict.keys())[0]]
            else:
                return samples_df_dict
        else:
            if len(samples_df_dict) == 1:
                return samples_df_dict[list(samples_df_dict.keys())[0]]
            else:
                return samples_df_dict

    def to_string(self):
        self.vcf_str = StringIO()
        for info_key, info_list in self.meta_info.items():
            for info in info_list:
                self.vcf_str.write('##{}={}\n'.format(info_key, info))
        self.vdf.rename(columns={'CHROM': '#CHROM'}, inplace=True)
        self.vdf.to_csv(self.vcf_str, sep='\t', index=False)
        return self.vcf_str

    def to_vcf(self, output_file):
        with open(output_file, 'w') as f:
            f.write(self.to_string().getvalue())

    def to_table(self, output_file):
        pass
