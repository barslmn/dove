# -*- coding: utf-8 -*-
__author__ = 'bars'

import os
import sys
import json
import pandas as pd
from pandas.io.json import json_normalize
from dove.utils import updateomim
from dove.utils import updatehpo
from dove.utils import updateclinvar
from dove.utils.vcf import Vcf
from dove.utils.ensemblapi import EnsemblApi

cache_path = os.path.join(os.path.expanduser('~'), '.cache', 'dove', 'data')


def chunks(l, n):
    '''Yield successive n-size chunks from l.'''
    for i in range(0, len(l), n):
        yield l[i:i + n]


class OnlineVariantAnnotation:

    def __init__(self, vcf_file=None, mode='file', output_file=None, json_file=None, omim=False, hpo=False, clinvar=False, resume=False, query_chunks=199, read_chunks=None, table_cols=None):
        self.vcf_file = vcf_file
        self.mode = mode
        self.output_file = output_file
        self.json_file = json_file
        self.omim = omim
        self.hpo = hpo
        self.clinvar = clinvar
        self.resume = resume
        self.query_chunks = query_chunks
        self.read_chunks = read_chunks
        self.table_cols = table_cols
        self.Eapi = EnsemblApi()

        self.default_cols = [
            'CHR', 'POS', 'REF', 'ALT', 'GT', 'AD', 'DP',
            'gene_symbol', 'transcript_id', 'id',
            'minor_allele_freq', 'gnomad', 'canonical',
            'transcript_consequence_terms', 'impact', 'biotype',
            'polyphen_prediction', 'polyphen_score',
            'sift_prediction', 'sift_score',
            'hgvsc', 'hgvsp', 'clin_sig', 'pubmed'
        ]

        if self.table_cols is None:
            self.table_cols = self.default_cols
        else:
            if self.table_cols[0] == '+':
                self.table_cols = self.default_cols + self.table_cols[1:]
            else:
                pass

        if self.omim:
            if not os.path.exists(cache_path + '/omim.csv'):
                # This functions resets self.omim to True or False
                self.omim = self.try_download('omim')
                if self.omim:  # Thats why there is a second check.
                    self.table_cols.insert(
                        len(self.table_cols), 'omim_phenotypes')
            else:
                self.table_cols.insert(len(self.table_cols), 'omim_phenotypes')

        if self.hpo:
            if not os.path.exists(cache_path + '/hpo.csv'):
                self.hpo = self.try_download('HPO')
                if self.hpo:
                    self.table_cols.insert(
                        len(self.table_cols), 'HPO_term_names')
            else:
                self.table_cols.insert(len(self.table_cols), 'HPO_term_names')

        if self.clinvar:
            if not os.path.exists(cache_path + '/clinvar.vcf.gz'):
                self.clinvar = self.try_download('clinvar')
                if self.clinvar:
                    self.read_clinvar_vcf()
            else:
                self.read_clinvar_vcf()

    def read_clinvar_vcf(self):
        with Vcf(cache_path + '/clinvar.vcf.gz') as vcf:
            self.cln_df = vcf.get_variant_info(concat=True, drop=True)
        self.cln_df = self.cln_df[[col for col in list(
            self.cln_df) if col.startswith('CLN')] + ['CHROM', 'POS', 'ALT']]
        self.table_cols = self.table_cols + \
            [col for col in list(self.cln_df) if col not in vcf.generic_header]
        self.cln_df.columns = ['CHR'] + list(self.cln_df)[1:]

    def process_vcf(self):
        if not self.json_file:
            self.df_vcf['input'] = self.df_vcf['CHROM'].map(str) + ' ' + \
                self.df_vcf['POS'].map(str) + ' ' + self.df_vcf['ID'].map(str) + ' ' + self.df_vcf['REF'].map(str) + \
                ' ' + self.df_vcf['ALT'].map(str)  # create vep input
        else:
            self.df_vcf['input'] = self.df_vcf[self.header].apply(
                lambda x: '\t'.join(x.astype(str)),
                axis=1
            )

    def parse_colocated(self):
        popkeys = ['assembly_name', 'strand', 'id']
        for popkey in popkeys:  # loop over keys to pop
            for variant in self.data:  # loop over all of the data
                if popkey in variant.keys():  # check if key exists
                    variant.pop(popkey, None)  # pop unnecessary keys

        for variant in self.data:  # loop over all of the data
            if 'colocated_variants' in variant.keys():  # check if key exists
                # loop over colocated variants
                for eachdict in variant['colocated_variants']:
                    if 'frequencies' in eachdict.keys():
                        for variant_allele in eachdict['frequencies']:
                            for key, value in eachdict['frequencies'][variant_allele].items():
                                variant[key] = value
                    if 'rs' in eachdict['id']:  # if its in dbSNP
                        for key, value in eachdict.items():  # loop over dict
                            variant[key] = value  # use rs code
                    else:
                        for key, value in eachdict.items():  # loop over dict
                            # else use HGMD, EST, COSMIC vs...
                            variant[key] = value

        for variant in self.data:
            if 'colocated_variants' in variant.keys():
                variant.pop('colocated_variants', None)

    def parse_table(self):
        self.table_data = []
        consequence_list = ['transcript_consequences', 'intergenic_consequences',
                            'regulatory_feature_consequences', 'motif_feature_consequences']
        for variant in self.data:
            for c in consequence_list:
                if c in variant.keys():
                    consequences = variant.pop(c)
                    for consequence in consequences:
                        consequence['{}_terms'.format(c[:-1])] = '; '.join(
                            consequence['consequence_terms'])
                        consequence.update(variant)
                        self.table_data.append(consequence)

    def gene_base_annotation(self, path_to_table):
        '''This function is used for omim&hpo annotation.'''
        gene_table = pd.read_csv(path_to_table)
        self.df_anno_table = pd.merge(self.df_anno_table, gene_table,
                                      on='gene_symbol', how='left')

    def position_base_annotation(self, vdf):
        '''This function is used for clinvar annotation.'''
        self.df_anno_table = pd.merge(self.df_anno_table, vdf,
                                      on=['CHR', 'POS', 'ALT'], how='left')

    def write_table(self):
        if os.path.exists(self.output_file):
            self.df_anno_table.to_csv(
                self.output_file, index=False, header=False, mode='a')
        else:
            self.df_anno_table.to_csv(self.output_file, index=False,)

    def create_table(self):
        self.df_anno_table = pd.DataFrame.from_dict(
            json_normalize(self.table_data), orient='columns')

        if self.df_anno_table.empty:
            pass
        else:
            self.df_anno_table = pd.merge(
                self.df_vcf[['input', 'POS', 'GT', 'AD', 'DP']], self.df_anno_table, on='input')
            if self.json_file:
                self.df_anno_table['CHR'], self.df_anno_table['start_1'], self.df_anno_table['id_1'], self.df_anno_table['REF'], self.df_anno_table['ALT'], self.df_anno_table['rest'] = self.df_anno_table['input'].str.split(n=5).str
            else:
                self.df_anno_table['CHR'], self.df_anno_table['start_1'], self.df_anno_table['id_1'], self.df_anno_table['REF'], self.df_anno_table['ALT'] = self.df_anno_table['input'].str.split(n=4).str
            if self.omim:
                self.gene_base_annotation(os.path.join(cache_path, 'omim.csv'))
            if self.hpo:
                self.gene_base_annotation(os.path.join(cache_path, 'hpo.csv'))
            if self.clinvar:
                self.position_base_annotation(self.cln_df)

            self.add_missing_cols()
            self.df_anno_table = self.df_anno_table[[
                col for col in self.table_cols]]
        '''
        Check every cell if there is list in it:
        map everything to a string:
        then join list elements
        '''
        self.df_anno_table = self.df_anno_table.applymap(
            lambda x: ';'.join(list(map(str, x))) if isinstance(x, list) else x)
        self.df_anno_table.fillna(value={'minor_allele_freq': 0,
                                         'gnomad': 0}, inplace=True)

    def add_missing_cols(self):
        diff_cols = set(self.table_cols) - set(list(self.df_anno_table))
        if bool(diff_cols):
            for col in diff_cols:
                self.df_anno_table[col] = ""

    def try_download(self, down_term):
        sys.stdout.write('{} file does not exists!'.format(down_term))
        while True:
            userinput = input('Do you want to download it?yes/no?')
            if userinput not in ['yes', 'no', 'y', 'n']:
                sys.stdout.write(
                    'Enter yes to download or no to cancel...')
                continue
            if userinput in ['no', 'n']:
                sys.stdout.write(
                    'Setting {} to False and continuing...'.format(down_term))
                return False
            if userinput in ['yes', 'y']:
                if down_term == 'omim':
                    apikey = input('Please enter omim api key:')
                    if updateomim.main(apikey):
                        return True
                if down_term == 'HPO':
                    if updatehpo.main():
                        return True
                if down_term == 'clinvar':
                    if updateclinvar.main():
                        return True
                else:
                    sys.stdout.write('Something gone wrong. Retry?')
                    continue

    def try_resume(self):
        if os.path.exists(self.output_file):
            df_existing = pd.read_csv(self.output_file)
            df_existing['ID'] = '.'
            df_existing['input'] = df_existing['CHR'].map(str) + ' ' + \
                df_existing['POS'].map(str) + ' ' + df_existing['ID'] + ' ' + df_existing['REF'] + \
                ' ' + df_existing['ALT']  # create vep input
            return list(set(self.df_vcf['input'].tolist()) - set(df_existing['input'].tolist()))
        else:
            pass

    def annotate(self):
        with Vcf(vcf_file=self.vcf_file, mode=self.mode, chunksize=self.read_chunks) as vcf_parser:
            self.header = vcf_parser.header
            vcf_parser.vdf = vcf_parser.get_sample_format(concat=True)
            if self.read_chunks is None:
                self.df_vcf = vcf_parser.vdf
                if self.output_file is not None:
                    self.annotate_to_file()
                else:
                    return self.annotate_to_df()
            else:
                for chunk in vcf_parser.vdf:
                    self.df_vcf = chunk
                    if self.output_file is not None:
                        self.annotate_to_file()
                    else:
                        return self.annotate_to_df()

    def annotate_to_file(self):
        self.process_vcf()

        sys.stdout.write('Started ensembl VEP for {}\nNumber of variants {}\n'.format(
            self.header[-1], len(self.df_vcf)))

        if not self.resume and os.path.exists(self.output_file) and self.read_chunks is None:
            os.remove(self.output_file)

        if self.resume and os.path.exists(self.output_file):
            sys.stdout.write('Resuming...\n')
            query_variants = self.try_resume()
        else:
            query_variants = self.df_vcf['input'].tolist()

        if self.json_file is None:
            for chunk in chunks(query_variants, self.query_chunks):
                self.data = self.Eapi.ensemblapi(chunk)
                self.parse_colocated()
                self.parse_table()
                self.create_table()
                self.write_table()
        else:
            with open(self.json_file) as f:
                self.jdata = json.loads('[{}]'.format(
                    ','.join(f.read().splitlines())))
            for chunk in chunks(self.jdata, self.query_chunks):
                self.data = chunk
                self.parse_colocated()
                self.parse_table()
                self.create_table()
                self.write_table()
        sys.stdout.write(
            'Finished annotation for {}\n'.format(self.header[-1]))

    def annotate_to_df(self):
        self.process_vcf()
        query_variants = self.df_vcf['input'].tolist()
        df_anno_tables = []

        if self.json_file is None:
            for chunk in chunks(query_variants, self.query_chunks):
                self.data = self.Eapi.ensemblapi(chunk)
                self.parse_colocated()
                self.parse_table()
                self.create_table()
                df_anno_tables.append(self.df_anno_table)
        else:
            with open(self.json_file) as f:
                self.jdata = json.loads('[{}]'.format(
                    ','.join(f.read().splitlines())))
            for chunk in chunks(self.jdata, self.query_chunks):
                self.data = chunk
                self.parse_colocated()
                self.parse_table()
                self.create_table()
                df_anno_tables.append(self.df_anno_table)
        return pd.concat(df_anno_tables)


def main(args):
    OVA = OnlineVariantAnnotation(
        vcf_file=args.input,
        output_file=args.output,
        json_file=args.json,
        omim=args.omim,
        hpo=args.hpo,
        clinvar=args.clinvar,
        resume=args.resume,
        query_chunks=args.query_chunks,
        read_chunks=args.read_chunks,
        table_cols=args.table_cols
    )
    OVA.annotate()
