# -*- coding: utf-8 -*-
__author__ = 'bars'

import os
import sys
import pandas as pd

hpo_link = 'http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt'


def GetHPO(hpopath):
    df = pd.read_csv(hpo_link, sep='\t', comment='#', names=[
                     'entrez-gene-id', 'gene_symbol', 'HPO_term_names', 'HPO_term_ID'])

    df = df.groupby('gene_symbol')['HPO_term_names'].apply(';'.join).reset_index()
    df.to_csv(os.path.join(hpopath, 'hpo.csv'), index=False)


def main(args):
    hpopath = os.path.join(os.path.expanduser('~'), '.cache', 'dove', 'data')
    if not os.path.exists(hpopath):
        sys.stdout.write('Downloading HPO data...\n')
        os.makedirs(hpopath, exist_ok=True)
    elif not os.path.exists(os.path.join(hpopath, 'hpo.csv')):
        sys.stdout.write('Downloading HPO data...\n')
    else:
        sys.stdout.write('Updating HPO data...\n')

    GetHPO(hpopath)

    if os.path.exists(os.path.join(hpopath, 'hpo.csv')):
        sys.stdout.write('Finished downloading hpo data.\n')
        return True
