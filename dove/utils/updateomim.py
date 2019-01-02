# -*- coding: utf-8 -*-
__author__ = 'bars'

import os
import sys
import json
import urllib
import pandas as pd
from pandas.io.json import json_normalize


def chunks(l, n):
    '''Yield successive n-size chunks from l.'''
    for i in range(0, len(l), n):
        yield l[i:i + n]


class OmimApi():

    """Docstring for OmimApi. """

    def __init__(self, apikey):
        """TODO: to be defined1. """
        self.apikey = apikey

    def mimnumbers(self):
        try:
            dfmim2gene = pd.read_table(
                'https://omim.org/static/omim/data/mim2gene.txt', skiprows=4)
        except urllib.error.URLError:
            raise SystemExit('Check internet connection.')
        return dfmim2gene['# MIM Number'].tolist()

    def omimapi(self):
        omimoutput = []
        for chunk in chunks(self.mimnumbers(), 100):
            request_data = {}
            request_data['apiKey'] = self.apikey
            request_data['format'] = 'json'
            request_data['mimNumber'] = ','.join(str(i) for i in chunk)

            url = 'http://api.omim.org/api/geneMap'
            # add parameters to url string
            url_values = urllib.parse.urlencode(request_data)
            url = url + '?' + url_values
            # query OMIM
            try:
                response = urllib.request.urlopen(url)
            except urllib.error.HTTPError:
                raise SystemExit("Failed to access OMIM API, valid key?")
            # read in response
            result = json.loads(response.read().decode('utf-8'))
            omimoutput = omimoutput + \
                result['omim']['listResponse']['geneMapList']
        return omimoutput


def parseomim(omimdata):
    omimpheno = [genemap['geneMap']
                 for genemap in omimdata if 'phenotypeMapList' in genemap['geneMap'].keys()]

    for geneMap in omimpheno:
        for key in geneMap['phenotypeMapList'][0]['phenotypeMap']:
            if key in ['phenotype']:
                geneMap[key] = ';'.join([phenotypeMap['phenotypeMap'][key]
                                         for phenotypeMap in geneMap['phenotypeMapList']])

    genesympheno = [{k: v for k, v in genemap.items(
    ) if k in ['geneSymbols', 'phenotype']} for genemap in omimpheno]

    genepheno = []
    for genemap in genesympheno:
        genemap['gene_symbol'] = genemap.pop('geneSymbols')
        genemap['phenotypes'] = genemap.pop('phenotype')
        genepheno.append(genemap)
    return genepheno


def saveomim(genepheno, omimpath):
    df = pd.DataFrame.from_dict(json_normalize(genepheno), orient='columns')
    df = pd.concat([pd.Series(row['phenotypes'], row['gene_symbol'].split(', '))
                    for _, row in df.iterrows()]).reset_index()
    df.columns = ['gene_symbol', 'omim_phenotypes']
    df['omim_phenotypes'] = df['omim_phenotypes'].str.replace('{|}', '')
    df['omim_phenotypes'] = df['omim_phenotypes'].str.replace('[|]', '')
    df.to_csv(os.path.join(omimpath, 'omim.csv'), index=False)


def main(args):
    if type(args) == str:
        apikey = args
    else:
        apikey = args.apikey

    omimpath = os.path.join(os.path.expanduser('~'), '.cache', 'dove', 'data')
    if not os.path.exists(omimpath):
        sys.stdout.write('Downloading omim data...')
        os.makedirs(omimpath, exist_ok=True)
    elif not os.path.exists(os.path.join(omimpath, 'omim.csv')):
        sys.stdout.write('Downloading omim data...')
    else:
        sys.stdout.write('Updating omim data...')

    omimdata = OmimApi(apikey).omimapi()
    genepheno = parseomim(omimdata)
    saveomim(genepheno, omimpath)

    if os.path.exists(os.path.join(omimpath, 'omim.csv')):
        sys.stdout.write('Finished downloading omim data.\n')
        return True
