# -*- coding: utf-8 -*-
__author__ = 'bars'

import json
import requests
import sys


class EnsemblApi():

    """Docstring for EnsemblApi. """

    def __init__(self):
        """TODO: to be defined1. """
        pass

    def ensemblapi(self, variantlist):
        server = "http://grch37.rest.ensembl.org"
        ext = "/vep/homo_sapiens/region"
        headers = {"Content-Type": "application/json",
                   "Accept": "application/json"}
        # print('Sending variant request of size {} to ensembl API'.format(
        r = requests.post(server + ext, headers=headers,
                data='{ "variants" : ' + json.dumps(variantlist) + ', "Blosum62": 1, "CSN": 1, "GeneSplicer": 1, "MaxEntScan": 1, "appris": 1, "canonical": 1, "ccds": 1, "dbNSFP": "LRT_pred,MutationTaster_pred", "dbscSNV": 1, "domains": 1, "failed": 1, "hgvs": 1, "miRNA": 1, "numbers": 1, "protein": 1, "refseq":1, "tsl": 1, "uniprot": 1, "variant_class": 1}')

        if not r.ok:
            r.raise_for_status()
            sys.exit('Something wrong. Check connection maybe??')
        return r.json()

