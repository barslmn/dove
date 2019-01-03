# -*- coding: utf-8 -*-
__author__ = 'bars'

import os
import sys
import urllib.request

clinvar_link = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
clinvar_path = os.path.join(os.path.expanduser('~'), '.cache', 'dove', 'data')


def get_clinvar():
    urllib.request.urlretrieve(clinvar_link, os.path.join(clinvar_path, 'clinvar.vcf.gz'))


def main(args):
    if not os.path.exists(clinvar_path):
        sys.stdout.write('Downloading Clinvar data...\n')
        os.makedirs(clinvar_path, exist_ok=True)
    elif not os.path.exists(os.path.join(clinvar_path, 'clinvar.vcf.gz')):
        sys.stdout.write('Downloading Clinvar data...\n')
    else:
        sys.stdout.write('Updating Clinvar data...\n')

    get_clinvar()

    if os.path.exists(os.path.join(clinvar_path, 'clinvar.vcf.gz')):
        sys.stdout.write('Finished downloading clinvar data.\n')
        return True
