#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from dove.utils import updateomim
from dove.utils import updatehpo
from dove.utils import updateclinvar
from dove.utils import copyexamplefiles
from dove.core import variantfilter
from dove.core import multigenomeanalysis
from dove.core import onlinevariantannotation
from dove import __version__


def get_args():
    parser = argparse.ArgumentParser()

    # Not functional
    parser.add_argument('-v', '--verbose', required=False, action='store_true',
                        help='WIP.Write what is happening to stdout.')

    parser.add_argument('--version',
                        action='version', version=__version__)

    subparsers = parser.add_subparsers(dest='tool')

    parser_OVA = subparsers.add_parser(
        'OVA',
        formatter_class=RawTextHelpFormatter,
        help='Online Variant Annotation tool. For help type: dove OVA --help',
        epilog='''Usage examples:
        dove OVA -i sample.vcf -o sample.csv -q 199
        dove OVA -i sample.vcf -o sample.csv -q 199 --omim --resume -c + exon genesplicer'''
    )
    parser_OVA.add_argument('-i', '--input', required=True,
                            default=None, type=str, help='input vcf file.')
    parser_OVA.add_argument('-o', '--output', required=True,
                            default=None, type=str, help='output annotation file.')
    parser_OVA.add_argument('-q', '--query-chunks', required=True, type=int,
                            help='''Number of variants to be requested with each post.
                            Maximum POST size here:https://rest.ensembl.org/documentation/info/vep_region_post''')
    parser_OVA.add_argument('-k', '--read-chunks', required=False, type=int,
                            help='''Read vcf in chunks. For large files.''')
    parser_OVA.add_argument('-c', '--table-cols', required=False, type=str, nargs='+',
                            help='''Columns to be included in the annotation table.
                            For additional columns to default columns use + before entering columns.
                            Possible columns:
                            CHR POS LOC REF ALT GT AD DP
                            amino_acids biotype blosum62 canonical
                            ccds cdna_end cdna_start cds_end cds_start
                            clin_sig codons transcript_consequence_terms csn distance domains
                            exon flags gene_id gene_symbol gene_symbol_source genesplicer
                            minor_allele minor_allele_freq
                            aa_allele aa amr_allele amr
                            afr_allele afr allele_string
                            ea_allele ea eas_allele eas
                            eur_allele eur sas_allele sas 
                            gnomad_allele gnomad 
                            gnomad_afr_allele gnomad_afr gnomad_amr_allele gnomad_amr
                            gnomad_asj_allele gnomad_asj gnomad_eas_allele gnomad_eas
                            gnomad_fin_allele gnomad_fin gnomad_nfe_allele gnomad_nfe
                            gnomad_oth_allele gnomad_oth gnomad_sas_allele gnomad_sas
                            hgnc_id hgvsc hgvsp id impact intron
                            motif_feature_consequence_terms motif_feature_id
                            motif_name motif_pos motif_score_change
                            maxentscan_alt maxentscan_diff maxentscan_ref
                            most_severe_consequence phenotype_or_disease
                            polyphen_prediction polyphen_score
                            protein_end protein_id protein_start
                            pubmed regulatory_feature_consequence_terms
                            seq_region_name sift_prediction sift_score
                            somatic start end strand swissprot
                            transcript_id trembl uniparc variant_allele variant_class''')
    parser_OVA.add_argument('-m', '--omim', required=False, dest='omim', action='store_true',
                            help='Enables gene annotation to omim diseases.')
    parser_OVA.add_argument('-p', '--hpo', required=False, dest='hpo', action='store_true',
                            help='Enables gene annotation to HPO.')
    parser_OVA.add_argument('-l', '--clinvar', required=False, dest='clinvar', action='store_true',
                            help='Enables variant annotation to Clinvar.')
    parser_OVA.add_argument('-r', '--resume', required=False, dest='resume', action='store_true',
                            help='Resumes interrupted annotation process.')

    parser_MGA = subparsers.add_parser(
        'MGA', help='Multi Genome Analysis tool. For help type: dove MGA --help')
    parser_MGA.add_argument('-i', '--index', type=str,
                            required=False, default=None,
                            help='enter index annotation file location.')
    parser_MGA.add_argument('-d', '--pedigree', type=str,
                            required=False, default=None,
                            help='enter pedigree file.')
    parser_MGA.add_argument('-o', '--output', type=str,
                            required=True, default=None,
                            help='enter index annotation file location.')
    parser_MGA.add_argument('-m', '--mother', type=str,
                            required=False, default=None,
                            help='enter parent annotation file location.')
    parser_MGA.add_argument('-f', '--father', type=str,
                            required=False, default=None,
                            help='enter parent annotation file location.')
    parser_MGA.add_argument('-a', '--analysistype', type=str,
                            required=True, default=None,
                            help='analysis type, siblings, trio, or both.')
    parser_MGA.add_argument('-p', '--pattern', type=str,
                            required=False, default=None,
                            help='enter possible inheritance pattern.Options are:AR, AD, DN, CH, XLD, XLR, YL')

    parser_VF = subparsers.add_parser(
        'VF',
        formatter_class=RawTextHelpFormatter,
        help='''Variant filtering tool.
        Filtering columns is done with -c flag.
        -c can take multiple arguments.
        The first argument it takes must be the column name of interest.
        The second argument must be filtering option.
        Any argument after first two will be used as filtering criteria.
        There are four types of filtering options:

        'Equal' and 'Not Equal'
        shorthen as 'eq' and 'ne'.
        These can be used to filter both string and numeric fields.

        'Less Than', 'Less Than or Equal', 'Greater Than', and 'Greater Than or Equal'
        shorthen 'as', 'lt', 'le', 'gt', and 'ge'
        These can be used to filter numeric fields.

        'Includes' and 'Excludes'
        shorthen as 'in' and 'ex'
        These can be used to search a column that includes or excludes a specific substring.

        Lastly there is the bed option.
        Right now it returns variants in between bed file intervals.
        Note: Annotation table must have LOC column formatted as start:end.
        ''',

        epilog='''Usage examples:
        dove VF -i sample.csv -o output_sample.csv -c gnomad lt 0.01 -c GT ne 0/1 -c biotype eq protein_coding -c HPO_term_names ex tumor
        dove VF -i sample.csv -o output_sample.csv -d CHR LOC gene_symbol transcript_consequence_terms
        dove VF -i sample.vcf -o output_sample.vcf -c FILTER in PASS -c QUAL ge 100
        dove VF -i sample.vcf -o output_sample.vcf -b bed_file.bed'''
    )
    parser_VF.add_argument('-i', '--input', required=True,
                           default=None, type=str,
                           help='input ensembl annotation file.')
    parser_VF.add_argument('-o', '--output', required=True,
                           default=None, type=str,
                           help='output annotation file.')
    parser_VF.add_argument('-b', '--bed-file', required=False,
                           default=None, type=str,
                           help='Bed file.')
    parser_VF.add_argument('-c', '--column', nargs='+', action='append')
    parser_VF.add_argument('-d', '--drop', nargs='+')

    parser_examplefile = subparsers.add_parser(
        'examplefile', help='copies given example utility file to current directory. For help type: examplefile --help')
    parser_examplefile.add_argument('-e', '--examplefile', required=True,
                                    default=None, type=str,
                                    help='pedigree, filterfile or both')

    parser_updateomim = subparsers.add_parser(
        'updateomim', help='install/update omim table')
    parser_updateomim.add_argument('-k', '--apikey', required=True,
                                   default=None, type=str,
                                   help='Please enter your api key. For how to get one: https://omim.org/api')

    parser_updatehpo = subparsers.add_parser(
        'updatehpo', help='install/update hpo table')

    parser_updateclinvar = subparsers.add_parser(
        'updateclinvar', help='install/update clinvar vcf')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    if args.tool == 'OVA':
        onlinevariantannotation.main(args)
    if args.tool == 'MGA':
        multigenomeanalysis.main(args)
    if args.tool == 'VF':
        variantfilter.main(args)
    if args.tool == 'examplefile':
        copyexamplefiles.main(args)

    if args.tool == 'updateomim':
        if updateomim.main(args):
            sys.stdout.write('Successfully downloaded omim data.\n')
        else:
            sys.stdout.write('Something gone wrong. Please try again.\n')
    if args.tool == 'updatehpo':
        if updatehpo.main(args):
            sys.stdout.write('Successfully downloaded hpo data.\n')
        else:
            sys.stdout.write('Something gone wrong. Please try again.\n')
    if args.tool == 'updateclinvar':
        if updateclinvar.main(args):
            sys.stdout.write('Successfully downloaded clinvar vcf.\n')
        else:
            sys.stdout.write('Something gone wrong. Please try again.\n')


if __name__ == '__main__':
    #min_version = (3, 6, 0)
    #if sys.version_info < min_version:
    #    sys.exit('Python {}.{}.{} or later is required.\n'.format(min_version))
    main()
