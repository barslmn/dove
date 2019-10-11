#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from dove.utils import updateomim
from dove.utils import updatehpo
from dove.utils import updateclinvar
from dove.core import variantfilter
from dove.core import infamilyanalysis
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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='Online Variant Annotation tool. For help type: dove OVA -h',
        description='''
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
        transcript_id trembl uniparc variant_allele variant_class''',
        epilog='''Usage examples:
        dove OVA -i sample.vcf -o sample.csv -q 199
        dove OVA -i sample.vcf -o sample.csv -q 199 --omim --resume -c + exon genesplicer'''
    )
    parser_OVA.add_argument('-i', '--input', required=True,
                            default=None, type=str, help='input vcf file.')
    parser_OVA.add_argument('-o', '--output', required=True,
                            default=None, type=str, help='output annotation file.')
    parser_OVA.add_argument('-j', '--json', required=False,
                            default=None, type=str, help='Localy pre VEP annotated json file')
    parser_OVA.add_argument('-q', '--query-chunks', required=True, type=int,
                            help='''Number of variants to be requested with each post.
                            Maximum POST size here:https://rest.ensembl.org/documentation/info/vep_region_post''')
    parser_OVA.add_argument('-k', '--read-chunks', required=False, type=int,
                            help='''Read vcf in chunks. For large files.''')
    parser_OVA.add_argument('-c', '--table-cols', required=False, type=str, nargs='+',
                            help='''Columns to add. For option type: dove OVA -h''')
    parser_OVA.add_argument('-m', '--omim', required=False, dest='omim', action='store_true',
                            help='Enables gene annotation to omim diseases.')
    parser_OVA.add_argument('-p', '--hpo', required=False, dest='hpo', action='store_true',
                            help='Enables gene annotation to HPO.')
    parser_OVA.add_argument('-l', '--clinvar', required=False, dest='clinvar', action='store_true',
                            help='Enables variant annotation to Clinvar.')
    parser_OVA.add_argument('-r', '--resume', required=False, dest='resume', action='store_true',
                            help='Resumes interrupted annotation process.')

    parser_IFA = subparsers.add_parser(
        'IFA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='''In Family Analysis Tool. For more type: dove IFA -h''',
        description=''' In Family Analysis tool.
        -r  Relatives are specified using -r option. This option requires 4 kinds of information.

            1. Relatedness info
            --------------------
            Relatedness is given using two numbers.
                1.1 First number defines the generation.
                Index's generations is taken as 0. generation.
                Generations below index is incremented by one and
                generations above are decremented by one.
                For example, index's siblings and cousins(no removed) are 0th generation,
                index parents and uncles/aunts are -1. generation, and
                index's children and nephews/nieces are +1. generations.

                1.2 Second number defines the degree of relationship between
                index and given relative. See:
                https://en.wikipedia.org/wiki/Consanguinity#/media/File:Table_of_Consanguinity_showing_degrees_of_relationship.svg

            2. Sex
            --------------------
            Options are:
                2.1 male
                2.2 female
                2.3 undefined.

            3. Condition
            --------------------
            Options are:
                3.1 normal
                3.2 affected

            4. Annotation File
            --------------------
            Full path to the given relatives annotation file.

            Input of normal parent would look like,
            for father,
                -r -1 1 male normal father.csv
            for mother,
                -r -1 1 female normal mother.csv
            Input of affected sibling would look like,
                -r 0 1 female affected sister.csv
            Input of affected great grand uncle would look like,
                -r -3 5 male affected greatgranduncle.csv

        -i  Since generations and degrees are based on index and index's condition assumed to be affected index option only requires sex and annotation file.
            Example:
                -i male index.csv
        ''',

        epilog='''Usage examples:
        For trio with recessive trait,
        dove IFA -i male index.csv -o index_AR.csv -r -1 1 male normal father.csv -r -1 1 female normal mother.csv -a trio -p AR
        For affected sibling,
        dove IFA -i female index.csv -o index_shared.csv -r 0 1 male affected brother.csv -a siblings
        '''
    )
    parser_IFA.add_argument('-i', '--index', type=str, nargs='+',
                            required=False, default=None,
                            help='enter index annotation file location.')
    parser_IFA.add_argument('-o', '--output', type=str,
                            required=True, default=None,
                            help='enter index annotation file location.')
    parser_IFA.add_argument('-r', '--relative', nargs='+', action='append')
    parser_IFA.add_argument('-a', '--analysis', type=str,
                            required=True, default=None,
                            help='analysis type, siblings, trio, or both.')
    parser_IFA.add_argument('-p', '--pattern', type=str,
                            required=False, default=None,
                            help='enter possible inheritance pattern.Options are:AR, AD, DN, CH, XLD, XLR, YL')

    parser_VF = subparsers.add_parser(
        'VF',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help='''Variant Filtering tool. For more info type: dove VF -h''',
        description='''
        Variant filtering tool.

        Filtering columns is done with -c flag.
        -c can take multiple arguments.
        The first argument it takes must be the column name of interest.
        The second argument must be filtering option.
        Any argument after first two will be used as filtering criteria.
        There are four types of filtering options:

        'Equal' and 'Not Equal'
        -----------------------
        Shorthen as 'eq' and 'ne'.
        These can be used to filter both string and numeric fields.
        To filter empty fields 'null' can be used. see examples

        'Less Than', 'Less Than or Equal', 'Greater Than', and 'Greater Than or Equal'
        -----------------------
        Shorthen as, 'lt', 'le', 'gt', and 'ge'.
        These can be used to filter numeric fields.

        'Includes' and 'Excludes'
        -----------------------
        Shorthen as 'in' and 'ex'.
        These can be used to search a column that includes or excludes a specific substring.
        Any string containing sigle quote should be surrounded by double quotes. vice versa. see examples.

        Lastly there is the bed option.
        Right now it returns variants in between bed file intervals.
        Note: Annotation table must have LOC column formatted as start:end.
        ''',

        epilog='''Usage examples:
        dove VF -i sample.csv -o output_sample.csv -c gnomad lt 0.01 -c GT ne 0/1 -c biotype eq protein_coding -c HPO_term_names ex tumor
        dove VF -i sample.csv -o output_sample.csv -d CHR LOC gene_symbol transcript_consequence_terms
        dove VF -i sample.vcf -o output_sample.vcf -c FILTER in PASS -c QUAL ge 100
        dove VF -i sample.vcf -o novel_variants.vcf -c id eq null
        dove VF -i sample.vcf -o longsearchterm.vcf -c omim_phenotypes in "Bar's disease"
        dove VF -i sample.vcf -o output_sample.vcf -b bed_file.bed
        '''
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
    parser_VF.add_argument('-d', '--drop', nargs='+',
                           help='drop rows based on selected columns')
    parser_VF.add_argument('-k', '--keep', nargs='+',
                           help='keep columns drop other columns')

    parser_updateomim = subparsers.add_parser(
        'updateomim', help='install/update omim table')
    parser_updateomim.add_argument('-k', '--apikey', required=True,
                                   default=None, type=str,
                                   help='Please enter your api key. For how to get one: https://omim.org/api')

    subparsers.add_parser(
        'updatehpo', help='install/update hpo table')

    subparsers.add_parser(
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
    if args.tool == 'IFA':
        infamilyanalysis.main(args)
    if args.tool == 'VF':
        variantfilter.main(args)

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
    main()
