# -*- coding: utf-8 -*-
__author__ = 'bars'

from dove.utils.vcf import Vcf
import pandas as pd


def fzygo(df, zygosity):
    if zygosity == '0/1':
        return df[df['GT'] != '1/1']
    if zygosity == '1/1':
        return df[df['GT'] != '0/1']


def autosomal_recessive(index, mother, father):
    # Disease causing variant must be homozygote in index and heterozygote in parents
    index_hom = fzygo(index, '1/1')
    mother_het = fzygo(mother, '0/1')
    father_het = fzygo(father, '0/1')

    # Get homozygote variants only if they're in both parents
    index_hom_mother_het = pd.merge(index_hom, mother_het[[
                                    'Position', 'transcript_id']], on=['Position', 'transcript_id'])
    index_hom_mother_het_father_het = pd.merge(index_hom_mother_het, father_het[[
                                               'LOC', 'transcript_id']], on=['LOC', 'transcript_id'])

    # Remove position column and return
    index_hom_mother_het_father_het.drop('Position', axis=1, inplace=True)
    return index_hom_mother_het_father_het


def autosomal_dominant(index, mother, father):
    # Get heterozygote variants from all
    index_het = fzygo(index, '0/1')
    mother_het = fzygo(mother, '0/1')
    father_het = fzygo(father, '0/1')

    # Get heterozygote varaints from both parents and concat
    index_het_mother_het = pd.merge(index_het, mother_het[[
                                    'Position', 'transcript_id']], on=['Position', 'transcript_id'])
    index_het_mother_het.insert(7, 'Variant Source', 'm')
    index_het_father_het = pd.merge(index_het, father_het[[
                                    'Position', 'transcript_id']], on=['Position', 'transcript_id'])
    index_het_father_het.insert(7, 'Variant Source', 'f')
    index_het_mother_het_index_het_father_het = pd.concat(
        [index_het_mother_het, index_het_father_het],
        keys=['x', 'y'])

    # Remove position column and sort
    index_het_mother_het_index_het_father_het.sort_values(by=['CHR', 'LOC'], ascending=[
        True, True], inplace=True)
    index_het_mother_het_index_het_father_het.drop(
        'Position', axis=1, inplace=True)
    return index_het_mother_het_index_het_father_het


def compound_heterozygote(index, mother, father):
    # Get heterozygote variants of all samples
    index_het = fzygo(index, '0/1')
    mother_het = fzygo(mother, '0/1')
    father_het = fzygo(father, '0/1')

    # Get shared variants of index and mother and mark as m
    index_het_mother_het = pd.merge(
        index_het, mother_het[['Position', 'transcript_id']], on=['Position', 'transcript_id'])
    index_het_mother_het.insert(7, 'Variant Source', 'm')
    # Get shared variants of index and father and mark as f
    index_het_father_het = pd.merge(
        index_het, father_het[['Position', 'transcript_id']], on=['Position', 'transcript_id'])
    index_het_father_het.insert(7, 'Variant Source', 'f')
    # Get denovo variants
    index_denovo = de_novo(index, mother, father, drop_pos=False)

    # Get shared genes and drop duplicates. This creates a data frame with just gene names.
    MF_compound_genes = pd.merge(index_het_mother_het[[
        'gene_symbol', 'transcript_id']], index_het_father_het[['gene_symbol', 'transcript_id']], on=['gene_symbol', 'transcript_id'])
    MDN_compound_genes = pd.merge(index_het_mother_het[[
        'gene_symbol', 'transcript_id']], index_denovo[['gene_symbol', 'transcript_id']], on=['gene_symbol', 'transcript_id'])
    FDNcompound_genes = pd.merge(index_het_father_het[[
        'gene_symbol', 'transcript_id']], index_denovo[['gene_symbol', 'transcript_id']], on=['gene_symbol', 'transcript_id'])
    compound_genes = pd.concat(
        [MF_compound_genes, MDN_compound_genes, FDNcompound_genes], ignore_index=True)
    compound_genes.drop_duplicates(
        ['gene_symbol', 'transcript_id'], inplace=True)

    # Removo not compound genes
    index_het_mother_het = pd.merge(index_het_mother_het, compound_genes, on=[
                                    'gene_symbol', 'transcript_id'])
    index_het_father_het = pd.merge(index_het_father_het, compound_genes, on=[
                                    'gene_symbol', 'transcript_id'])
    index_denovo_variants = pd.merge(index_denovo, compound_genes, on=[
                                     'gene_symbol', 'transcript_id'])

    # concat two dataframes and sort by chromosome and gene name
    index_compound_het = pd.concat(
        [index_het_mother_het, index_het_father_het, index_denovo_variants], ignore_index=True)
    index_compound_het.sort_values(by=['CHR', 'LOC', 'gene_symbol'], ascending=[
                                   True, True, True], inplace=True)
    index_compound_het.drop('Position', axis=1, inplace=True)
    return index_compound_het


def de_novo(index, mother, father, drop_pos=True):
    # Get variants not contain in mother
    index_mother = index[~index.Position.isin(mother.Position)]
    # Get variants not contain in father
    index_mother_father = index_mother[~index_mother.Position.isin(
        father.Position)]

    # Remove position column
    index_mother_father.insert(7, 'Variant Source', 'dn')
    if drop_pos is not False:
        index_mother_father = index_mother_father.drop('Position', axis=1)
    return index_mother_father


def XL_recessive(index, mother, father):
    pass


def XL_dominant(index, mother, father):
    pass


def Y_linked(index, mother, father):
    pass


def preprocess_table(member):
    if member['file_path'].endswith('csv'):
        df = pd.read_csv(member['file_path'], low_memory=False)
        df['Position'] = df['CHR'].map(str) + ':' + df['LOC'].map(str)
    if member['file_path'].endswith(('vcf', 'vcf.gz',)):
        with Vcf(member['file_path']) as vcf:
            vcf.vdf = vcf.get_sample_format(concat=True)
            df = vcf.vdf
            df['Position'] = df['CHROM'].map(str) + ':' + df['POS'].map(str)
    return df


def trio(index, mother, father, pattern):
    df_index = preprocess_table(index)
    df_mother = preprocess_table(mother)
    df_father = preprocess_table(father)

    if pattern == 'AR':
        df = autosomal_recessive(df_index, df_mother, df_father)
    if pattern == 'AD':
        df = autosomal_dominant(df_index, df_mother, df_father)
    if pattern == 'CH':
        df = compound_heterozygote(df_index, df_mother, df_father)
    if pattern == 'DN':
        df = de_novo(df_index, df_mother, df_father)
    return df


def siblings(index, siblings_list, df_index=None):
    if df_index is None:
        df_index = preprocess_table(index)
    df_index_het = fzygo(df_index, '0/1')
    df_index_hom = fzygo(df_index, '1/1')

    for sibling in siblings_list:
        if sibling['condition'] == 'affected':
            df_asibling = preprocess_table(sibling)
            df_asibling_het = fzygo(df_asibling, '0/1')
            df_asibling_hom = fzygo(df_asibling, '1/1')

            df_index_hom = pd.merge(
                df_index_hom,
                df_asibling_hom[['Position', 'transcript_id']],
                on=['Position', 'transcript_id']
            )
            df_index_het = pd.merge(
                df_index_het,
                df_asibling_het[['Position', 'transcript_id']],
                on=['Position', 'transcript_id']
            )

            df_index = pd.concat(
                [df_index_hom, df_index_het], sort=False, ignore_index=True)

        if sibling['condition'] == 'normal':
            df_nsibling = preprocess_table(sibling)
            df_nsibling_het = fzygo(df_nsibling, '0/1')
            df_index_het = pd.merge(
                df_index_het,
                df_nsibling_het[['Position', 'transcript_id']],
                on=['Position', 'transcript_id']
            )

            df_index = df_index[~(df_index['Position'] + df_index['transcript_id']).isin(
                df_nsibling['Position'] + df_nsibling['transcript_id'])]
            df_index = df_index.append(df_index_het, ignore_index=True)

    df_index.drop('Position', axis=1, inplace=True)

    df_index.sort_values(
        by=['CHR', 'LOC'],
        ascending=[True, True],
        inplace=True
    )
    return df_index


def check_input(relatives):
    # Prettify
    relatives = [
        {
            'generation': relative[0],
            'degree': relative[1],
            'sex': relative[2],
            'condition': relative[3],
            'file_path': relative[4]
        } for relative in relatives
    ]
    return relatives


def get_parents(relatives):
    for relative in relatives:
        if relative['generation'] == '-1' and relative['degree'] == '1' and relative['sex'] == 'male':
            father = relative
        if relative['generation'] == '-1' and relative['degree'] == '1' and relative['sex'] == 'female':
            mother = relative
    return father, mother


def get_siblings(relatives):
    return [relative for relative in relatives if relative['generation'] == '0' and relative['degree'] == '2']


def main(args):
    index = {'sex': args.index[0], 'file_path': args.index[1]}
    output = args.output
    relatives = check_input(args.relative)
    analysis = args.analysis
    pattern = args.pattern

    if analysis == 'trio':
        father, mother = get_parents(relatives)
        trio(index, father, mother, pattern).to_csv(output, index=False)

    if analysis == 'siblings':
        siblings_list = get_siblings(relatives)
        siblings(index, siblings_list).to_csv(output, index=False)

    if analysis == 'both':
        pass
