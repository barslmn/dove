# -*- coding: utf-8 -*-
__author__ = 'bars'

'''
Tool for filtering ensembl annotation tables.
'''

__author__ = 'bars'

import pandas as pd
from dove.utils.parsepedigree import ParsePedigree

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


def preprocess_table(path):
    df = pd.read_csv(path, low_memory=False)
    df['Position'] = df['CHR'].map(str) + ':' + df['LOC'].map(str)
    return df


def trio(index_path, mother_path, father_path, pattern):
    dfindex = preprocess_table(index_path)
    dfmother = preprocess_table(mother_path)
    dffather = preprocess_table(father_path)

    if pattern == 'AR':
        df = autosomal_recessive(dfindex, dfmother, dffather)
    if pattern == 'AD':
        df = autosomal_dominant(dfindex, dfmother, dffather)
    if pattern == 'CH':
        df = compound_heterozygote(dfindex, dfmother, dffather)
    if pattern == 'DN':
        df = de_novo(dfindex, dfmother, dffather)
    return df


def siblings(index, family, dfindex=None):
    if dfindex is None:
        dfindex = preprocess_table(index['path'])
    dfindex_het = fzygo(dfindex, '0/1')
    dfindex_hom = fzygo(dfindex, '1/1')
    for i, member in enumerate(family):
        if member['index'] == 1:
            family.pop(i)

    for member in family:
        if member['mother'] == index['mother'] and member['father'] == index['father']:
            if member['affected'] == 1:
                dfasibling = preprocess_table(member['path'])
                dfasibling_het = fzygo(dfasibling, '0/1')
                dfasibling_hom = fzygo(dfasibling, '1/1')
                dfindex_hom = pd.merge(dfindex_hom, dfasibling_hom[['Position', 'transcript_id']], on=[
                    'Position', 'transcript_id'])
                dfindex_het = pd.merge(dfindex_het, dfasibling_het[['Position', 'transcript_id']], on=[
                    'Position', 'transcript_id'])
                dfindex = pd.concat(
                    [dfindex_hom, dfindex_het], ignore_index=True)

            if member['affected'] == 0:
                dfusibling = preprocess_table(member['path'])
                dfusibling_het = fzygo(dfusibling, '0/1')
                dfindex_het = pd.merge(dfindex_het, dfusibling_het[['Position', 'transcript_id']], on=[
                    'Position', 'transcript_id'])
                # Following line works but is stupid and should be changed.
                dfindex = dfindex[~(dfindex['Position'] + dfindex['transcript_id']).isin(
                    dfusibling['Position'] + dfusibling['transcript_id'])]
                dfindex = dfindex.append(dfindex_het, ignore_index=True)
    dfindex.drop('Position', axis=1, inplace=True)
    dfindex.sort_values(by=['CHR', 'LOC'], ascending=[
        True, True], inplace=True)
    return dfindex


def main(args):
    index_path = args.index
    pedigree = args.pedigree
    outputtable = args.output
    mother_path = args.mother
    father_path = args.father
    analysis = args.analysistype
    pattern = args.pattern

    if pedigree is not None:
        ped = ParsePedigree(pedigree)
        family = ped.parsepedigree()

        for member in family:
            if member['index'] == 1:
                index = member
                index_path = index['path']
        if analysis is not 'siblings':
            for member in family:
                if member['id'] == index['mother']:
                    mother_path = member['path']
            for member in family:
                if member['id'] == index['father']:
                    father_path = member['path']

    if analysis == 'trio':
        trio(index_path, mother_path, father_path, pattern).to_csv(
            outputtable, index=False)

    if analysis == 'siblings':
        siblings(index, family).to_csv(outputtable, index=False)

    if analysis == 'both':
        dfindex = trio(index_path, mother_path, father_path, pattern).to_csv(
            outputtable, index=False)
        siblings(index, family, dfindex).to_csv(outputtable, index=False)

class MultiGenomeAnalysis(object):

    def __init__(self):
        """TODO: to be defined1. """
        pass
        
    def siblings(self):
        pass

    def trio(self):
        pass
