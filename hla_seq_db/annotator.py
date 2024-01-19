#
# Copyright (c) 2024 NMDP.
#
# This file is part of ExPAT 
# (see https://github.com/nmdp-bioinformatics/hla-db).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
from .ref_data import RefData
from typing import List, Optional, Tuple
import pandas as pd
from py2neo import Graph
from .allotype import Allotype
import editdistance

class Annotator(object):

    def __init__(self, ref_data : RefData, loci : List[str]):
        self.ref_data = ref_data
        self.loci = loci

    def assign_motif(self, label : str, locus : str, gene_feature : str, 
                     positions : Tuple[List[int], range],
                     df : pd.core.frame.DataFrame = None,
                     remove_spacers : bool = False, 
                     remove_unknown : bool = False) -> pd.core.frame.DataFrame:
        """
        Extracts a motif from locus, gene feature, and positions provided
        and assigns it to a labelled column on RefData
        """
        if not isinstance(df, pd.core.frame.DataFrame):
            df = self.ref_data.alleles[locus]
        if gene_feature in df:
            df[gene_feature] = df[gene_feature].astype(str).str.replace(' ', '')
            #.str.replace(' ', '')
            col = pd.concat([df[gene_feature].str[pos]
                    for pos in positions], axis=1)\
                    .apply(lambda row: ''.join(row), axis=1)
            if remove_spacers:
                col = col.str.replace('.', '', regex=False)
            if remove_unknown:
                col = col.str.replace('*', '', regex=False)
            df[label] = col
        return df
    
    def assign_g_groups(self, locus : str) -> None:
        header = 'g_group'
        df_alleles = self.ref_data.alleles[locus]
        if header not in df_alleles:
            g_groups = self.ref_data.g_groups[locus]
            df_g_groups = pd.DataFrame([[g_group, allele] 
                    for g_group, alleles in g_groups.items()
                    for allele in alleles], columns=[header, 'allele'])\
                        .set_index('allele')
            df_allele_g_groups = df_alleles.join(df_g_groups, how='left')
            df_alleles[header] = df_allele_g_groups[header]

    def assign_tce(self) -> None:
        if self.ref_data.tce:
            locus = 'DPB1'
            header = 'tce'
            if locus in self.loci:
                df_alleles = self.ref_data.alleles[locus]
                if header not in df_alleles.columns:
                    tce_df = pd.DataFrame.from_dict(self.ref_data.tce, orient='index', columns=['tce'])
                    tce_df = df_alleles.join(tce_df, how='left')
                    df_alleles[header] = tce_df[header]

    def assign_ciwd(self, locus : str) -> None:
        """
        Downloads the CIWD 3.0.0 table from the NMDP Bioinformatics repository
        as a Pandas DataFrame.
        :return: CIWD table as Pandas Dataframe with 'allele' and 'CIWD' status columns.
        :rtype: pd.DataFrame
        """
        df_alleles = self.ref_data.alleles[locus]
        if 'CIWD_TOTAL' not in df_alleles:
            for i in range(1, 5):
                if i == 1:
                    df_ciwd = df_alleles.join(self.ref_data.ciwd, how='left')
                else:
                    num2word = {2 : 'two', 3 : 'three', 4 : 'four'}
                    num_word = num2word[i]
                    self.assign_field_typings(locus, i)
                    df_ciwd = pd.merge(df_alleles, 
                        self.ref_data.ciwd, how='left',
                        suffixes=("", "_{}_field".format(num_word)),
                        left_on='allele_{}_field'.format(num_word), right_index=True)
                self.ref_data.alleles[locus] = df_alleles = df_ciwd
        df_alleles['CIWD_TOTAL_combined'] = df_alleles['CIWD_TOTAL_four_field']\
            .fillna(df_alleles['CIWD_TOTAL_three_field'])\
            .fillna(df_alleles['CIWD_TOTAL_two_field'])
    
    def assign_field_typings(self, locus : str, num_fields : int):
        df_alleles = self.ref_data.alleles[locus]
        num2word = {1 : 'one', 2 : 'two', 3 : 'three', 4 : 'four'}
        num_word = num2word[num_fields]
        header = 'allele_{}_field'.format(num_word)
        df_alleles['allele'] = df_alleles.index
        if header not in df_alleles:
            # allele_typing = df_alleles.index.str.split(':').str[:num_fields].str.join(':')
            allele_typing = (df_alleles['allele'].str.split(':').str[:num_fields].str.join(':').str.extract('([A-Z]+[0-9]*\*[0-9:]+)')[0] +
                    df_alleles['allele'].str.extract('([A-Z]$)')[0].fillna(''))
            # allele_typing = (df_alleles.index.str.split(':').str[:num_fields].str.join(':').str.extract('([A-Z]+[0-9]*\*[0-9:]+)')[0]  +
            #                     df_alleles.index.str.extract('([A-Z]$)')[0].fillna(''))
            # print(allele_typing)
            self.ref_data.alleles[locus][header] = allele_typing

    def assign_expression(self, locus : str) -> None:
        df_alleles = self.ref_data.alleles[locus]
        if locus == 'DPB1':
            expr_motifs = {'low' : {'exon_3_motif' : 'GTTGTCT', 'rs9277534' : 'A'},
                        'high' : {'exon_3_motif' : 'ACCACTC', 'rs9277534' : 'G'},
                        'unknown' : {'exon_3_motif' : '*******'}}
            
            for expr, motifs in expr_motifs.items():
                for header, motif in motifs.items():
                    df_alleles.loc[df_alleles[header] == motif, 'expression'] = expr
            df_alleles.loc[:, 'expr_n_diffs'] = 0
            for i, row in df_alleles[df_alleles['expression'].isnull()].iterrows():
                exon_3_motif = row['exon_3_motif']
                if isinstance(exon_3_motif, str):
                    for expr, motifs in expr_motifs.items():
                        diff = editdistance.eval(motifs['exon_3_motif'], exon_3_motif)
                        if diff <= 2:
                            df_alleles.loc[i, 'expression'] = expr
                            df_alleles.loc[i, 'expr_n_diffs'] = diff
            self.assign_experimental_alleles('expression', ["DPB1*02:01", "DPB1*02:02", "DPB1*04:01", "DPB1*04:02",
                    "DPB1*17:01", "DPB1*23:01", "DPB1*40:01", "DPB1*46:01",
                    "DPB1*55:01", "DPB1*71:01", "DPB1*94:01", "DPB1*105:01",
                    "DPB1*128:01", "DPB1*01:01", "DPB1*05:01", "DPB1*11:01",
                        "DPB1*13:01", "DPB1*15:01", "DPB1*18:01", "DPB1*19:01",
                        "DPB1*85:01", "DPB1*03:01", "DPB1*06:01", "DPB1*09:01",
                        "DPB1*10:01", "DPB1*14:01", "DPB1*16:01", "DPB1*20:01"])
        mapping = {'Q' : 'questionable', 
                    'N' : 'unknown',
                    'L' : 'low',
                    'S' : 'secreted',
                    'C' : 'cytoplasm',
                    'A' : 'aberrant'}
        for code, expr in mapping.items():
            df_indices = df_alleles.index.str.contains(code + '$', regex = True)
            df_alleles.loc[df_indices, 'expression'] = expr
        # print(df_alleles['expression'].value_counts())

    def assign_experimental_alleles(self, label : str, alleles : List[str]) -> None:
        alleles = [Allotype(allele, expand=False) for allele in alleles]
        loci = {}
        for allele in alleles:
            if allele.locus not in loci:
                loci[allele.locus] = []
            loci[allele.locus].append(str(allele))
        for locus, alleles in loci.items():
            df_alleles = self.ref_data.alleles[locus]
            df_alleles[label + '_experimental'] = 'no'
            df_alleles.loc[df_alleles['allele_two_field'].isin(alleles), label + '_experimental'] = 'yes'

    def assign_b_leader(self) -> None:
        df_proteins = self.ref_data.proteins_aa['B']
        df_alleles = self.ref_data.alleles['B']
        df_alleles['P2'] = df_proteins['aa'].str[3]
    
    def annotate_seq_data(self, locus : str = 'DPB1',
            df : pd.core.frame.DataFrame = None) -> pd.core.frame.DataFrame:
        if locus == 'DPB1':
            self.assign_motif('exon_3_motif', 'DPB1', 'exon_3', [9, 16, 43, 80, 228, 236, 265],
                df=df)
            self.assign_motif('rs9277534', 'DPB1', 'utr3', [791],
                df=df)
            self.assign_motif('ctcf', 'DPB1', 'intron_2', [1954, 1971, 1983],
                df=df)
            self.assign_motif('str', 'DPB1', 'intron_2', range(3893, 3960),
                    remove_spacers=False, remove_unknown=False, df=df)
        return df

    def annotate_ref_data(self, modes : List[str], loci : List[str] = None):
        """
        Annotates reference data
        """
        if not loci:
            loci = self.loci
        if isinstance(modes, str):
            modes = [modes]
        for mode in modes:
            if mode == 'dpb1_expr':
                self.annotate_seq_data('DPB1')
                self.assign_tce()
                self.assign_ciwd('DPB1')
                self.assign_expression('DPB1')
            elif mode == 'tce':
                self.assign_tce()
            elif mode == 'g_groups':
                for locus in loci:
                    self.assign_g_groups(locus)
            elif mode == 'CIWD':
                for locus in loci:
                    self.assign_ciwd(locus)
            elif mode == 'b_leader':
                self.assign_b_leader()
            # elif mode == 'gfe':
            #     for locus in loci:
            #         self.assign_gfe(locus)
