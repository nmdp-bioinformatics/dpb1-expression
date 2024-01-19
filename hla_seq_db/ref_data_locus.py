#
# Copyright (c) 2024 NMDP.
#
# This file is part of ExPAT 
# (see https://github.com/nmdp-bioinformatics/dpb1-expression).
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
from hla_seq_db.sequence import Sequence
from hla_seq_db.allotype import Allotype
from hla_seq_db.features import Features
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from .util import get_two_field_allele, calc_hla_class
import pandas as pd
from typing import List, Tuple, Dict, Union
from pyard.ard import ARD
import editdistance
import inflect

class RefDataLocus(object):
    def __init__(self, df : pd.core.frame.DataFrame, locus : str,
            tce : Dict[str, str] = None,
            ciwd : pd.core.frame.DataFrame = None,
            ard : ARD = None,
            g_groups : Dict[str, Dict[str, List[str]]] = None,
            proteins : pd.core.frame.DataFrame = None,
            motifs : Dict[str, Tuple[str, Union[List[int], range]]] = None):
        self.locus = locus
        self.wildcard = None
        self.columns = {}
        self.gene_feat_names = self._calc_features(locus)
        self.ard = ard
        self.tce = tce
        self.proteins = proteins
        self.motifs = motifs #self._calc_motifs()
        # self.ciwd = self._filter_to_locus(ciwd)
        self.ciwd = ciwd
        self.df = self._annotate_data(df, locus, ciwd=ciwd, tce=tce, g_groups=g_groups)
        # print(self.locus, self.df.columns)
    
    # def _calc_motifs(self, motif_name : str) -> Dict[str, SeqFeature]:
    def get_motif(self, motif_name : str) -> SeqFeature:
        # motifs = {}
        # for motif_name, motif in self.motifs.items():
        if self.motifs and motif_name in self.motifs:
            feat_name, positions = self.motifs[motif_name]
            if isinstance(positions, list):
                if len(positions) == 1:
                    location = FeatureLocation(positions[0], positions[0] + 1)
                else:
                    location = CompoundLocation([FeatureLocation(pos, pos + 1) for pos in positions])
            elif isinstance(positions, range):
                location = FeatureLocation(positions.start, positions.stop)
            else:
                location = None
            return SeqFeature(location, type=feat_name)
            # motifs[feat_name] = SeqFeature(location, type=feat_name)
        # return motifs
    
    # def _filter_to_locus(self, pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    #     if isinstance(data, dict):
    #         print(data)
    def _calc_features(self, locus : str, db_type : str = 'gen') -> List[str]:
        """
        Calculates a list of feature names HLA loci.
        """
        if db_type in ['gen', 'nuc']:
            if locus == 'B':
                n_exons = 7
            elif locus in ['DRB1', 'DQB1']:
                n_exons = 6
            else:
                hla_class = calc_hla_class(locus)
                if hla_class == 1:
                    n_exons = 8
                elif hla_class == 2:
                    n_exons = 5
            feature_names = (['utr5'] + 
                            [feature + '_' + str(i)
                                for i in range(1, n_exons + 1) for feature in ['exon', 'intron'] 
                                if not (i == n_exons and feature == 'intron')] +
                            ['utr3'])
            if db_type=='nuc':
                feature_names = [feat for feat in feature_names if 'exon' in feat]
        elif db_type == 'prot':
            feature_names = ['aa']
        return feature_names
    
    def _update_columns(self, col_type : str, cols : List[str]) -> None:
        if isinstance(cols, str):
            cols = [cols]
        if col_type not in self.columns:
            self.columns[col_type] = []
        if not any([col in self.columns[col_type] for col in cols]):
            self.columns[col_type] += cols
    
    def _annotate_data(self, df : pd.core.frame.DataFrame, locus : str,
            ciwd : pd.core.frame.DataFrame = None,
            tce : Dict[str, str] = None,
            g_groups : Dict[str, Dict[str, List[str]]] = None) -> pd.core.frame.DataFrame:
        self._update_columns('gene_feats', list(df.columns))
        if locus == 'DPB1':
            df = self._annotate_ref_data(df,
             ['dpb1_expr', 'g_groups', 'tce'], 
             ciwd=ciwd, tce=tce, g_groups=g_groups,
             loci=['DPB1'])
        elif locus == 'B':
            df = self._annotate_ref_data(df, ['b_leader', 'g_groups'], loci=['B'], ciwd=ciwd, g_groups=g_groups)
        df = self._annotate_ref_data(df, ['CIWD'], ciwd=ciwd, g_groups=g_groups)
        df = self._annotate_seq(df)
        return df

    def _annotate_seq(self, df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
        df_gene_feats = df[self.columns['gene_feats']]
        df['seq'] = df_gene_feats.sum(1).str.replace('.', '', regex=False)
        return df
    
    def _annotate_ref_data(self, df : pd.core.frame.DataFrame,
            modes : List[str], loci : List[str] = None,
            ciwd : pd.core.frame.DataFrame = None,
            tce : Dict[str, str] = None,
            g_groups : Dict[str, Dict[str, List[str]]] = None) -> pd.core.frame.DataFrame:
        """
        Annotates reference data
        """
        if isinstance(modes, str):
            modes = [modes]
        for mode in modes:
            if mode == 'dpb1_expr':
                df = self._annotate_seq_data(df, 'DPB1')
                df = self._assign_ciwd(df, ciwd)
                df = self._assign_tce(df, tce)
                df = self._assign_expression(df, 'DPB1')
            elif mode == 'g_groups':
                # for locus in loci:
                df = self._assign_g_groups(df, g_groups)
            elif mode == 'CIWD':
                # for locus in loci:
                df = self._assign_ciwd(df, ciwd)
            elif mode == 'b_leader':
                self._assign_b_leader(df)
            # elif mode == 'gfe':
            #     for locus in loci:
            #         self._assign_gfe(locus)
        return df

    def _assign_b_leader(self, df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
        df['P2'] = self.proteins['aa'].str[3]
        self._update_columns('gene_ann_feats', 'P2')
        return df
    
    def _annotate_seq_data(self, df : pd.core.frame.DataFrame, 
            locus : str) -> pd.core.frame.DataFrame:
        motifs = self.motifs
        if locus == 'DPB1' and motifs:
            for motif, (gene_feat, positions) in motifs.items():
                df = self._assign_motif(df, motif, gene_feat, positions)
            # df = self._assign_motif(df, 'rs9277534', 'utr3', [791])
            # df = self._assign_motif(df, 'ctcf', 'intron_2', [1954, 1971, 1983])
            # df = self._assign_motif(df, 'str', 'intron_2', range(3893, 3960),
                    # remove_spacers=False, remove_unknown=False)
        self._update_columns('gene_ann_feats', list(motifs.keys()))
        return df

    def _assign_motif(self, df : pd.core.frame.DataFrame,
                     label : str, gene_feature : str, 
                     positions : Tuple[List[int], range],
                     remove_spacers : bool = False, 
                     remove_unknown : bool = False) -> pd.core.frame.DataFrame:
        """
        Extracts a motif from locus, gene feature, and positions provided
        and assigns it to a labelled column on RefData
        """
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
    
    def _assign_tce(self, df : pd.core.frame.DataFrame,
            tce : Dict[str, str]) -> pd.core.frame.DataFrame:
        header = 'tce'
        if header not in df.columns:
            tce_df = pd.DataFrame.from_dict(tce, orient='index', columns=['tce'])
            tce_df = df.join(tce_df, how='left')
            df[header] = tce_df[header]
        self._update_columns('ann_feats', header)
        return df
    
    def _assign_ciwd(self, df : pd.core.frame.DataFrame,
            df_ciwd : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
        """
        Downloads the CIWD 3.0.0 table from the NMDP Bioinformatics repository
        as a Pandas DataFrame.
        :return: CIWD table as Pandas Dataframe with 'allele' and 'CIWD' status columns.
        :rtype: pd.DataFrame
        """
        header = 'CIWD_TOTAL_combined'
        if header not in df:
            for i in range(1, 5):
                if i == 1:
                    df = df.join(df_ciwd, how='left')
                else:
                    num2word = {2 : 'two', 3 : 'three', 4 : 'four'}
                    num_word = num2word[i]
                    df = self._assign_field_typings(df, i)
                    df = pd.merge(df, 
                        df_ciwd, how='left',
                        suffixes=("", "_{}_field".format(num_word)),
                        left_on='allele_{}_field'.format(num_word), right_index=True)
            df[header] = df['CIWD_TOTAL_four_field']\
                .fillna(df['CIWD_TOTAL_three_field'])\
                .fillna(df['CIWD_TOTAL_two_field'])
        self._update_columns('ann_feats', header)
        return df
    
    def _assign_field_typings(self, df : pd.core.frame.DataFrame,
            num_fields : int):
        num2word = {1 : 'one', 2 : 'two', 3 : 'three', 4 : 'four'}
        num_word = num2word[num_fields]
        header = 'allele_{}_field'.format(num_word)
        df['allele'] = df.index
        if header not in df:
            # allele_typing = df.index.str.split(':').str[:num_fields].str.join(':')
            allele_typing = (df['allele'].str.split(':').str[:num_fields].str.join(':').str.extract('([A-Z]+[0-9]*\*[0-9:]+)')[0] +
                    df['allele'].str.extract('([A-Z]$)')[0].fillna(''))
            # allele_typing = (df.index.str.split(':').str[:num_fields].str.join(':').str.extract('([A-Z]+[0-9]*\*[0-9:]+)')[0]  +
            #                     df.index.str.extract('([A-Z]$)')[0].fillna(''))
            # print(allele_typing)
            df[header] = allele_typing
        self._update_columns('allele', header)
        return df

    def _assign_expression(self, df : pd.core.frame.DataFrame, locus : str) -> pd.core.frame.DataFrame:
        header_expr = 'expression'
        if locus == 'DPB1':
            expr_motifs = {'low' : {'exon_3_motif' : 'GTTGTCT', 'rs9277534' : 'A'},
                        'high' : {'exon_3_motif' : 'ACCACTC', 'rs9277534' : 'G'},
                        'unknown' : {'exon_3_motif' : '*******'}}
            
            for expr, motifs in expr_motifs.items():
                for header, motif in motifs.items():
                    df.loc[df[header] == motif, header_expr] = expr
            df.loc[:, 'expr_n_diffs'] = 0
            for i, row in df[df[header_expr].isnull()].iterrows():
                exon_3_motif = row['exon_3_motif']
                if isinstance(exon_3_motif, str):
                    for expr, motifs in expr_motifs.items():
                        diff = editdistance.eval(motifs['exon_3_motif'], exon_3_motif)
                        if diff <= 2:
                            df.loc[i, header_expr] = expr
                            df.loc[i, 'expr_n_diffs'] = diff
            df = self._assign_experimental_alleles(df, header_expr, ["DPB1*02:01", "DPB1*02:02", "DPB1*04:01", "DPB1*04:02",
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
            df_indices = df.index.str.contains(code + '$', regex = True)
            df.loc[df_indices, header_expr] = expr
        self._update_columns('ann_feats', header_expr)
        self._update_columns('ann_feats', 'expr_n_diffs')
        return df

    def _assign_experimental_alleles(self, 
            df : pd.core.frame.DataFrame,
            label : str, alleles : List[str]) -> pd.core.frame.DataFrame:
        label += '_experimental'
        alleles = [Allotype(allele, expand=False) for allele in alleles]
        loci = {}
        for allele in alleles:
            if allele.locus not in loci:
                loci[allele.locus] = []
            loci[allele.locus].append(str(allele))
        for locus, alleles in loci.items():
            df[label] = 'no'
            df.loc[df['allele_two_field'].isin(alleles), label] = 'yes'
        self._update_columns('ann_feats', label)
        return df

    def _assign_g_groups(self, df : pd.core.frame.DataFrame,
            g_groups : Dict[str, Dict[str, List[str]]]) -> pd.core.frame.DataFrame:
        header = 'g_group'
        df_g_groups = pd.DataFrame([[g_group, allele] 
                for g_group, alleles in g_groups.items()
                for allele in alleles], columns=[header, 'allele'])\
                    .set_index('allele')
        df_allele_g_groups = df.join(df_g_groups, how='left')
        df[header] = df_allele_g_groups[header]
        return df
    
    def filter_to_col_type(self, col_type : str) -> pd.core.frame.DataFrame:
        return self.df[self.columns[col_type]]

    def filter_to_allotype(self, allotype : Allotype) -> pd.core.frame.DataFrame:
        if 'XX:XX' in allotype:
            return self.df
        # if allotype.resolution in ['allelic']:
        #     regex = allotype.typing.replace('*', '\*')
        # else:
        alleles = allotype.get_potential_alleles()
        regex = '|'.join(["(?:%s)" % allele.typing.replace('*', '\*') for allele in alleles])
        if not regex:
            raise Exception('This allotype does not appear to contain any IMGT-defined alleles.', allotype.typing)
        return self.df[self.df.index.str.contains('^(?:%s)'  % regex)]
    
    def get_features(self, allotype : Allotype, feats: Union[str, List[str]] = None,
                field_level : int = 2,
                gene_feats : bool = False, consensus : bool = False,
                sire : str = None, verbose : bool = False) -> Features:
        """
        Obtains a DataFrame subset that contains the provided allele.
        Sorts based on CIWD.
        :param str allele: The allele name
        :return: pd.DataFrame
        """
        is_wildcard = 'XX:XX' in allotype or 'XX:XX' in allotype.typing
        if is_wildcard and self.wildcard:
            return self.wildcard
            # self.wildcard = self.generate_cons_gene_feats(verbose=verbose)
            # return self.wildcard
        if feats and isinstance(feats, str):
            feats = [feats]
        df = self.filter_to_allotype(allotype)
        allele_name_cols = []
        if field_level:
            p = inflect.engine()
            col = 'allele_{}_field'.format(p.number_to_words(field_level))
            if col not in df:
                raise Exception('This field level is not supported', field_level)
            allele_name_cols = [col]
        data = {col_type : df[[col for col in cols 
                if (feats == None) or
                   ('gene_feats' in feats) or
                   (col in feats) or 
                   ('CIWD' in col)
                   ] + allele_name_cols]
                for col_type, cols in self.columns.items() if col_type != 'allele'}
        # print([df.columns for df in data.values()])
        feats = Features(data, wildcard=is_wildcard)
        if is_wildcard:
            self.wildcard = feats
        return feats

    def __len__(self):
        return len(self.df)