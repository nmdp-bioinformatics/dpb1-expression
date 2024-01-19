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
from .annotator import Annotator
from .sequence_match import SeqMatch
from .allotype import Allotype
from typing import List, Dict, Union, Tuple
import pandas as pd
import numpy as np
from collections import Counter
from .sequence import Sequence
from hla_seq_db.sequence_feat_searched import SeqFeatSearched
import re
from pyard import ARD

class Extractor(object):
    def __init__(self, annotator : Annotator = None,
            ard : ARD = None):
        # self.ref_data = ref_data
        self.annotator = annotator
        self.ard = ard
        self.wildcards = {}

    # def annotate_allotype(self, allotype : Allotype, feats : List[str] = None,
    #             field_level : int = None, sire = None) -> Allotype:
    #     # if isinstance(allotype, str):
    #     #     allotype = Allotype(allotype, ref_data=self.ref_data)
    #     df = self.extract(str(allotype), field_level=field_level)
    #     allele_col = None
    #     if field_level:
    #         num2word = {2 : 'two', 3 : 'three', 4 : 'four'}
    #         allele_col = "allele_{}_field".format(num2word[field_level])
    #     ciwd_header = 'CIWD_{}_combined'.format(sire or 'TOTAL')
    #     cols = [ciwd_header]
    #     if feats:
    #         cols += feats
    #     if allele_col:
    #         cols += [allele_col]
    #     df_alleles = self.group_and_consolidate(df[cols],
    #                                         allele_col)
    #     for allele in (allotype.alleles_hi_res 
    #                     if field_level == 2
    #                     else allotype.get_potential_alleles()):
    #         allele.annotation = df_alleles.loc[allele.typing].to_dict()
    #     allotype.alleles_hi_res = self.sort_alleles(allotype.alleles_hi_res)
    #     if feats:
    #         # TODO: Standardize this
    #         unk_val = None
    #         if "P2" in feats:
    #             unk_val = "*"
    #         annotation, alleles = self.calc_common_annotation(df[cols], cols, unk_val=unk_val)
    #         allotype.annotation = annotation
    #         for allele in alleles:
    #             df_annotation = None
    #             if allele.typing in df_alleles.index:
    #                 df_annotation = df_alleles
    #             elif allele.typing in df.index:
    #                 df_annotation = df[cols]
    #             if isinstance(df_annotation, pd.core.frame.DataFrame):
    #                 allele.annotation = df_annotation.loc[allele.typing].fillna('unknown').to_dict()
    #         allotype.alleles = self.sort_alleles(alleles)
    #         diffs = {}
    #         feats_annotated = {}
    #         df_seq = pd.DataFrame.from_dict({feat_name : [feat_seq.seq_two.raw_seq()]
    #                 for feat_name, feat_seq in allotype.gene_feats['sequences'].features.items()})
    #         df_seq = self.annotator.annotate_seq_data('DPB1', df_seq)
    #         df_seq_ref = self.extract(allotype.typing, feats=feats, consensus=True) # Takes a while
    #         for feat_name in feats:
    #             if feat_name in df_seq:
    #                 feat = Sequence(df_seq[feat_name], feat_name)
    #                 feat_ref = df_seq_ref[feat_name]
    #                 feats_annotated[feat_name] = SeqFeatSearched(feat_ref, feat)
    #                 print(feat_name)
    #             else:
    #                 feats_annotated[feat_name] = str(annotation[feat_name])
    #         allotype.set_features(feats_annotated)
    #         allotype.calc_and_set_snps()
    #     return allotype
        
    def sort_alleles(self, alleles : List[Allotype]) -> List[Allotype]:
        """
        Sorts list of Allotype objects based of most common CIWD code, 
        then HLA field nomenclature (first field, second, third, then fourth).
        """
        ciwd_order = ['C', 'I', 'WD', 'R', 'unknown']
        cols_fields = ['f1', 'f2', 'f3', 'f4']
        header_ciwd_highest = 'ciwd_highest_index'
        df_sorter = pd.DataFrame([[allele] + [int(f) if f else 0 for f in allele.fields] +
                           [allele.annotation['CIWD_TOTAL_combined']]
                            for allele in alleles], 
                         columns=['allele'] + cols_fields + ['ciwd'])
        df_sorter[header_ciwd_highest] = df_sorter['ciwd'].apply(lambda x : 
                                [i for i, ciwd in enumerate(ciwd_order)
                                if ciwd in x.split('/') or ciwd == 'unknown'][0])
        df_sorter = df_sorter.sort_values([header_ciwd_highest] + cols_fields)
        return list(df_sorter['allele'])

    # def extract(self, typing : str, feats: Union[str, List[str]] = None,
    #             field_level : int = None,
    #             gene_feats : bool = False, consensus : bool = False,
    #             sire : str = None, verbose : bool = False) -> Union[pd.core.frame.DataFrame,
    #                                                                  Sequence,
    #                                                                  Dict[str, str]]:
    #     """
    #     Obtains a DataFrame subset that contains the provided allele.
    #     Sorts based on CIWD.
    #     :param str allele: The allele name
    #     :return: pd.DataFrame
    #     """
    #     if 'XX:XX' in typing:
    #         locus = typing.split('*')[0]
    #         if locus in self.wildcards:
    #             return self.wildcards[locus]
    #         df = self.ref_data.alleles[locus].df
    #         result = self.generate_cons_gene_feats(df, typing=typing, verbose=verbose)
    #         if locus not in self.wildcards:
    #             self.wildcards[locus] = result
    #     else:
    #         if feats and isinstance(feats, str):
    #             feats = [feats]
    #         if isinstance(typing, Allotype):
    #             allele = typing
    #         else:
    #             allele = Allotype(str(typing), ref_data=self.ref_data, ard=self.ard)
    #         if allele.resolution in ['allelic']:
    #             regex = allele.typing.replace('*', '\*')
    #         else:
    #             # if field_level == 2:
    #             #     alleles = allele.alleles_hi_res
    #             # else:
    #             alleles = allele.get_potential_alleles()
    #             regex = '|'.join(["(?:%s)" % allele.typing.replace('*', '\*') for allele in alleles])
    #         alleles = self.ref_data.alleles[allele.locus]
    #         if not regex:
    #             raise Exception('This allotype does not appear to contain any IMGT-defined alleles.', allele.typing)
    #         df = alleles[alleles.index.str.contains('^(?:%s)'  % regex)].copy()
    #     if gene_feats:
    #         df = df[self.ref_data._calc_features(allele.locus)] #.to_dict('records')
    #         return self.generate_cons_gene_feats(df, typing=typing, verbose=verbose)
    #     if consensus:
    #         # if feats:
    #         #     return df[feats]
    #         df = self.generate_cons_gene_feats(df, typing=typing, verbose=verbose)
    #     if feats:
    #         if len(feats) == 1 and consensus:
    #             return Sequence(df[feats[0]], name=feats[0])
    #         return df
    #     ciwd_header = 'CIWD_{}_combined'.format(sire or 'TOTAL')
    #     if ciwd_header in df:
    #         df['sort_col'] = pd.Categorical(
    #             df[ciwd_header],
    #             categories=['C', 'I', 'WD'],
    #             ordered=True
    #         )
    #         # print(df[[feat, 'sort_col', allele_col]].groupby(allele_col))
    #         df = df.sort_values('sort_col').drop('sort_col', axis=1)
    #     return df
        
    def generate_cons_gene_feats(self, 
            df : pd.core.frame.DataFrame, typing : str = None,
            verbose : bool = False) -> pd.core.frame.DataFrame:
        """
        Iterates through the gene feature columns of a dataframe to
        create the consensus sequence for each gene feature. Ignores
        empty typing '*' if other alleles have typing. If there is more than
        one nucleotide at a given position, then an "X" is inserted.
        """
        result = {}
        for i, gene_feat in enumerate(df.columns):
            if ('intron' in gene_feat) or ('exon' in gene_feat) or ('utr' in gene_feat) or ('ctcf' in gene_feat) or ('str' in gene_feat):
                seqs = df[gene_feat]
                seq = Sequence(seqs, name=gene_feat, index=i)
                if seq.seq:
                    result[gene_feat] = seq
                # seqs = seqs.replace(r'^[.*]+$', np.nan, regex=True).dropna()
                # seq = ""
        #         for i in range(len(seqs[0])):
        #             nts = set(seqs.str[i])
        #             for char in ['*', np.nan]:
        #                 if char in nts:
        #                     nts.remove(char)
        #             if len(nts) > 1:
        #                 nt = 'X'
        #                 # nt = '[{}]'.format(''.join(sorted(nts)))
        #             else:
        #                 nt = nts.pop()
        #             seq += nt
        #         result[gene_feat] = [seq]
        # df = pd.DataFrame(result)
        # if typing:
        #     df.index = [typing]
        return result

    def group_and_consolidate(self, df : pd.core.frame.DataFrame,
                             group_col : str = None) -> pd.core.frame.DataFrame:
        unk_val = "unknown"
        df_out = df.fillna(unk_val).groupby(group_col or (lambda x : True))\
            .agg(lambda x: '/'.join(list(dict.fromkeys(
                sorted(filter(lambda y : y != unk_val, x),
                    key=lambda z: -Counter(filter(lambda y : y != unk_val, z))[z])
                                                      )
                                         )
                                    )
                 ).replace('', unk_val)
        if group_col:
            df_out = df_out.reset_index()
        df_out['field_one'] = df_out['allele_two_field'].str.split('\*').str[1].str.split(':')\
                                        .str[0].astype(int) #.str.split('*')) #.str[1])
        df_out['field_two'] = df_out['allele_two_field'].str.split('\*').str[1].str.split(':')\
                                        .str[1].str.replace('[A-Z]', '', regex=True).astype(int)
        df_out['CIWD_one'] = df_out['CIWD_TOTAL_combined'].str.split('/').str[0]
        df_out = df_out.sort_values(['CIWD_one', 'field_one', 'field_two'])\
                .drop(['CIWD_one', 'field_one', 'field_two'], axis=1)\
                .set_index(group_col)
        return df_out

        
        
    # def calc_common_annotation(self, df : pd.core.frame.DataFrame,
    #         feats : List[str],
    #         unk_val : str = None) -> Tuple[Dict[str, str], Union[List[Allotype], None]]:
    #     annotations = {}
    #     unk_val = unk_val or "unknown"
    #     df = df.fillna(unk_val)
    #     null_val = None
    #     alleles = []
    #     for attrib in feats:
    #         if attrib == 'tce':
    #             null_val = "0"
    #         elif attrib == 'expression':
    #             null_val = "null"
    #         attribs = set(df[attrib])
    #         attribs = [attrib for attrib in attribs
    #                     if attrib != unk_val and attrib != null_val]
    #         if len(attribs) == 1:
    #             common_ann = attribs.pop()
    #         elif len(attribs) == 0:
    #             attribs = set(df[attrib])
    #             if (len(attribs) == 1) and (attribs.pop() == null_val):
    #                 common_ann = null_val
    #             else:
    #                 common_ann = unk_val
    #         else:
    #             # m = 
    #             # if m and ('exon' in )
    #             # print(attribs, m)
    #             if re.match('[ACGT.*]+', attribs[0]) or ('exon' in attrib) or ('intron' in attrib):
    #                 common_ann = self.ref_data.get_cons_seq(gene_feat=attrib,
    #                     df=df)
    #             else:
    #                 # Add alleles with minor features
    #                 df_minor = df[~df[attrib].isin(df[attrib].mode())]
    #                 if not df_minor.empty and ('CIWD' not in attrib):
    #                     alleles += list(df_minor.index)

    #                 for ciwd in ['C', 'I', 'WD', 'R', None]:
    #                     if ciwd:
    #                         df_ciwd = df[df['CIWD_TOTAL_combined'].str.contains('^' + ciwd, regex=True, na=False)]
    #                     else:
    #                         df_ciwd = df[df['CIWD_TOTAL_combined'].isnull()]
    #                     if not df_ciwd.empty:
    #                         attribs_ciwd = set(df_ciwd[attrib])
    #                         attribs_ciwd = [attrib for attrib in attribs_ciwd
    #                                                 if attrib != unk_val and attrib != null_val]
    #                         if not bool == set([type(attrib_ciwd) for attrib_ciwd in attribs_ciwd]).pop():
    #                             if len(attribs_ciwd) == 1:
    #                                 common_ann = '~' + attribs_ciwd.pop()
    #                             elif (len(df_ciwd[attrib]) == 2) and (len(set(df_ciwd[attrib])) == 2):
    #                                 common_ann = '/'.join(sorted(df[attrib]))
    #                             elif len(attribs_ciwd) > 1:
    #                                 common_ann = '?' + df[attrib].mode().iloc[0]
    #                             else:
    #                                 common_ann = unk_val
    #                             break
    #         annotations[attrib] = common_ann
    #     alleles = [Allotype(allele, expand=False, ref_data=self.ref_data) for allele in set(alleles)]
    #     return annotations, alleles