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
from typing import Tuple, Union, List, Dict
import pandas as pd
import numpy as np
import collections
from .util import sort_df

class Annotation(object):
    
    def __init__(self,
                 data : pd.core.frame.DataFrame,
                 unk_val : any = None,
                 field_level : int = 2):
        self.data = data
        self.field_level = field_level
        self.unk_val = unk_val
        self.df = self._parse_data()

    def _parse_data(self) -> pd.core.frame.DataFrame:
        df = self.data
        unk_val = self.unk_val
        
        # annotations = {}
        unk_val = unk_val or "unknown"
        df = df.fillna(unk_val)
        # alleles = []
        allele_col = [col for col in df.columns if 'allele' in col]
        if len(allele_col) == 1:
            # raise Exception('There needs to be exactly one allele column present in this dataset.', df)
            allele_col = allele_col[0]
            df = df.groupby(allele_col).agg(lambda x : 
                '/'.join((sorted(x.replace(unk_val, np.nan)\
                                .dropna().unique().astype(str),
                                key=lambda y: -collections.Counter(x)[y])   
                        ))
            )
            df = sort_df(df.replace('', np.nan).fillna(unk_val))
        return df

        # for attrib in df.columns:
        #     if attrib == 'tce':
        #         null_val = "0"
        #     elif attrib == 'expression':
        #         null_val = "null"
        #     attribs = set(df[attrib])
        #     attribs = [attrib for attrib in attribs
        #                 if attrib != unk_val and attrib != null_val]
        #     if len(attribs) == 1:
        #         common_ann = attribs.pop()
        #     elif len(attribs) == 0:
        #         attribs = set(df[attrib])
        #         if (len(attribs) == 1) and (attribs.pop() == null_val):
        #             common_ann = null_val
        #         else:
        #             common_ann = unk_val
        #     else:
                # m = 
                # if m and ('exon' in )
                # print(attribs, m)
                # if re.match('[ACGT.*]+', attribs[0]) or ('exon' in attrib) or ('intron' in attrib):
                #     common_ann = self.get_cons_seq(gene_feat=attrib,
                #         df=df)
                # else:

                # if True:
                #     # Add alleles with minor features
                #     df_minor = df[~df[attrib].isin(df[attrib].mode())]
                #     if not df_minor.empty and ('CIWD' not in attrib):
                #         alleles += list(df_minor.index)

                    # for ciwd in ['C', 'I', 'WD', 'R', None]:
                    #     if ciwd:
                    #         df_ciwd = df[df['CIWD_TOTAL_combined'].str.contains('^' + ciwd, regex=True, na=False)]
                    #     else:
                    #         df_ciwd = df[df['CIWD_TOTAL_combined'].isnull()]
                    #     if not df_ciwd.empty:
                    #         attribs_ciwd = set(df_ciwd[attrib])
                    #         attribs_ciwd = [attrib for attrib in attribs_ciwd
                    #                                 if attrib != unk_val and attrib != null_val]
                    #         types = set([type(attrib_ciwd) for attrib_ciwd in attribs_ciwd])
                    #         if types and (not bool == types.pop()):
                    #             if len(attribs_ciwd) == 1:
                    #                 common_ann = '~' + str(attribs_ciwd.pop())
                    #             elif (len(df_ciwd[attrib]) == 2) and (len(set(df_ciwd[attrib])) == 2):
                    #                 common_ann = '/'.join(sorted(df[attrib].astype(str)))
                    #             elif len(attribs_ciwd) > 1:
                    #                 common_ann = '?' + str(df[attrib].mode().iloc[0])
                    #             else:
                    #                 common_ann = unk_val
                    #             break
            #     print(df.columns)
            #     print(attrib)
            #     print(attribs)
            #     common_ann = attribs
            # annotations[attrib] = common_ann
        # # alleles = [Allotype(allele, expand=False, ref_data=self.ref_data) for allele in set(alleles)]
        # return annotations

    def calc_common_annotations(self, df : pd.core.frame.DataFrame):
        unk_val = self.unk_val
        unk_val = unk_val or "unknown"
        null_val = None
        alleles = []
        annotations = {}
        df = df.replace('', np.nan)
        df = df.fillna(value=unk_val)
        for attrib in df.columns:
            if attrib == 'tce':
                null_val = "0"
            elif attrib == 'expression':
                null_val = "null"
            elif attrib == 'P2':
                null_val = '*'
            attribs = set(df[attrib])
            attribs = [attrib for attrib in attribs
                        if attrib != unk_val and attrib != null_val]
            if len(attribs) == 1:
                common_ann = attribs.pop()
            elif len(attribs) == 0:
                attribs = set(df[attrib])
                if (len(attribs) == 1) and (attribs.pop() == null_val):
                    common_ann = null_val
                else:
                    common_ann = unk_val
            else:
                # m = 
                # if m and ('exon' in )
                # print(attribs, m)
                # if re.match('[ACGT.*]+', attribs[0]) or ('exon' in attrib) or ('intron' in attrib):
                #     common_ann = self.get_cons_seq(gene_feat=attrib,
                #         df=df)
                # else:
                if True:
                    # Add alleles with minor features
                    df_minor = df[~df[attrib].isin(df[attrib].mode())]
                    if not df_minor.empty and ('CIWD' not in attrib):
                        alleles += list(df_minor.index)

                    for ciwd in ['C', 'I', 'WD', 'R', None]:
                        if ciwd:
                            df_ciwd = df[df['CIWD_TOTAL_combined'].str.contains('^' + ciwd, regex=True, na=False)]
                        else:
                            df_ciwd = df[df['CIWD_TOTAL_combined'].isnull()]
                        if not df_ciwd.empty:
                            attribs_ciwd = set(df_ciwd[attrib])
                            attribs_ciwd = [attrib for attrib in attribs_ciwd
                                                    if attrib != unk_val and attrib != null_val]
                            types = set([type(attrib_ciwd) for attrib_ciwd in attribs_ciwd])
                            if types and (not bool == types.pop()):
                                if len(attribs_ciwd) == 1:
                                    prefix = '~'
                                    attrib_ciwd = attribs_ciwd.pop()
                                    if isinstance(attrib_ciwd, str) and '/' in attrib_ciwd:
                                        prefix = ''
                                    common_ann = prefix + str(attrib_ciwd)
                                elif (len(df_ciwd[attrib]) == 2) and (len(set(df_ciwd[attrib])) == 2):
                                    common_ann = '/'.join(sorted(df[attrib].astype(str)))
                                elif len(attribs_ciwd) > 1:
                                    common_ann = '?' + str(df[attrib].mode().iloc[0])
                                else:
                                    common_ann = unk_val
                                break
                # print(df.columns)
                # print(attrib)
                # print(attribs)
            annotations[attrib] = common_ann
        return annotations

    def __repr__(self):
        return str(self.data)

    def serialize(self):
        # return {feat_name : str(seq) for feat_name, seq in self.feats.items()}
        return self.calc_common_annotations(self.data)
        