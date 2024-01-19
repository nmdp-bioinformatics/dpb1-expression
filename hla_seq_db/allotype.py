#
# Copyright (c) 2024 NMDP.
#
# This file is part of DPB1 Expression 
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
import re
import requests
import json
from .gene_features import GeneFeatures
from .sequence_features import SeqFeatures
from .features import Features
from hla_seq_db.sequence_search import SeqSearch
from hla_seq_db.sequence_feat_searched import SeqFeatSearched
import hla_seq_db
# from .annotation import Annotation
from .sequence import Sequence
from typing import List, Dict, Tuple, Union
from pyard import ARD
import pandas as pd
import time
import regex
from .dataset import Dataset

class Allotype(object):
    
    def __init__(self, hla: Union[str, pd.core.frame.DataFrame], ref_data = None, ard : ARD = None,
                 seq_diffs : bool = False,
                 seq : Sequence = None, db_version : str = None,
                 expand : bool = True, verbose : bool = False) -> None:
        """
        Represents an allele. In this context, only DPB1 alleles are allowed.
        Valid allele names start with DPB1* followed by two multi-digit integers
        separated by colons.
        """
        self.hla = hla
        self.seq = seq
        self.feats = self.feats_searched = self.feats_expanded = None
        self.ref_data = ref_data
        self.alleles_called = self.alleles = None
        self.p_group = self.alleles_hi_res = self.fields = self.resolution = self.var_expression = self.g_group = None
        self.ard = ard
        self.expand = expand
        self.typing, self.seq = self._process_hla()
        # self.feats = feats
        self.verbose = verbose
        self.db_version = db_version
        self.seq_diffs = seq_diffs
        self.gene_feats = {}
        self.snps = None
        self.annotation = self.annotation_seq = self.annotation_diffs = None
        self.gene_feat_diffs = []
        self.locus = self._calc_locus()
        self.is_gl, self.alleles = self._is_gl()
        self.wildcard = False
        self._get_info()
        self.exact_alleles_called = True
        # self.alleles_called =  self.call_alleles() if self.seq else None
        if not self.is_gl and expand:
            self.alleles = self.get_potential_alleles()
        self.hla_class = self._calc_hla_class()
        if expand:
            self.alleles_hi_res = self._calc_alleles_hi_res()
            self.p_group = self.calc_group('p')
            self.g_group = self.calc_group('g')
            self.allele_family = self._calc_allele_family()
            # self.features = self.get_features()
        self.tce = self.expr_level = self.expr_annot_type = None
        
        self.matched = False
        self.potent_matched = False
    #     # if self.feats or self.seq:
        # self.features = self.get_features()
    #     if self.seq and self.seq_diffs:
    #         self.features_searched = self._parse_gene_features()
    #         self.calc_and_set_snps()

    def get_features(self, feats : List[str] = [], field_level : int = 2) -> Features:
        # feats = self.feats or feats
        if isinstance(feats, str):
            feats = [feats]
        if feats and ('exon_3_motif' in feats) and ('exon_3' not in feats):
            feats += ['exon_3']
        # if feats and 'gene_feats' not in feats:
        #     feats += ['gene_feats']
        # if not self.features and not self.alleles_called:
        features = self.ref_data.get_features(self, feats=feats, field_level=field_level)
        self.feats = features
        df = self.feats.anns.df
        if len(df) > 1:
            if feats and ('CIWD_TOTAL_combined' not in feats):
                feats += ['CIWD_TOTAL_combined']
                df = df[[feat for feat in feats if feat in df.columns]]
            self.feats_expanded = df.reset_index().rename(columns={'index' :  'typing'}).to_dict(orient='records')
        return self.feats

    # def set_features_searched(self, seq_features : SeqFeatures) -> None:
    #     self.features_searched = seq_features

    def call_alleles(self) -> str:
        df_locus_seqs = self.ref_data.alleles[self.locus].df
        alleles = df_locus_seqs[df_locus_seqs['seq'].str.contains(self.seq)].index
        self.alleles_called = '/'.join(alleles)
        return self.alleles_called

    def call_similar_alleles(self) -> str:
        df = self.ref_data.alleles[self.locus].filter_to_col_type('gene_feats')
        df = df[df.columns[df.columns.str.contains('exon')]]
        alleles = set()
        for col in df.columns:
            df_result = df[df[col].apply(lambda x : 
                                        bool(regex.search(x.replace('\.', '').replace('*', ''),
                                                        self.seq)))]
            if not df_result.empty:
                alleles.update(df_result.index)
        self.alleles_called = '/'.join(alleles)
        self.exact_alleles_called = False
        return self.alleles_called

    def parse_gene_features(self, feats : List[str] = None, verbose : bool = False) -> SeqSearch:
        # feats = self.feats or feats
        seq_searches = None
        verbose = self.verbose or verbose
        for seq in self.seq.split('|'):
            if True:
            # try:
                if self.verbose:
                    print('Parsing gene features for', self.typing)
                    start_time = time.time()
                seq_search = SeqSearch(seq, self.feats.seqs,
                                locus=self.locus, feats=feats, exact=bool(self.alleles_called),
                                ref_data=self.ref_data, verbose=verbose)
                if self.verbose:
                    print('Finished parsing after {0:.2f}s\n'.format(time.time() - start_time))
                if seq_searches:
                    seq_searches += seq_search
                else:
                    seq_searches = seq_search
            # except Exception as e:
            #     print(self.typing, str(e))
        self.feats_searched = seq_searches
        return self.feats_searched
    
    def _process_hla(self) -> Tuple[str, str]:
        seq = self.seq
        typing = None
        if isinstance(self.hla, pd.core.frame.DataFrame):
            df = self.hla
            if 'ALLOTYPE' not in df:
                ds = Dataset(df)
                # gene_feats = ds.get_dataframe(headers_allele=ds.headers_allele, headers_type=['GENE_FEAT'], headers_gene_feats=['EXON_2'])
                loci = ds.headers['loci']
                locus = None
                if len(loci) != 1:
                    raise Exception('There needs to be exactly one locus present in the allotype input.')
                locus = loci[0]
                df_typing = ds.get_dataframe(headers_allele=ds.headers_allele, headers_gene_feats=[])
                if not typing and (len(df_typing) == 1) and (len(df_typing.columns) == 1):
                    typing = df_typing.iloc[0][0]
                if locus not in typing:
                    typing = locus + '*' + typing.strip()
                self.typing = typing
                self.locus = locus or self._calc_locus()
                feats = {}
                for header_gene_feat in ds.headers_gene_feats:
                    df_gene_feats = ds.get_dataframe(headers_gene_feats=header_gene_feat,
                                headers_type='GENE_FEAT')
                    df_gene_feats = df_gene_feats[df_gene_feats.columns[df_gene_feats.columns.str.contains('GENE_FEAT')]]
                    if not df_gene_feats.empty:
                        if not ((len(df_gene_feats) == 1) and (len(df_gene_feats.columns) == 1)):
                            raise Exception('There were multiple gene features and columns detected.', df_gene_feats)
                        gene_feat = df_gene_feats.iloc[0][0]
                        if not self.feats and self.expand:
                            self.get_features(feats='gene_feats')
                        name = header_gene_feat.lower()
                        seq_feat = SeqFeatSearched(self.feats.seqs.feats[name],
                                    gene_feat, name=name)
                        feats[name] = seq_feat
                if self.feats:
                    feats = {'seqs' : SeqFeatures(feats)}
                    self.feats_searched = SeqSearch(seq, self.feats.seqs,
                        seq_features=feats)
                seq = None
            else:
                if ('ALLOTYPE' not in df) and ('SEQ' not in df):
                    raise InvalidAllotypeError(str(df), 'Ensure that the provided dataframe has "ALLOTYPE" and "SEQ" headers.')
                df = df.reset_index().groupby('SEQ').agg(
                        lambda x : '/'.join(set(x)) if isinstance(x, str) else min(x)
                    ).reset_index().sort_values('index')
                seq = '|'.join(df['SEQ'])
                typing = '/'.join(sorted(set(df['ALLOTYPE'])))
        else:
            typing = self.hla
        return str(typing).replace('HLA-', ''), seq

    def _calc_locus(self) -> str:
        """
        Calculates the HLA locus.
        """
        return self.typing.split('*')[0]

    def _calc_hla_class(self) -> str:
        """
        Calculate the HLA class of this particular allele typing.
        """
        if self.locus in ['A', 'C','B']:
            return '1'
        elif 'D' in self.locus:
            return '2'
        else:
            raise InvalidAllotypeError(self.typing,
             "Check this locus. Cannot determine its HLA class")

    def _calc_alleles_hi_res(self) -> List:
        """
        Calculates the hi-res (two-field) alleles of this particular allele typing.
        """
        # if (self.resolution not in ['intermediate', 'low', 'serology'] and not self.is_gl): # or 
        if self.resolution == 'high' and len([field for field in self.fields if field]) <= 2:
            return [self]
        alleles = self.alleles or [self]
        alleles_hi_res = self._redux(self._redux(self.typing, 'W'), 'U2').split('/')
        return [Allotype(allele, expand=False, ref_data=self.ref_data, ard=self.ard) for allele in alleles_hi_res]
    
    def _redux(self, typing : str, mode : str) -> str:
        """
        Calls the reduction functions from py-ard.
        """
        if not self.ard:
            raise Exception('No ARD active. Please address this.')
        if '/' not in typing and self.resolution in ['high', 'allelic']:
            result = self.ard.redux(typing, mode)
        else:
            try:
                result = self.ard.redux_gl(typing, mode)
            except Exception as e:
                if self.resolution == 'g_group':
                    return self.typing
                print(e)
        return result
    
    def _get_two_field_allele(self, allele_name : str) -> str:
        """
        Obtains allele name with only first two fields but keep suffix
        """
        # if self.ard:
        return self.ard.redux(allele_name, 'U2')
        # else:
        #     # print('No ARD')
        #     suffix = ""
        #     if re.search('[A-Z]$', allele_name) and allele_name.count(':') > 2:
        #         suffix = allele_name[-1]
        #     return ':'.join(str(allele_name).split(':')[:2]) + suffix

    def _calc_allele_family(self) -> str:
        """
        Calculates the allele family of this particular allele typing.
        """
        if not self.alleles_hi_res:
            return None
        families = set([str(int(allele.fields[0])) for allele in self.alleles_hi_res])
        return '/'.join(families)

    def calc_and_set_snps(self) -> None:
        results = {}
        for feat_type, seq_features in self.gene_feats.items():
            results[feat_type] = {}
            for feat_name, feat in seq_features.features.items():
                print('a', feat_name, feat, type(feat))
                if feat.snps:
                    results[feat_type][feat_name] = [str(snp) for snp in feat.snps]
        self.snps = results

    def assign_gene_features(self, ref_sequence : str, sequence : str, cds : bool = False):
        """
        Creates gene features object
        """
        if not cds or (cds and not self.sequence and not self.gene_feats):
            self.sequence = Sequence(sequence, ref_sequence=ref_sequence)
            self.gene_feats = GeneFeatures(ref_sequence, sequence, self.locus, cds=cds)

    def add_prot(self, sequence : str, ref_sequence : str) -> None:
        """
        Assigns protein sequence and leader type as well.
        """
        self.protein = Sequence(sequence, seq_type="aa", ref_sequence=ref_sequence)

    def _is_gl(self):
        """ 
        Detects if allele_name is in a GL format
        :return: boolean and list of Allotype objects
        """
        try:
            if '/' in self.typing:
                alleles = []
                prefix = None
                if self.typing.count(':') == 1:
                    prefix = self.typing.split(':')[0]
                pot_alleles = self.typing.split('/')
                self.locus = self.typing.split('*')[0]
                for i, allele in enumerate(pot_alleles):
                    if ':' not in allele and prefix:
                        allele = prefix + ':' + allele
                    if i != 0 and self.locus + '*' not in allele:
                        if ':' not in allele:
                            allele = ':'.join(pot_alleles[0].split(':')[:-1]) + ':' + allele
                        else:
                            allele = self.locus + '*' + allele
                    allele = Allotype(allele, expand=False)
                    alleles.append(allele)
                return True, alleles
            else:
                return False, None
        except Exception as e:
            return False, None

    def get_potential_alleles(self):
        """
        Returns the potential high-resolution alleles.
        If typing is intermediate, expands the alleles via a MAC service.
        :return: list of Allotype objects
        """
        alleles = set(self.alleles_called.split('/')) if self.alleles_called else set()
        if self.alleles and not self.alleles_called:
            return self.alleles
        alleles.update(self._redux(self.typing, 'W').split('/'))
        alleles = sorted(alleles)
        return [Allotype(allele, expand=False, ref_data=self.ref_data) for allele in alleles]

    def _call_mac_api(self) -> List[any]:
        url = "https://hml.nmdp.org/mac/api/expand?typing="
        try:
            response = requests.get(url + 'HLA-' + str(self.typing))
            data = json.loads(response.content)
            expanded_list = [result['expanded'].replace('HLA-','') for result in data]
            return [Allotype(hla_name, expand=False, ref_data=self.ref_data) for hla_name in expanded_list]
        except Exception as e:
            raise InvalidAllotypeError(self.typing, 'Error retrieving expanded alleles: ' + str(e))


    def calc_group(self, group_type : str) -> str:
        """
        Calculates the P group of this allele.
        """
        groups = set()
        if group_type == 'p':
            groups_ref = self.ref_data.p_groups
            if self.p_group:
                return self.p_group
            return self._redux(self.typing, 'lgx') + 'P'
        elif group_type == 'g':
            groups_ref = self.ref_data.g_groups
            if self.g_group:
                return self.g_group
            return self._redux(self.typing, 'G')
        else:
            raise Exception("Ensure the 'group_type' is a 'p' or 'g'. Current value: ", group_type)
        # for group, alleles in groups_ref[self.locus].items():
        #     if Allotype(group, expand=False).locus == self.locus:
        #         for allele in alleles:
        #             if group_type == 'p':
        #                 allele_hi_res = self._get_two_field_allele(allele)
        #                 alleles_hi_res = self.alleles_hi_res if self.alleles_hi_res else [self.typing]
        #                 if allele_hi_res in [str(allele) for allele in alleles_hi_res]:
        #                     groups.add(group)
        #             elif group_type == 'g':
        #                 if allele in [str(allele) for allele in self.alleles or self.get_potential_alleles()]:
        #                     groups.add(group)
        # if len(groups) > 1 and 'XX' not in self.typing:
        #     print(self.typing, 'Multiple {} groups: '.format(group_type), groups)
        #     # raise InvalidAllotypeError(self.typing, 'Multiple P groups were found for this allele.')            
        # elif len(groups) == 1:
        #     return groups.pop()
        # else:
        #     return None
        # return groups
        # return None
    
    def _get_info(self) -> None:
        """
        Extracts the locus name, a list of fields (integers), any variant expression,
        and the resolution level of the allele and assigns it to class variables.
        Raises error if format is invalid.
        """
        if self.is_gl:
            self.resolution = None
            return
        if self.ref_data and self.typing in self.ref_data.ser_antigens:
            m = re.search('([A-Za-z]+)(\d+)', self.typing)
            self.locus, self.fields = m.groups()
            self.resolution = 'serology'
            return
        if ('*' not in self.typing) or (self.typing.count('*') != 1):
            raise InvalidAllotypeError(self.typing, "Typing is neither a known serology value nor contains valid allelic formatting (needs to contain exactly one asterisk '*').")
        locus, fields = self.typing.split('*')
        
        if not re.match('[\dX]+(\:[A-Z\d]+)+[A-Z]?$', fields):
            raise InvalidAllotypeError(self.typing, "The fields are incorrectly formatted. "
                                            "Please include at least two fields (integers). "
                                            "Additional fields are appended with a preceding semicolon ':'. "
                                            "For example, DPB1*07:02 is a valid format but A*07:02: is not.")

        g_group =  bool(re.search('\d+G$', self.typing))
        p_group =  bool(re.search('\d+P$', self.typing))
        var_expression = ((re.search('\d+[A-Z]$', self.typing) and not (g_group or p_group)) and
                          fields[-1] or None)
        if var_expression or g_group or p_group: fields = fields[:-1]

        fields = fields.split(':')

        if self.ref_data and self.typing in self.ref_data.alleles[self.locus].df.index:
            resolution = 'allelic'
        elif g_group:
            resolution = 'g_group'
        elif p_group:
            resolution = 'p_group'
        else:
            resolution = (len(fields) == 1 and 'low' or
                            (g_group or re.match('^[A-Z]+$', fields[1])) and 'intermediate' or
                            re.match('^[0-9]+$', fields[1]) and 'high'  or None)

        while len(fields) < 4:
            fields.append(None)
        if not resolution:
            raise InvalidAllotypeError(self.typing, "The level of typing resolution cannot be determined.")
        if fields[1] == 'XX':
            self.wildcard = True
            resolution = 'low'
        self.fields = fields
        self.resolution = resolution
        self.var_expression = var_expression

    def __str__(self) -> str:
        return self.typing
    
    def __iter__(self):
        yield 'name', self.typing

    def __eq__(self, allotype2 : any) -> bool:
        allotype1 = self
        if allotype1.seq and allotype2.seq:
            return allotype1.seq == allotype2.seq
        return allotype1.typing == allotype2.typing

    def __repr__(self):
        return self.typing

    def serialize(self) -> Dict[str, str]:
        output = {'typing' : self.typing}
        if self.alleles_called:
            output['alleles_called'] = self.alleles_called
        if self.expand:
            output['exact_alleles_called'] = self.exact_alleles_called
            output['resolution'] = self.resolution
        if self.annotation:
            output['annotation'] = self.annotation
        if self.g_group:
            output['g_group'] = self.g_group
        if self.p_group:
            output['p_group'] = self.p_group
        if self.alleles_hi_res and len(self.alleles_hi_res) > 1:
            output['alleles_hi_res'] = [str(allele) for allele in self.alleles_hi_res]
        if self.feats_expanded:
            output['feats_expanded'] = self.feats_expanded
        # if self.alleles and self.expand and len(self.alleles) > 1:
        #     output['alleles'] = [allele.serialize() for allele in self.alleles]
        if self.annotation_seq:
            output['annotation_seq'] = self.annotation_seq
        if self.annotation_diffs:
            output['annotation_diffs'] = self.annotation_diffs
        if self.seq_diffs:
            output['seq_diffs'] = self.seq_diffs
        if self.snps:
            output['snps'] = self.snps
        if self.gene_feat_diffs:
            output['gene_feat_diffs'] = self.gene_feat_diffs
        if self.feats:
            output['feats'] = self.feats.serialize()
        if self.feats_searched:
            output['feats_searched'] = self.feats_searched.serialize()
        return output

class InvalidAllotypeError(Exception):

    def __init__(self, allele_name, message) -> None:
        self.allele_name = allele_name
        self.message = message