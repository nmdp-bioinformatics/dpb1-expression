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
# from .sequence_search import SeqSearch
from .sequence_features import SeqFeatures
from .annotation import Annotation
# from hla_seq_db.allotype import Allotype
import pandas as pd
from typing import Tuple, Dict, List, Union
import re
import inflect

class Features(object):
    def __init__(self, dfs : Dict[str, pd.core.frame.DataFrame], wildcard : bool = False,
            field_level : int = 2):
        self.wildcard = wildcard
        self.field_level = field_level
        self.seqs, self.seq_anns, self.anns = self._parse_data(dfs)

    def _parse_data(self, dfs : Dict[str, pd.core.frame.DataFrame]) -> Tuple[SeqFeatures, SeqFeatures, Annotation]:
        seqs = anns = seq_anns =None
        for feat_type, df in dfs.items():
            if 'gene_feats' == feat_type:
                seqs = self.generate_cons_gene_feats(df)
            elif 'gene_ann_feats' == feat_type:
                seq_anns = self.generate_cons_gene_feats(df)
            elif 'ann_feats' == feat_type:
                anns = self.calc_common_annotation(df)
        return seqs, seq_anns, anns

    def set_searched_features(self, searched_feats : any):
        self.gene_feats = searched_feats

    def generate_cons_gene_feats(self, 
            df : pd.core.frame.DataFrame) -> SeqFeatures:
        """
        Iterates through the gene feature columns of a dataframe to
        create the consensus sequence for each gene feature. Ignores
        empty typing '*' if other alleles have typing. If there is more than
        one nucleotide at a given position, then an "X" is inserted.
        """
        return SeqFeatures(df, wildcard=self.wildcard)

    def calc_common_annotation(self, df : pd.core.frame.DataFrame,
            unk_val : str = None) -> Tuple[Dict[str, str], Union[List[any], None]]:
        return Annotation(df, unk_val=unk_val)

    # def get_cons_seq(self, typing : str = None, gene_feat : str = None,
    #         df : pd.core.frame.DataFrame = None,
    #         alleles : List[str] = None, spacers : bool = True,
    #         most_freq_nt : bool = False, cds_only : bool = False) -> str:
    #     """
    #     Get consensus sequence of the provided HLA typing.
    #     Any ambiguity is denoted with an 'X'. You may also
    #     give a locus as a parameter i.e., 'B' to obtain 
    #     the consensus sequence for the whole locus.
    #     """
    #     # TODO: Refactor this
    #     whole_locus = False
    #     if isinstance(df, pd.core.frame.DataFrame) and not df.empty:
    #         df = df
    #     elif typing in self.loci and not alleles:
    #         df = self.alleles[typing]
    #     else:
    #         if typing in self.loci:
    #             locus = typing
    #             if locus in self.cons_seq_whole_locus:
    #                 return self.cons_seq_whole_locus[locus]
    #             whole_locus = True
    #         elif typing:
    #             # if not isinstance(typing, Allotype):
    #             #     allele = Allotype(typing, ref_data=self)
    #             # else:
    #             allele = typing
    #             locus = allele.locus
    #         if not alleles:
    #             alleles = [str(allele) for allele in allele.get_potential_alleles()]
    #         df = self.alleles[locus]
    #         df = df[df.index.isin(alleles)]
    #     seqs_feat_out = []
    #     for col in df.columns:
    #         if not ('utr' in col or 'exon' in col or 'intron' in col or 'str' in col):
    #             continue
    #         if cds_only and 'exon' not in col:
    #             continue
    #         if gene_feat and gene_feat != col:
    #             continue
    #         seqs_feat = df[col]
    #         seq_feat_out = ""
    #         seqs_feat = seqs_feat.dropna()
    #         if seqs_feat.empty:
    #             break
    #         if col == 'expression_experimental':
    #             continue #TODO: Remove this
    #         for i in range(len(seqs_feat.iloc[0])):
    #             nts = []
    #             nt = None
    #             for seq in seqs_feat:
    #                 if isinstance(seq, str):
    #                     if not most_freq_nt:
    #                         if seq[i] != '*':
    #                             nts.append(seq[i])
    #                         # nts.add(seq[i])
    #                         # if len(set(nts)) >= 2:
    #                         #     nt = 'X'
    #                         #     break
    #                     else:
    #                         if seq[i] != '*':
    #                             nts.append(seq[i])
    #             if not nt:
    #                 if not len(nts):
    #                     nt = '*'
    #                 else:
    #                     if most_freq_nt:
    #                         nt = max(set(nts), key = nts.count)
    #                     else:
    #                         nt_counts = pd.Series(nts).value_counts(normalize=True).to_dict()
    #                         # if col == 'exon_2':
    #                         #     print(nt_counts)
    #                         if (len(nt_counts) == 1) or list(nt_counts.values())[0] >= .99:
    #                             nt = list(nt_counts.keys())[0]
    #                         else:
    #                             nt = 'X'
    #                         # if len(set(nts)) == 1:
    #                         #     nt = nts[0]
    #                         # else:
    #                         #     nt = 'X'
    #             seq_feat_out += nt
    #         seqs_feat_out.append(seq_feat_out)
    #     seq_joined = '|'.join(seqs_feat_out)
    #     if not spacers:
    #         seq_joined = seq_joined.replace('.', '')
    #     if whole_locus:
    #         self.cons_seq_whole_locus[locus] = seq_joined
    #     return seq_joined

    def __repr__(self):
        return str(self.__dict__)

    def serialize(self):
        # return {feat_name : str(feat) 
        results = {'seqs' : self.seqs,
                   'seq_anns' : self.seq_anns,
                   'anns' : self.anns}
        for label, ann in results.items():
            if ann:
                results[label] = ann.serialize()
        return results
        # return {'gene_feats' : self.gene_feats.serialize(),
        #     'gene_feat_annotations' : self.gene_feat_annotations.serialize(),
        #      'annotations' : self.annotations}
        #             for feat_name, feat in self.gene_feats.items()}