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
from Bio.Seq import Seq
from Bio import Align
from hla_seq_db.snp import Snp
from typing import Tuple, Union, List, Dict
from Bio.SeqFeature import SeqFeature, FeatureLocation
from hla_seq_db.sequence_feat_searched import SeqFeatSearched
import pandas as pd
import re

class SeqFeatures(object):
    
    def __init__(self,
                 features : Union[Dict[str, Dict[str, str]],
                             pd.core.frame.DataFrame, Dict[str, SeqFeatSearched]],
                 wildcard : bool = False):
        self.wildcard = wildcard
        self.introns = []
        self.exons = []
        self.utrs = []
        self.feats = self._parse_features(features)

    def _parse_features(self, features : Union[Dict[str, Dict[str, str]], pd.core.frame.DataFrame]) -> Dict[str, Union[Sequence, SeqFeatSearched]]:
        if isinstance(features, pd.core.frame.DataFrame):
            features = features.to_dict()
        features_parsed = {}
        for i, (feature, seq) in enumerate(features.items()):
            if 'allele' not in feature:
                if isinstance(seq, str) or isinstance(seq, dict):
                    seq = Sequence(seq, name=feature, index=i, 
                    wildcard=self.wildcard if feature != 'exon_2' else False,
                    protein=feature == 'P2')
                if 'exon' in feature:
                    self.exons.append(seq)
                elif 'intron' in feature:
                    self.introns.append(seq)
                elif 'utr' in feature:
                    self.utrs.append(seq)
                features_parsed[feature] = seq
        return features_parsed

    # def ordered_features(self) -> Dict[str, Sequence]:
    #     order = ['exon_2', 'intron_2', 'exon_1', 'intron_1', 'exon_3', 'exon_4', 'exon_5', 'exon_6']
    #     for feat, seq in self.feats.items():
            

    def get_info(self, print_out : bool = False) -> Union[str, None]:
        display = None
        for j, (feature, seq) in enumerate(self.feats.items()):
            rows = seq.get_info().split('\n')
            print(feature, seq, rows)
            if not display:
                display = [[] for i in range(len(rows))]
            for i in range(len(rows)):
                print(i, rows[i])
                display[i].append(rows[i])
        display = '\n'.join(['|'.join(row) for row in display])
        if print_out:
            print(display)
        else:
            return display
            

    def __repr__(self):
        return str(self.feats)

    def serialize(self):
        return {feat_name : str(seq) for feat_name, seq in self.feats.items()}
        