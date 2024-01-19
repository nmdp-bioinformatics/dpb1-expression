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
from __future__ import annotations
from hla_seq_db.sequence_match import SeqMatch
from Bio.SeqFeature import FeatureLocation
import hla_seq_db
from hla_seq_db.sequence import Sequence
from typing import List, Union, Tuple

class SeqFeatSearched(SeqMatch):
    
    def __init__(self,
                 seq_ref : Sequence, seq : Sequence,
                 location : FeatureLocation = None, searched : bool = None,
                 name : str = None,
                 imputed : bool = False, index : int = None):
        self.seq_ref = seq_ref
        self.seq = seq
        self.imputed = imputed
        super().__init__(seq_ref, seq, index=index, location=location, name=name, searched=searched)
        self.alleles, self.novel = self._determine_matched_alleles(seq_ref, seq)
    
    def _determine_matched_alleles(self, seq_ref : Sequence, seq : Sequence) -> Tuple[Union[None, List[str]], bool]:
        if seq_ref and seq in seq_ref:
            if not isinstance(seq, str):
                seq_raw = seq.raw_seq()
            else:
                seq_raw = seq
            return [seq_allele.allotype 
                    for seq_allele in seq_ref.seqs
                     if (seq_allele == seq_raw) and not isinstance(seq_allele, str)], False
        return None, True

    def extract(self, seq_feat : hla_seq_db.sequence_features.SeqFeatures, name : str) -> SeqFeatSearched:
        return SeqFeatSearched(Sequence(seq_feat.extract(self.seq_one.raw_seq()), name=name),
                               Sequence(seq_feat.extract(self.seq_two.raw_seq()), name=name))