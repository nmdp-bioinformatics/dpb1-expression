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
from .sequence import Sequence

class GeneFeatures(object):
    def __init__(self, ref_sequence : str, delimited_sequence : str, locus : str, cds : bool):
        self.utr5 = None
        self.exons = []
        self.introns = []
        self.utr3 = None
        self.locus = locus
        self.cds = cds
        self.seq = Sequence(delimited_sequence, ref_sequence=ref_sequence)
        self._parse_sequence(ref_sequence, delimited_sequence)
    
    def _parse_sequence(self, ref_sequence : str, delimited_sequence : str) -> None:
        """
        Splits IMGT sequence delimited by '|' and assigns gene features based on
        this order: 5' UTR, exon 1, intron 1, exon 2, intron 2, exon 3, intron 4,
        exon 5, intron 5, exon 6, intron 6, exon 7, 3' UTR.
        """
        gene_features = delimited_sequence.split('|')
        ref_gene_features = ref_sequence.split('|')
        gene_features.reverse()
        ref_gene_features.reverse()
        if self.cds:
            n = len(gene_features)
        else:
            n = int((len(gene_features) - 3) / 2) + 1
            self.utr5 = self._create_seq(gene_features, ref_gene_features)
        for i in range(1, n):
            self.exons.append(self._create_seq(gene_features, ref_gene_features))
            if not self.cds:
                self.introns.append(self._create_seq(gene_features, ref_gene_features))
        self.exons.append(self._create_seq(gene_features, ref_gene_features)) # exon 7
        if not self.cds:
            self.utr3 = self._create_seq(gene_features, ref_gene_features)
    
    def _create_seq(self, gene_features : list, ref_gene_features : list) -> Sequence:
        """
        Creates Sequence object from a list of gene features, 
        popping the last item in the list,
        in conjunction with a reference list.
        """
        return Sequence(gene_features.pop(), ref_sequence=ref_gene_features.pop())
