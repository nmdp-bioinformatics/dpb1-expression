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
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import Seq
from Bio import Align
from hla_seq_db.snp import Snp
from typing import Tuple, Union, List, Dict

class SeqMatch(Sequence):
    
    def __init__(self,
                 seq_one : Sequence, seq_two : Sequence,
                name : str = None, pos : int = 1, pos_outer : int = None,
                searched : bool = False, index : int = None, location : FeatureLocation = None):
        if seq_one and isinstance(seq_one, str) or seq_one == None:
            seq_one = Sequence(seq_one, name=name)
        if seq_two and isinstance(seq_two, str) or seq_two == None:
            seq_two = Sequence(seq_two, name=name)
        self.seq_one = seq_one
        self.seq_two = seq_two
        self.cds = self.seq_one.cds and self.seq_two.cds if self.seq_one and self.seq_two else None
        self.pos = pos
        self.pos_outer = pos_outer
        self.seq, self.snps, self.distance, self.partial = self._compare_seqs()
        if not name:
            name = self.seq_one.name or self.seq_two.name
        if index == None:
            index = self.seq_one.index or self.seq_two.index
        if not location:
            location = self.seq_one.location or self.seq_two.location
        super().__init__(str(self.seq), name=name, pos=self.pos,
             index=index, searched=searched, location=location)
        
    def _compare_seqs(self) -> Tuple[str, List[Snp], Dict[str, int], bool]:
        # if self.seq_one == self.seq_two:
        #     return '-' * len(self.seq_one.seq), None
#         aligner = Align.PairwiseAligner()
#         score = aligner.score(str(self.seq_one), str(self.seq_two))
        snps = []
        fill_chars = ['.', '*']
        if (self.cds and
            self.seq_one.name != 'codon' and self.seq_two.name != 'codon'):
                alignment = []
                # pos = 1
                for j, (codon_one, codon_two) in enumerate(zip(self.seq_one.codons, self.seq_two.codons)):
                    cod_match = SeqMatch(codon_one, codon_two, 
                                         name='codon', pos=codon_one.pos, pos_outer=codon_one.pos_outer)
                    # pos += len(str(codon_one))
                    alignment.append(cod_match.seq)
                    snps += cod_match.snps
                alignment = ' '.join(alignment)
        else:
            alignment = ""
            for i, (a, b) in enumerate(zip(str(self.seq_one), str(self.seq_two))):
                snp = False
                if a in fill_chars and b in fill_chars:
                    alignment += b
                elif a == '*' or b == '*':
                    alignment == b
                elif a == '.':
                    alignment += a
                    snp = True
                elif b == '.':
                    alignment += b
                    snp = True
                elif a == b:
                    alignment += "-"
                else:
                    snp = True
                    alignment += b
#                     snp.translate_codons()
#                     print(snp)
                if snp:
                    snps.append(Snp(a, b,
                            cds=self.cds,
                            codon_one=self.seq_one if self.seq_one.name == 'codon' else None,
                            codon_two=self.seq_two if self.seq_two.name == 'codon' else None,
                            pos=self.pos + i if self.pos else i, 
                            pos_outer=self.pos_outer, seq=self))
        distances = {'nt' : len(snps)}
        if self.cds:
            distances['aa'] = len(set([snp.pos_outer for snp in snps if snp.missense]))
        partial = '*' in self.seq_one or '*' in self.seq_two
        return alignment, snps, distances, partial

    def get_info(self, print_out : bool = False) -> Union[str, None]:
        results_one = self.seq_one.get_info().split('\n')
        results_two = self.seq_two.get_info().split('\n')
        results_aligned = Sequence(self.seq, name=self.name).get_info(translate=False).split('\n')
        print(results_one)
        print(results_two)
        if self.cds:
            results = [results_one[0], results_one[2], results_one[3], results_aligned[3],
                       results_two[3], results_two[2], results_two[4]]
        else:
            results = [results_one[0], results_one[3], str(self.seq), results_two[3]]
        results = '\n'.join(results)
        if print_out:
            print(results)
        else:
            return results
    
    def __repr__(self):
        return self.seq_two.__repr__()
        # return '\n'.join([str(self.seq_one),
        #                   self.seq,
        #                   str(self.seq_two)])

    def __add__(self, seq_match : any) -> any:
        self.seq_two += seq_match.seq_two
        return self