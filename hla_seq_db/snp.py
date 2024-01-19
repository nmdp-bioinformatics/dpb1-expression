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
from typing import Tuple

class Snp(object):
    def __init__(self, nt_one : str, nt_two : str, cds : bool = None,
                 pos : int = None, pos_outer : int = None,
                 seq : Sequence = None,
                 codon_one : Sequence = None, codon_two : Sequence = None):
        self.nt_one = nt_one
        self.nt_two = nt_two
        self.type = self._determine_type(nt_one, nt_two)
        self.cds = cds
        self.pos = pos
        self.pos_outer = pos_outer
        self.codon_one = codon_one
        self.seq = seq
        self.codon_two = codon_two
        self.aa_one, self.aa_two, self.missense = self._calc_AAs()

    def _determine_type(self, nt_one : str, nt_two : str) -> str:
        if nt_one != '.' and nt_two != '.':
            return 'sub'
        elif nt_one == '.':
            return 'ins'
        elif nt_two == '.':
            return 'del'            

    
    def _calc_AAs(self) -> Tuple[str, str, bool]:
        if self.codon_one and self.codon_two:
            aa_one = self.codon_one.translate()
            aa_two = self.codon_two.translate()
            missense = aa_one != aa_two
            return aa_one, aa_two, missense
        return None, None, False
    
    def serialize(self) -> tuple:
        output = {'nt_one' : self.nt_one,
                  'nt_two' : self.nt_two,
                  'pos' : self.pos}
        if self.codon_one:
            output['codon_one'] = self.codon_one
        if self.codon_two:
            output['codon_two'] = self.codon_two
        if self.aa_one:
            output['aa_one'] = self.aa_one
        if self.aa_two:
            output['aa_two'] = self.aa_two
#         if self.pos_codon:
#             output['pos_codon'] = self.pos_codon
        return output
        
    def __repr__(self):
        if self.type == 'sub':
            print_str = "c.{pos}{nt_one}>{nt_two}".format(
                pos=self.pos, nt_one=self.nt_one, nt_two=self.nt_two)
        elif self.type == 'ins':
            print_str = "c.{pos}ins{nt_two}".format(pos=self.pos, nt_two=self.nt_two)
        elif self.type == 'del':
            print_str = "c.{pos}del{nt_one}".format(pos=self.pos, nt_one=self.nt_one)
        if self.cds:
            coding_str = " (p.{aa_one}{aa_pos}{aa_two})".format(
                aa_one=self.codon_one.translate(3) if self.codon_one else None,
                aa_two=self.codon_two.translate(3) if self.codon_two else None,
                aa_pos=self.pos_outer)
            print_str += coding_str
        return print_str

# class SNP(object):
#     def __init__(self, nt_pat : str, 
#                 nt_don : str, pos : int, gene_feat : str,
#                 codon_pat : str = None, codon_don : str = None,
#                 pos_codon : int = None):
#         self.nt_pat = nt_pat
#         self.nt_don = nt_don
#         self.pos = pos
#         self.ambiguous = False
#         self.gene_feat = gene_feat
#         self.codon_pat = codon_pat
#         self.codon_don = codon_don
#         self.aa_pat = self.aa_don = None
#         self.pos_codon = None
#         if ('X' == nt_pat) or ('X' == nt_don):
#             self.ambiguous = True
    
#     def calc_amino_acids(self) -> None:
#         if self.codon_pat and self.codon_don:
#             try:
#                 self.aa_pat = str(Seq(self.codon_pat, IUPAC.unambiguous_dna).transcribe().translate())
#                 self.aa_don = str(Seq(self.codon_don, IUPAC.unambiguous_dna).transcribe().translate())
#             except Exception as e:
#                 print(e)
        
#     def serialize(self) -> tuple:
#         output = {'nt_pat' : self.nt_pat,
#                   'nt_don' : self.nt_don,
#                   'pos' : self.pos}
#         if self.codon_pat:
#             output['codon_pat'] = self.codon_pat
#         if self.codon_don:
#             output['codon_don'] = self.codon_don
#         if self.aa_pat:
#             output['aa_pat'] = self.aa_pat
#         if self.aa_don:
#             output['aa_don'] = self.aa_don
#         if self.pos_codon:
#             output['pos_codon'] = self.pos_codon
#         return output
    
#     def __repr__(self):
#         return "{}>{} ({})".format(self.nt_pat,
#                      self.nt_don, self.pos)