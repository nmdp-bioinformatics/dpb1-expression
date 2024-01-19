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
from Bio import motifs, Alphabet
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped
from typing import Union, List, Tuple, Dict
import pandas as pd
from Bio.SeqUtils import seq3
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re

class Sequence(Seq):
    def __init__(self, seqs : Union[list, str],
                       name : str = None, allotype : str = None,
                       pos : int = 1, pos_outer : int = 1,
                       protein : bool = False,
                       searched : bool = False, wildcard : bool = False,
                       index : int = None, location : FeatureLocation = None):
        self.name = name # Change to label
        self.wildcard = wildcard
        self.protein = protein
        self.gene_feat, self.cds, self.ard = self._parse_name(self.name)
        self.allotype = allotype
        self.cds, self.frame = self._determine_frame()
        self.motif = None
        self.seqs, self.seq, self.empty = self._parse_seqs(seqs)
        self.codons = None
        self.pos = pos
        self.pos_outer = pos_outer
        self.protein = None
        self.index = index
        self.searched = searched
        self.location = location
        if not self.empty and self.cds and self.name != 'codon' or (self.seq and ' ' in self.seq):
            self.seq, self.codons = self.parse_into_codons()
        super().__init__(str(self.seq))
    
    def _parse_name(self, name : str) -> Union[Tuple[bool, bool, bool], None]:
        if not name:
            return None, None, None
        gene_feat = bool(re.match('^(((exon)|(intron))_\d|utr[35])$', name))
        cds = gene_feat and 'exon' in name
        ard = name == 'exon_2'
        return gene_feat, cds, ard
    
    def _determine_frame(self) -> int:
        """
        Determines how many nucleotides to skip based on the
        name i.e., gene feature, provided, e.g., 'exon_1'
        """
        #DPB1
        if self.name:
            mapper = {'exon_1' : 0} #,
                    #   'exon_2' : 2,
                    #   'exon_3' : 2,
                    #   'exon_4' : 2}
            if self.name in mapper:
                return True, mapper[self.name]
            if self.cds:
                return True, 2
            if 'codon' in self.name:
                return True, 3
        return False, 0
        
    def _parse_seqs(self, seqs : Union[list, str, 
            Dict, pd.core.series.Series]) -> Tuple[List[Union[str, any]], str, bool]:
        seq = None
        if isinstance(seqs, Seq):
            seqs = str(seqs)
        if isinstance(seqs, str):
            if re.match('^[.*]*$', seqs):
                return [], seqs, True
            seq = seqs
            seqs = [seqs]
        elif isinstance(seqs, pd.core.series.Series) or isinstance(seqs, dict):
            seqs = [Sequence(seq, allotype=allotype)
                    for allotype, seq in seqs.items()]
            seqs = [seq for seq in seqs if seq.seq]
        elif isinstance(seqs, list) and len(seqs) and isinstance(seqs[0], str):
            seqs = [Sequence(seq) for seq in seqs]
        if not seq and not seqs:
            return [], '', True
        if not isinstance(seqs, list):
            seqs = list(seqs)
        if not seq and seqs:
            seq = self.get_consensus(seqs)
        return seqs, seq, False

    def parse_into_codons(self) -> Tuple[str, List[str]]:
        seq = None
        codon = ''
        codons = []
        j = 0
        if not seq or not isinstance(seq, str):
            seq = str(self.seq)
        if ' ' in seq:
            seqs = seq.split(' ')
            return seq, [Sequence(c, name='codon', 
                            pos=sum([len(seq) for seq in seqs[:i]]) + 1,
                            pos_outer=i+1)
                            for i, c in enumerate(seqs)]
        codon_pos = 1
        n_spacers = 0
        for i in range(len(seq)):
            nt = seq[i]
            # print(nt, j, codons, codon)
            codon += nt
            if nt != '.':
                j += 1
            else:
                n_spacers += 1
            if j and (j % 3 == self.frame) and codon or (i == (len(seq) - 1)):
                if nt != '.' or (i == len(seq) - 1):
                    # codons.append(codon)
                    codons.append(Sequence(codon, 
                            name='codon', 
                            pos=codon_pos, pos_outer=(len(codons)+1)))
                    codon_pos = j + 1 + n_spacers
                    codon = ''
            # spaced_seq += 
        # codons.append(codon)
        # codons.append(Sequence(codon, name='codon', pos=j))
        # self.codons = codons
        # print(codons[0].__dict__)
        # return ' '.join(codons), codons
        return ' '.join([str(c) for c in codons]), codons

    def get_consensus(self, seqs : List[Sequence]) -> Union[str, None]:
        len_init = len(seqs[0])
        seqs_empty = seqs
        seqs = [seq for seq in seqs if not seq.empty]
        if seqs:
            if self.wildcard:
                df = pd.DataFrame(seqs)
                consensus = ""
                for i in range(len(seqs[0])):
                    nts = set(df[0].str[i])
                    if '*' in nts:
                        nts.remove('*')
                    if len(nts) == 1:
                        consensus += nts.pop()
                    # elif len(nts) == 0:
                    #     consensus += 'N'
                    elif '.' in nts:
                    #     nts.remove('.')
                    #     consensus += '[{}]?'.format(''.join(nts))
                        consensus += '.'
                    else:
                        consensus += 'N' #'[{}]'.format(''.join(nts))
                return consensus
            else:
                seqs = [str(seq).replace('*', 'N') for seq in seqs]
                m = self.get_motif(seqs)
                consensus = str(m.degenerate_consensus)
                if not len_init == len(consensus):
                    consensus = ""
                    for i in range(len_init):
                        stats = m.counts[:, i]
        #                 if '-' in stats:
        #                     print(stats)
                        nt = max(stats, key=stats.get)
                        consensus += nt
        #             consensus = Seq(consensus)
                consensus = consensus.replace('-', '.')\
                                    .replace('N', '*')
                return consensus
        else:
            return str(seqs_empty[0])
    
    def get_motif(self, seqs : List[Sequence] = None) -> motifs.Motif:
        if not seqs:
            seqs = self.seqs
        if seqs and not isinstance(seqs[0], Sequence):
            seqs = [Sequence(seq) for seq in seqs]
        seqs = [seq for seq in seqs if not seq.empty]
        alph = Gapped(IUPAC.extended_protein 
                if self.protein else IUPAC.ambiguous_dna, '-')
        m = motifs.create([Seq(str(seq).replace('.', '-').replace('*', '-'), alph) 
                       for seq in seqs], 
                      alphabet=alph)
        self.motif = m
        return m
    
    def expand(self) -> str:
        # SeqUtils.nt_search(seq, seq)
        ambig_mapper = IUPAC.IUPACData.ambiguous_dna_values
        ambig_mapper = {k : v for k, v in ambig_mapper.items() if k != v}
        seq = self.seq
        seq = seq.replace('N', '.')
        for code, values in ambig_mapper.items():
            seq = seq.replace(code, '[{}]'.format(values))
        return seq

    def get_info(self, print_out : bool = False,
                 translate : bool = True) -> Union[str, None]:
        if self.codons:
            results = []
            for i, codon in enumerate(self.codons):
                results.append(codon.get_info(translate=translate).split('\n')[1:])
            result_len = len(results[0])
            results = [' '.join([result[i] 
                        for result in results]) 
                        for i in range(result_len)]
            results = [self.name.ljust(len(results[0]))] + results
            results = '\n'.join(results)
        else:
            els = [self.name,
                   self.pos_outer if self.cds else '', 
                   self.translate(3) if self.cds and translate else '', 
                   self.seq,
                   self.pos if self.cds else '']
            els = [str(el) for el in els if el != None]
            max_len = max([len(el) for el in els])
            results = [el.ljust(max_len) for el in els]
            results = '\n'.join(results)
        if print_out:
            print(results)
        else:
            return results

    def export_weblogo(self, filename : str) -> None:
        if not self.motif:
            print('Motif not generated already')
            self.get_motif()
        self.motif.weblogo(filename)
        
    def raw_seq(self, fillers : bool = True) -> str:
        seq = self.seq.replace(' ', '')
        if not fillers:
            seq = seq.replace('*', '').replace('.', '')
        return seq
        
    def get(self, pos : int, mode : str = 'nt') -> Union[str, Seq]:
        if mode == 'codon':
            if not self.codons:
                raise Exception('No codons available.')
            return self.codons[pos]
        elif mode == 'nt':
            return Sequence(self.raw_seq()[pos])
        elif mode == 'aa':
            return self.codons[pos].translate()
        
    def translate(self, n_letters : int = 1) -> str:
        if self.codons:
            return ''.join([str(codon.translate(n_letters=n_letters))
                     for codon in self.codons])
        if self.protein and len(self.protein) == n_letters:
            return self.protein
        seq = Seq(self.seq.replace('.', '')).translate()
        if n_letters == 3:
            seq = seq3(seq)
        self.protein = seq
        return seq
        # try:
        #     is_empty = all([nt == '*' for nt in sequence])
        #     if len(sequence) == 3 and sequence[2] == '*' and not is_empty:
        #         possible_AAs = set([self.translate(sequence[:2] + nt) for nt in ['A', 'T', 'G', 'C']])
        #         if len(possible_AAs) == 1:
        #             return possible_AAs.pop()
        #         else:
        #             return '*'
        #     if len(sequence) <= 3:
        #         if is_empty:
        #             return '*'
        #         return str(Seq(sequence) \
        #                             .transcribe().translate())
        #     else:
        #         return self.translate(sequence[:3]) + self.translate(sequence[3:])
        # except Exception as e:
        #     print(e)
        #     return '*'

    def __repr__(self) -> str:
        name = self.name.capitalize().replace('_', '') if self.name else 'Seq'
        seq = (self.seq[:54] + '...' + self.seq[-3:]) if len(self.seq) > 60 else self.seq
        seq_str = "{}('{}')".format(name, seq)
        return seq_str
        
    def __contains__(self, seq2 : any) -> bool:
        seq1 = self.seq
        seq2 = str(seq2)
        return (seq2 in seq1) or (seq2 in self.seqs)

    def __add__(self, seq2 : any) -> Sequence:
        if (self.name != seq2.name):
            raise Exception('Cannot add sequences with difference names:', self.name, seq2.name)
        if (len(self.seq) != len(seq2.seq)):
            raise Exception('Cannot add sequences with different lengths:', len(self.seq), len(seq2.seq), self.name, self.allotype)
        seq = ""
        for i in range(len(self.seq)):
            nt1 = self.seq[i]
            nt2 = seq2.seq[i]
            if (nt1 != '*') and (nt2 != '*'):
                if (nt1 == nt2):
                    seq += nt1
                else:
                    raise Exception('These sequences overlap at position {} with different nucleotides'.format(i), nt1, nt2)
            if (nt1 != '*'):
                seq += nt1
            else:
                seq += nt2
        self.seq = seq
        return self