#
# Copyright (c) 2024 NMDP.
#
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
from .allotype import Allotype
from .ref_data import RefData
from typing import Union, Dict, List, Tuple
from pyard import ARD
import pandas as pd
from .dataset import Dataset

class Genotype(object):

    def __init__(self, 
             hla : Union[str, List[Allotype], pd.core.frame.DataFrame],
             id=None, verbose : bool = False,
             ref_data : RefData = None, ard : ARD = None, sire : str = None) -> None:
        """
        Represents a pair of HLA-DPB1 alleles. Needs to be delimited by either a
        '+' (unphased) or '~' (phased) character. The order of the alleles is stored
        based on a numerical sort.

        :param typing:  Two allele names delimited by either '+' or '~' or a
            Dataframe with an 'allotype' and 'seq' column and an optional 'allotype' column.
        """
        self.verbose = verbose
        self.hla = hla
        self.id = id
        self.sire = sire
        self.ref_data = ref_data
        self.ard = ard
        self.delimiter = '+'
        self.glstring, self.allotypes = self._process_hla()
        self.annotation = {}

        if len(self.allotypes) == 1:
            self.allotypes += self.allotypes
        if len(self.allotypes) != 2:
            raise InvalidGenotypeError(self.allotypes,
                    "Genotype does not contain exactly two alleles")
        self.homozygous = self.allotypes[0].calc_group('p') == self.allotypes[1].calc_group('p')
        self.min_tce = None

        self._set_name_and_alleles(self.allotypes)
        self.locus = self._validate_locus()

    def _process_hla(self) -> Tuple[str, List[Allotype]]:
        if isinstance(self.hla, pd.core.frame.DataFrame):
            if 'ALLOTYPE' in self.hla:
                df = self.hla.copy()[['ALLOTYPE', 'PHASE', 'SEQ']].drop_duplicates()
                df.loc[:, 'LOCUS'] = df['ALLOTYPE'].str.split('*').str[0]
                allotypes = dict(tuple(df.groupby('PHASE'))).values()
                # allotypes = [Allotype(df_allotype, ref_data=self.ref_data, ard=self.ard, verbose=self.verbose)
                #         for df_allotype in df_allotypes]
            else:
                df = self.hla.copy()
                ds = Dataset(df)
                allotypes = [ds.get_dataframe(headers_num=header_num) for header_num in ds.headers_num]
        else:
            typing = self.hla
            if isinstance(typing, list):
                allotypes = typing
            else:
                if '+' in typing:
                    self.delimiter = '+'
                elif '~' in typing:
                    self.delimiter = '~'
                else:
                    raise InvalidGenotypeError(typing,
                            "Genotype is not split by '+' or '~'")
                allotypes = typing.split(self.delimiter)
        allotypes = [Allotype(allotype, ref_data=self.ref_data, ard=self.ard, verbose=self.verbose)
                                        for allotype in allotypes]
        glstring = '+'.join(sorted([allotype.typing for allotype in allotypes]))
        return glstring, allotypes


    def _validate_locus(self):
        """
        Validates of alleles have the same loci.
        """
        loci = set([allele.locus for allele in self.allotypes])
        if len(loci) == 1:
            return loci.pop()
        raise InvalidGenotypeError(self.glstring,
             'There is an appropriate amount of loci. Please include only one HLA locus.')
    
    def flip(self):
        alleles = self.allotypes
        alleles.reverse()
        # self.flip_matched = not self.flip_matched
        self.flip_matched = True
        self._set_name_and_alleles(alleles)
    
    def get_annotations(self):
        if not self.annotation:
            for allele in self.allotypes:
                if allele.feats:
                    ann = allele.feats.anns.serialize()
                    if ann:
                        for feat, value in ann.items():
                            if feat not in self.annotation:
                                self.annotation[feat] = []
                            self.annotation[feat].append(value)

    def _set_name_and_alleles(self, alleles):
        self.first_allele, self.second_allele = self.allotypes = alleles
        self.name = str(self.first_allele) + self.delimiter + str(self.second_allele)

    def _sort(self, alleles):
        """
        Sorts the two alleles based on a numeric sorting of the fields, with
        priority on the first (left-most) field.
        Returns the sorted list of alleles.

        Empty fields have lower priority over populated, numeric fields.

        :param alleles: 
        :type alleles: List of Allotype objects
        """
        if not alleles[0].fields or not alleles[1].fields:
            return alleles
        fields_a, fields_b = alleles[0].fields, alleles[1].fields
        i = -1
        reversed = False
        while i < len(fields_a) - 1 and i < len(fields_b) - 1 and not reversed:
            i += 1
            if re.match('\d+', str(fields_a[i])) and re.match('\d+', str(fields_b[i])):
                field_a, field_b = int(fields_a[i]), int(fields_b[i])
            else:
                field_a, field_b = str(fields_a[i]), str(fields_b[i])
            if field_a > field_b:
                reversed = True
            elif field_a < field_b:
                break
        if reversed or (not fields_a[i] and fields_b[i]):
            alleles.reverse()
            # self.flip_sorted = True
        return alleles

    def first(self):
        return self.allotypes[0]

    def second(self):
        return self.allotypes[1]

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name
    
    def __iter__(self):
        yield 'id', self.id
        yield 'name', self.name
        # yield 'flip_sorted', self.flip_sorted
        # yield 'flip_matched', self.flip_matched
        yield 'allele_one', dict(self.first_allele)
        yield 'allele_two', dict(self.second_allele)

    def serialize(self) -> Dict[str, str]:
        serialization = {'genotype' : self.name,
                'annotation' : self.annotation,
                'allotype_one' : self.allotypes[0].serialize(),
                'allotype_two' : self.allotypes[1].serialize()}
        if self.id:
            serialization['id'] = self.id
        # if self.sire:
        #     serialization['sire'] = self.sire
        return serialization

class InvalidGenotypeError(Exception):

    def __init__(self, typing, message) -> None:
        self.typing = typing
        self.message = message