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
from pyard.ard import ARD
import pandas as pd
from .dataset import Dataset

class Haplotype(object):

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
        self.delimiter = '^'
        self.glstring, self.allotypes = self._process_hla()
        self.annotation = {}

    def _process_hla(self) -> Tuple[str, List[Allotype]]:
        if isinstance(self.hla, list):
            allotypes = self.hla
        elif isinstance(self.hla, str):
            allotypes = self.hla.split(self.delimiter)
        else:
            raise Exception("Please input a list or string haplotype input.")
        allotypes = [Allotype(allotype, ref_data=self.ref_data, ard=self.ard, verbose=self.verbose)
                                        for allotype in allotypes]
        glstring = self.delimiter.join(sorted([allotype.typing for allotype in allotypes]))
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
    
    def __str__(self) -> str:
        return self.glstring

    def __repr__(self) -> str:
        return self.glstring

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
