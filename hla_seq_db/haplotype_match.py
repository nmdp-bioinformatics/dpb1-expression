#
# Copyright (c) 2024 NMDP.
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
from .allotype import Allotype
# from .genotype import Genotype
from .haplotype import Haplotype
from .ref_data import RefData
from .allotype_match import AllotypeMatch
from pyard.ard import ARD
from typing import Tuple, List, Dict

class HaplotypeMatch(object):

    def __init__(self, haplotype_patient: Haplotype, haplotype_donor: Haplotype, ref_data : RefData = None,
                id : str = None,
                ard : ARD = None) -> None:
        self.ref_data = ref_data
        self.ard = ard
        self.id = id
        self.haplotype_patient = haplotype_patient
        self.haplotype_donor = haplotype_donor
        self._process_hla()
        self.matches, self.score, self.matched = self._get_matches()

    def _process_hla(self):
        if isinstance(self.haplotype_patient, str) or isinstance(self.haplotype_patient, list):
            self.haplotype_patient = Haplotype(self.haplotype_patient, ard=self.ard, ref_data=self.ref_data)
        if isinstance(self.haplotype_donor, str) or isinstance(self.haplotype_donor, list):
            self.haplotype_donor = Haplotype(self.haplotype_donor, ard=self.ard, ref_data=self.ref_data)

    def _get_matches(self) -> Tuple[List[AllotypeMatch], int, bool]:
        """
        Compares the allotype matches between the forward and reverse versions of the haplotypes
        to get the best match grade combination, which is whatever contains the highest match grade.
        """
        allotypes_pat = {allotype_pat.locus : allotype_pat for allotype_pat in self.haplotype_patient.allotypes}
        allotypes_don = {allotype_don.locus : allotype_don for allotype_don in self.haplotype_donor.allotypes}
        
        matches = []
        score = 0
        matched = True
        for locus, allotype_pat in allotypes_pat.items():
            if locus in allotypes_don:
                allotype_don = allotypes_don[locus]
                match = AllotypeMatch(allotype_pat, allotype_don, ref_data=self.ref_data)
                matches.append(match)
                score += match.score
                if not match.matched:
                    matched = False
        return matches, score, matched

    def __repr__(self) -> str:
        return '{} - {} ({})'.format(self.haplotype_patient,
                        self.haplotype_donor,
                        'matched' if self.matched else 'mismatched')