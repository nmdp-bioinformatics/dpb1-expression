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
from .genotype_match import GenotypeMatch
from .genotype import Genotype
from .ref_data import RefData
from typing import List, Union, Dict, Tuple
import pandas as pd
from pyard import ARD

class MatchList(object):
    
    def __init__(self, geno_recip : Union[Genotype, str], 
                 geno_donors : List[Union[Genotype, str, Dict[str, str]]],
                 matches : List[GenotypeMatch] = None,
                 ref_data : RefData = None, verbose : bool = False,
                 ard : ARD = None):
        """
        Represents a list of GenotypeMatch objects between a recipient
        and potential donor stem cell sources.
        """
        self.ref_data = ref_data
        self.ard = ard
        self.verbose = verbose
        self.df_sorter = None
        self.geno_recip, self.geno_donors = self._process_genotypes(geno_recip, geno_donors)
        self.matches = self._assign_matches(matches)

    def _process_genotypes(self, geno_recip : Union[Genotype, str], 
                 geno_donors : List[Union[Genotype, str, Dict[str, str]]]) -> Tuple[Genotype, List[Genotype]]:
        if isinstance(geno_recip, str):
            geno_recip = Genotype(geno_recip, ref_data=self.ref_data, ard=self.ard)
        if geno_donors and isinstance(geno_donors, list):
            if isinstance(geno_donors[0], str):
                geno_donors = [Genotype(geno_donor, ref_data=self.ref_data, ard=self.ard)
                                            for geno_donor in geno_donors]
            elif isinstance(geno_donors[0], dict):
                geno_donors = [Genotype(geno_donor['genotype'], 
                                    id=geno_donor['id'], ref_data=self.ref_data, ard=self.ard)
                                    for geno_donor in geno_donors]
            # else:
            #     raise Exception('Please provide correctly formatted genotypes for the donor')
        return geno_recip, geno_donors

    def _assign_matches(self, matches : List[GenotypeMatch]):
        if matches:
            return matches
        matches = []
        if self.geno_recip and self.geno_donors:
            for geno_donor in self.geno_donors:
                matches.append(GenotypeMatch(self.geno_recip,
                                     geno_donor, ref_data=self.ref_data, ard=self.ard))
            return matches
        raise Exception("No matches provided.")

    def sort(self, 
            sort_feat_order : List[str] = ['expression', 'tce']) -> List[GenotypeMatch]:
        # match_grades = ['AA', 'AP', 'AL', 'AM', 'PP', 'LP', 'MP', 'LL', 'LM', 'MM']
        num_matches = [2, 1, 0]
        expr_matches = ['matched', 'favorable', 'unfavorable', 'mismatched', 'Unknown', None]
        tce_matches = ['Allotype', 'Permissive', 'nonpermissive', 'Unknown', None]
        # match_directionality = ['HvG', 'bidirectional', 'GvH', None]
        matches = [match for match in self.matches if isinstance(match, GenotypeMatch)]
        errors = [match for match in self.matches if not isinstance(match, GenotypeMatch)]
        df_sorter = pd.DataFrame([[match, match.id,
                                # match_grades.index(''.join(sorted(match.grade))),
                                num_matches.index(match.num_matches),
                                expr_matches.index(match.annotation['expr_match']), 
                                tce_matches.index(match.annotation['tce_match'].replace('HvG_', '').replace('GvH_', '')),
                                ]
                                    for match in matches],
                                columns=['match', 'id', 'num_matches', 'expression', 'tce'])
        sort_feat_order = ['num_matches'] + sort_feat_order
        df_sorter = df_sorter.sort_values(sort_feat_order)
        # Rows with the same sorting features get the same rank.
        df_sorter['rank'] = pd.factorize(df_sorter[sort_feat_order].apply(tuple, axis=1))[0] + 1
        if self.verbose:
            print(df_sorter)
        matches_sorted = []
        for i, row in df_sorter.iterrows():
            match = row['match']
            match.rank = row['rank']
            match.index = i + 1
            matches_sorted.append(match)
        # df_sorter['match'] = df_sorter['match'].astype(str)
        # df_sorter = df_sorter[sort_feat_order]
        # return df_sorter
        self.df_sorter = df_sorter
        self.matches_sorted = matches_sorted + errors
        return self.matches_sorted
            
    def __repr__(self):
        return "{}-:-{}".format(self.geno_recip, self.geno_donors)

    def serialize(self) -> List[Dict[str, any]]:
        return [match.serialize() for match in self.matches_sorted]