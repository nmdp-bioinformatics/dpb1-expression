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
from .genotype import Genotype
from .ref_data import RefData
from .allotype_match import AllotypeMatch
from pyard.ard import ARD
from typing import Tuple, List, Dict

class GenotypeMatch(object):

    def __init__(self, genotype_patient: Genotype, genotype_donor: Genotype, ref_data : RefData = None,
                id : str = None,
                ard : ARD = None) -> None:
        """
        Represents the match status between two HLA-B genotypes as a 
        two-element list of AllotypeMatches. Also, aligns the genotypes
        with the match grades. The first genotype represents the patient/recipient while
        the second genotype represents the potential donor.

        :param genotype_patient: Genotype object
        :param genotypes_donors: List of Genotype objects
        """
        self.ref_data = ref_data
        self.ard = ard
        self.id = id
        if isinstance(genotype_patient, str) and isinstance(genotype_donor, str):
            self.ref_data = ref_data or RefData()
            genotype_patient = Genotype(str(genotype_patient), ref_data=self.ref_data, ard=self.ard)
            genotype_donor = Genotype(str(genotype_donor), ref_data=self.ref_data, ard=self.ard)
        self.genotype_patient = genotype_patient
        self.genotype_donor = genotype_donor
        self.locus = self._validate_locus()

        self.genotype_donor_flipped = False
        self.grade, self.allele_matches, self.num_matches = self._get_match()
        self.allotype_mismatches = None
        self.directionality = self._determine_directionality()
        self.annotation = None
        self.rank = self.index = None

        self.matched_alleles_pat, self.matched_alleles_don, \
             self.mismatched_alleles_pat, self.mismatched_alleles_don = self.organize_alleles()


    def _validate_locus(self) -> None:
        loci = set([self.genotype_patient.locus, self.genotype_donor.locus])
        if len(loci) == 1:
            return loci.pop()
        raise Exception('There is an appropriate amount of loci. Please include only one HLA locus.')

    def annotate(self) -> None:
        """
        Annotates the match with known mismatching models.
        """
        annotation = {}
        if self.locus == 'DPB1':
            annotation['tce_match'] = self.determine_TCE_match()
            annotation['expr_match'] = self.determine_expr_match()
        elif self.locus == 'B':
            annotation['b_leader'] = self.determine_leader_match()
        if annotation:
            self.annotation = annotation

    def _get_match(self) -> Tuple[str, List[AllotypeMatch], int]:
        """
        Compares the allotype matches between the forward and reverse versions of the genotypes
        to get the best match grade combination, which is whatever contains the highest match grade.

        (highest) A > P > L > M (lowest)

        The match grade combination is returned with the higher grade on the second element (right-most).
        The aligned genotypes are returned as well, in accordance with the match grade combination.
        
        Returns two HLA allotype match grades:
            AA - A genotype match occurs if both allotype pairs are allele matches.
            MA - A single mismatch constitutes one allele match and one mismatch.
            MM - A double mismatch occurs when both allotype pairs are mismatches.
            PA/PP/LP/MP - A potential single/double genotype match occurs when at least one pair has a potential match.
            LA/LP/LL/ML - A single/double allele genotype mismatch occurs when at least one pair has an allele mismatch.
        """
        g1_a, g1_b = self.genotype_patient.allotypes
        g2_a, g2_b = self.genotype_donor.allotypes
        match_for = [AllotypeMatch(g1_a, g2_a, ref_data=self.ref_data), AllotypeMatch(g1_b, g2_b, ref_data=self.ref_data)]
        match_rev = [AllotypeMatch(g1_a, g2_b, ref_data=self.ref_data), AllotypeMatch(g1_b, g2_a, ref_data=self.ref_data)]

        rank_for = [match_for[0].score, match_for[1].score]
        rank_rev = [match_rev[0].score, match_rev[1].score]

        matches = match_for
        if min(rank_for) < min(rank_rev):
            match_code = str(match_for[0]) + str(match_for[1])
        elif min(rank_for) > min(rank_rev):
            self.genotype_donor.flip()
            self.genotype_donor_flipped = True
            matches = match_rev
            match_code = str(match_rev[0]) + str(match_rev[1])
        else: #TODO: Think of examples for this. Might need to look at total score.
            match_code = str(match_for[0]) + str(match_for[1])
        num_matches = len([match for match in matches if match.matched])
        return match_code, matches, num_matches

    def _determine_directionality(self) -> str:
        """
        Determines the directionality of a match for single mismatches.
        GvH indicates a graft-versus-host mismatch vector while
        HvG indicates a host-versus-graft mismatch vector.
        :return: Match directionality
        """
        if self.grade in ['AA', 'PP', 'PA', 'AP']:
            return None
        if (not self.genotype_patient.homozygous and not self.genotype_donor.homozygous):
            return 'bidirectional'
        elif not self.genotype_patient.homozygous and self.genotype_donor.homozygous:
            return 'GvH'
        elif self.genotype_patient.homozygous and not self.genotype_donor.homozygous:
            return 'HvG'

    def determine_TCE_match(self) -> str:
        """
        Sets the TCE (T-Cell Epitope) match category. Assumes
        TCE has been assigned in the alleles already.
        """
        if self.grade in ['AA', 'AP', 'PA', 'PP']:
            return 'Allele'
        tce_pat = [allele.feats.anns.serialize()['tce'] for allele in self.genotype_patient.allotypes]
        tce_don = [allele.feats.anns.serialize()['tce'] for allele in self.genotype_donor.allotypes]
        unk_val = 'unknown'
        if (unk_val in tce_pat) or (unk_val in tce_don):
            return 'Unknown'
        if not any(tce_pat + tce_don):
            raise InvalidMatchError(self.genotype_patient,
                self.genotype_donor,
                'No TCE groups are assigned. Please assign first.')
        min_tce_pat = min(tce_pat)
        min_tce_don = min(tce_don)
        if min_tce_pat == min_tce_don:
            return 'Permissive'
        elif min_tce_pat < min_tce_don:
            return 'GvH_nonpermissive'
        else:
            return 'HvG_nonpermissive'

    def organize_alleles(self) -> Tuple[List[Allotype], List[Allotype], List[Allotype], List[Allotype]]:
        """
        Organizes alleles into four outputs:
            1: alleles that are matched in the patient.
            2: alleles that are matched in the donor.
            3: alleles that are mismatched in the patient.
            4: alleles that are mismatched in the donor.
        """
        matched_grades = ['A', 'P']
        matched_alleles_pat = [match.allele_one for match in self.allele_matches if match.match_grade in matched_grades]
        matched_alleles_don = [match.allele_two for match in self.allele_matches if match.match_grade in matched_grades]
        mismatched_alleles_pat = [match.allele_one for match in self.allele_matches if match.match_grade not in matched_grades]
        mismatched_alleles_don = [match.allele_two for match in self.allele_matches if match.match_grade not in matched_grades]
        return matched_alleles_pat, matched_alleles_don, mismatched_alleles_pat, mismatched_alleles_don

    def parse_as_str(self, el_list : List[Allotype]) -> List[str]:
        """
        Returns list of Allotype object as list of string.
        """
        return [str(el) for el in el_list]

    def determine_expr_match(self) -> str:
        """
        Determines the patient's mismatched allele's expression level.
        """
        mismatched_allele_pat_expr_level = None
        if len(self.mismatched_alleles_pat) == 1:
            mismatched_allele_pat_expr_level = self.mismatched_alleles_pat[0].feats.anns.serialize()['expression']
        # elif self.grade in ['MA', 'AM', 'MP', 'PM', 'LA', 'AL', 'LP', 'PL']:
        #     mismatched_allele_pat_expr_level = self.matched_alleles_pat[0].annotation['expression']
        if not mismatched_allele_pat_expr_level:
            if self.grade in ['AA', 'AP', 'PA', 'PP']:
                return 'matched'
            elif self.grade in ['MM', 'LM', 'ML', 'LL', 'MM']:
                return 'mismatched'
            return None
        cats = {'high' : 'unfavorable',
                'low' : 'favorable',
                'unknown' : 'unknown'}
        mismatched_allele_pat_expr_level = mismatched_allele_pat_expr_level.replace('?', '').replace('~', '')
        if mismatched_allele_pat_expr_level in cats:
            self.expr_match = cats[mismatched_allele_pat_expr_level]
        else:
            self.expr_match = mismatched_allele_pat_expr_level
        return self.expr_match

    def determine_leader_match(self) -> str:
        if self.grade in ['AA', 'AP', 'PA', 'PP']:
            return 'matched'
        elif self.grade in ['MM', 'LM', 'ML', 'LL', 'MM']:
            return 'mismatched'
        leader_match_status = ''.join(
            [allotype[0].feats.seq_anns.serialize()['P2'].replace('?', '')
                for allotype in [self.mismatched_alleles_pat,
                                 self.mismatched_alleles_don,
                                 self.matched_alleles_pat]
            ])
        self.leader_match_status = leader_match_status
        return leader_match_status

    def get_genotype_diffs(self, align_mismatches_only : bool = True) -> List[AllotypeMatch]:
        allotype_mismatches = []
        # diffs = None
        for allele_match in self.allele_matches:
            # try:
            diffs = allele_match.get_allotype_diffs(align_mismatches_only=align_mismatches_only)
            # except Exception as e:
            #     print(allele_match.__repr__(), e)    
            if diffs:
                allotype_mismatches.append(allele_match)
        self.allotype_mismatches = allotype_mismatches
        return self.allotype_mismatches
    
    def serialize(self) -> Dict[str, str]:
        """
        Export dictionary for exposing information in the API.
        """
        output = {'genotype_patient' : self.genotype_patient.serialize(),
                    'genotype_donor' : self.genotype_donor.serialize(),
                    'directionality' : self.directionality,
                    'genotype_donor_flipped' : self.genotype_donor_flipped,
                    'grade' : self.grade,
                    'matched_alleles_pat' : self.parse_as_str(self.matched_alleles_pat),
                    'matched_alleles_don' : self.parse_as_str(self.matched_alleles_don),
                    'mismatched_alleles_pat' : self.parse_as_str(self.mismatched_alleles_pat),
                    'mismatched_alleles_don' : self.parse_as_str(self.mismatched_alleles_don)}
        if self.annotation:
            output['annotation'] = self.annotation
        if self.rank:
            output['rank'] = self.rank
        if self.index:
            output['index'] = self.index
        if self.allele_matches:
            output['allele_matches'] = [match.serialize() for match in self.allele_matches]
        return output
    
    def __repr__(self) -> str:
        return "{}-{}".format(self.genotype_patient, self.genotype_donor)