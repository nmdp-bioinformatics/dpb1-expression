#
# Copyright (c) 2024 NMDP.
#
# This file is part of Aggregate Matching Tool
# (see https://github.com/nmdp-bioinformatics/agg-match-tool).
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
from .genotype import Genotype
from .haplotype import Haplotype
from typing import List
from .ref_data import RefData
from pyard import ARD
from typing import Union, Tuple, Dict
import pandas as pd
from .dataset import Dataset

class Individual(object):
    
    def __init__(self, data : Union[str, pd.core.frame.DataFrame], id : str = None,
            ref_data : RefData = None, ard : ARD = None, verbose : bool = False) -> None:
        """
        An individual can be represented via a GL string (str) 
        or a table (pd.core.frame.DataFrame) of allotype, phase, and seq.
        """
        self.verbose = verbose
        self.data = data
        self.ref_data = ref_data
        self.loci = ['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPB1']
        self.ard = ard
        self.id = id
        self.glstring, self.genos = self._parse_data()
        self.annotation = None
    
    def _parse_data(self) -> Tuple[str, List[Genotype]]:
        if isinstance(self.data, pd.core.frame.DataFrame):
            df = self.data.copy()
            if 'ALLOTYPE' in df:
                if 'id' in df:
                    id = set(df['id'])
                    if len(id) != 1:
                        raise Exception("You need to supply exactly one unique identifier in the 'id' column of the supplied DataFrame.")
                    self.id = id.pop()
                df.loc[:, 'LOCUS'] = df['ALLOTYPE'].str.split('*').str[0].str.replace('HLA-', '')
                df = df[df['LOCUS'].isin(self.loci)]
                df_genotypes = dict(tuple(df.groupby('LOCUS'))).values()
                genotypes = df_genotypes
            else:
                ds = Dataset(df)
                genotypes = [ds.get_dataframe(headers_locus=locus) for locus in ds.headers_locus]
        elif isinstance(self.data, str):
            glstring = self.data
            genotypes = glstring.split('^')
        genotypes = [Genotype(genotype, verbose=self.verbose, ref_data=self.ref_data, ard=self.ard) for genotype in genotypes]
        glstring = '^'.join([str(genotype) for genotype in genotypes])
        self.loci = [locus for locus in self.loci if locus in [geno.locus for geno in genotypes]]
        return glstring, genotypes

    def annotate(self) -> Dict[str, any]:
        if self.annotation: return self.annotation
        loci = {geno.locus : geno for geno in self.genos}
        results = {}
        if ('DQA1' in loci) and ('DQB1' in loci):
            results.update(self._annotate_dq_heterodimers(loci))
        self.annotation = results
        return self.annotation

    def _annotate_dq_heterodimers(self, loci : Dict[str, Genotype]) -> Dict[str, any]:
        dqa1 = loci['DQA1']
        dqb1 = loci['DQB1']
        g1a = ['2', '3', '4', '5', '6']
        g1b = ['2', '3', '4']
        g2a = ['1']
        g2b = ['5', '6']
        group_geno = []
        dq_molecules = {}
        for dqa1_allo in dqa1.allotypes:
            # For calculating the HLA-DQ Group genotype
            group = None
            if (dqa1_allo.allele_family in g1a):
                group = 'G1'
            elif (dqa1_allo.allele_family in g2a):
                group = 'G2'
            else:
                raise Exception('The DQA1 typing does not fall into known HLA-DQ groupings.')
            group_geno.append(group)

            for dqb1_allo in dqb1.allotypes:
                # Group calculation for checking G1, G2 combinations
                group = None
                if (dqa1_allo.allele_family in g1a) and (dqb1_allo.allele_family in g1b):
                    group = 'G1'
                elif (dqa1_allo.allele_family in g2a) and (dqb1_allo.allele_family in g2b):
                    group = 'G2'
                # else:
                #     raise Exception('This combination is neither G1 or G2', 
                #                 (dqa1_allo.allele_family, dqb1_allo.allele_family))
                if group:
                    # Num DQ molecule calculation
                    # if (len(dqa1_allo.alleles_hi_res) != 1) or (len(dqb1_allo.alleles_hi_res) != 1):
                    #     raise Exception('The calculated HLA-DQ molecules are ambiguous', 
                    #         (dqa1_allo.alleles_hi_res, dqb1_allo.alleles_hi_res))
                    dq_molecule = Haplotype([dqa1_allo, dqb1_allo],
                                    ref_data=self.ref_data, ard=self.ard)
                    dq_molecule.annotation['dq_group'] = group
                    if dq_molecule.glstring not in dq_molecules:
                        dq_molecules[dq_molecule.glstring] = dq_molecule
        result = {'num_dqs' : len(dq_molecules),
                  'dq_groups' : '+'.join(group_geno),
                  'dq_heterodimers' : list(dq_molecules.values())}
        return result

    def serialize(self) -> Dict[str, any]:
        output = {'glstring' : self.glstring}
        output['genotypes'] = [geno.serialize() for geno in self.genos]
        if self.annotation:
            output['annotation'] = self.annotation
        return output

    def __repr__(self):
        return self.glstring