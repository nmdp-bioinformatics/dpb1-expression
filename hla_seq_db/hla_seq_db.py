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
from urllib.request import urlopen
import re
import hla_seq_db.hla_seq_db
from .ref_data import RefData
from .hml_parser import HmlParser
from .genotype_match import GenotypeMatch
from .match_list import MatchList
from .sequence import Sequence
from .allotype import Allotype
from .genotype import Genotype
from .annotator import Annotator
from .extractor import Extractor
from .individual import Individual
from .features import Features
from .individual_match import IndividualMatch
from .dataset import Dataset
import numpy as np
from Bio import pairwise2
from typing import Dict, List, Union
import cProfile, pstats, io
from pstats import SortKey
import time
import pandas as pd
from pyard import ARD
from hla_seq_db.sequence_search import SeqSearch

# Antigen-level
class HlaDB(object):
    def __init__(self, db_version : str = '3470', 
                       loci : List[str] = ['A', 'C', 'B', 'DRB1', 'DQB1', 'DQA1', 'DPB1'],
                       ard : ARD = None,
                       verbose : bool = False,
                    profile : bool = False):
        self.verbose = verbose
        if self.verbose:
            start_time = time.time()
        if profile:
            self.init_profile()
        if isinstance(loci, str):
            loci = [loci]
        self.loci = loci
        self.ref_allele : Allotype = None
        self.alleles = {}
        self.headers = []
        self.db_version = db_version.replace('.', '')
        self.ard = ard or ARD(self.db_version)
        self.ref_data = RefData(loci=loci, 
                                db_version=db_version, ard=self.ard,
                                verbose=verbose)
        self.ref_data_other = {}
        self.hml_parser = HmlParser()
        self.annotator = Annotator(self.ref_data, self.loci)
        self.extractor = Extractor(annotator=self.annotator, ard=self.ard)
        # self.annotate_data()
        if profile:
            self.print_stats()
        if self.verbose:
            print("Initialized in {:.2} seconds.".format(time.time() - start_time))

    def init_profile(self) -> None:
        self.pr = cProfile.Profile()
        self.pr.enable()
    
    def print_stats(self) -> None:
        """
        Prints profiling data
        """
        self.pr.disable()
        s = io.StringIO()
        sortby = SortKey.TIME
        ps = pstats.Stats(self.pr, stream=s).sort_stats(sortby)
        ps.print_stats(10)
        print(s.getvalue())


    # def get_seq(self, typing : str, feats: Union[str, List[str]] = None,
    #         gene_feats : bool = False,
    #         consensus : bool = False,
    #         verbose : bool = False) -> Union[pd.core.frame.DataFrame, str]:
    #     return self.extractor.extract(typing, 
    #         feats=feats,
    #         gene_feats=gene_feats,
    #         consensus=consensus,
    #         verbose=verbose)
    
    def get_map(self, locus : str,
            feature : str,
            group_col : str = 'allele_two_field',
            unk_val : str = 'unknown') -> Dict[str, str]:
        return self.ref_data.get_map(locus, group_col, feature, unk_val)

    def _load_ref_data(self, db_version : str) -> RefData:
        db_version = db_version.replace('.', '')
        if db_version == self.db_version:
            return self.ref_data
        if db_version not in self.ref_data_other:
            ref_data = RefData(loci=self.loci, 
                             db_version=db_version,
                             verbose=self.verbose)
            self.ref_data_other[db_version] = ref_data
        return self.ref_data_other[db_version]

    def get_features(self, allotype : str) -> Features:
        return self.ref_data.get_features(allotype)

    def annotate_allotype(self, allotype : Union[Allotype, str, pd.core.frame.DataFrame], 
            seq : str = None, db_version : str = None, call_alleles : bool = True,
            feats : Union[List[str], bool] = None, verbose : bool = False,
            field_level : int = 2) -> Allotype:
        db_version = db_version or self.db_version
        if isinstance(allotype, str) or isinstance(allotype, pd.core.frame.DataFrame):
            allotype = Allotype(allotype, seq=seq, verbose=verbose,
                         db_version=db_version,
                         ref_data=self._load_ref_data(db_version),
                            ard=self.ard, seq_diffs=True)
        # if feats or seq:
        if allotype.seq and call_alleles:
            allotype.call_alleles()
            # if not allotype.alleles_called:
            #     allotype.call_similar_alleles()
        allotype.get_features(feats=feats, field_level=field_level)
        if allotype.seq:
            # for seq in seq.split('|'):
            allotype.parse_gene_features(feats=feats)
            # if not allotype.feats_searched.seqs.feats:
            #     allotype.call_alleles()
            #     allotype.get_features(feats=feats)
            #     allotype.parse_gene_features(feats=feats)
            allotype.calc_and_set_snps()
        return allotype
    
    def parse_gene_features(self, allotype : Allotype, seq : str, feats : List[str] = None, verbose : bool = False,
            features_ref : Features = None) -> SeqSearch:
        # annotation = Annotation(allele,
        #      str(sequence), ref_data=self.ref_data, verbose=self.verbose)
        # annotation.annotate()
        if not isinstance(allotype, Allotype):
            allotype = self.create_allotype(allotype)
        if not features_ref or True:
            features_ref = self.ref_data.get_features(allotype)
        seq_search = SeqSearch(seq, 
                features_ref.seqs, ref_data=self.ref_data, verbose=verbose, 
                locus=allotype.locus,
                feats=feats)
        return seq_search

    def annotate_genotype(self, genotype : Union[Genotype, List[Allotype]],
            feats : List[str] = []):
        if (isinstance(genotype, str) or (isinstance(genotype, list) and len(genotype) == 2) or
            isinstance(genotype, pd.core.frame.DataFrame)):
            genotype = Genotype(genotype, ref_data=self.ref_data, ard=self.ard)
        for allele in genotype.allotypes:
            self.annotate_allotype(allele, feats=feats)
        genotype.get_annotations()
        return genotype

    def annotate_match(self, geno_recip : str = None, geno_donor : str = None,
                             match : GenotypeMatch = None, diffs : bool = False,
                             feats : List[str] = []) -> dict:
        """
        Annotates and returns matching information.
        :param geno_recip: genotype of the recipient
        :param geno_donor: genotype of the donor
        """
        if geno_recip and geno_donor:
            match = GenotypeMatch(geno_recip, geno_donor, ref_data=self.ref_data, ard=self.ard)
        elif not match:
            raise Exception("Please provide the genotypes of recipient and donor or a GenotypeMatch object.")
        feats = []
        if match.locus == 'DPB1':
            feats = ['tce', 'expression', 'expression_experimental']
        elif match.locus == 'B':
            feats = ['P2']
        match.genotype_patient = self.annotate_genotype(match.genotype_patient, feats=feats)
        match.genotype_donor = self.annotate_genotype(match.genotype_donor, feats=feats)
        match.annotate()
        if diffs:
            match.get_genotype_diffs()
        return match

    def annotate_individual(self, data : Union[str, pd.core.frame.DataFrame]) -> Individual:
        """
        An individual can be represented via a GL string (str) 
        or a table (pd.core.frame.DataFrame) of allotype, phase, and seq.
        """
        individual = Individual(data, ref_data=self.ref_data, ard=self.ard)
        individual.annotate()
        return individual

    def annotate_individual_match(self, data : Union[List[str], pd.core.frame.DataFrame] = None,
            patient : str = None, donor : str = None,
            loci : List[str] = ['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPB1'],
            align_mismatches_only : bool = True,
            verbose : bool = False) -> IndividualMatch:
        match = IndividualMatch(data=data, patient=patient, donor=donor,
             loci=loci, ref_data=self.ref_data, ard=self.ard, verbose=verbose)
        match.annotate(align_mismatches_only=align_mismatches_only)
        return match

    def annotate_match_list(self, geno_recip : str,
                                  geno_donors : List[Union[str, Dict[str, str]]],
                                  verbose : bool = False):
        """
        Annotates and sorts through stem cell sources for a recipient.
        """
        # self.annotate_match(match=GenotypeMatch)
        match_list = MatchList(geno_recip, geno_donors, ref_data=self.ref_data, ard=self.ard,
                                verbose=verbose)
        for match in match_list.matches:
            self.annotate_match(match=match)
        return match_list

    def annotate_dataset(self, dataset : Union[Dataset, pd.core.frame.DataFrame], seqs : pd.core.frame.DataFrame = None,
            annotations : List[str] = [], stack_loci : bool = False,
            loci : List[str] = None, nrows : int = None, n_print_out : int = None,
            output : str = None, verbose : bool = False):
        """
        Specific annotations can be provided in the output such as:
            'genotypes' : Columns of formatted genotypes such as 'R_DPB1_geno'.
            'tce_match' : Includes TCE matching in the 'tce_match' column.
            'tce' : TCE group assignments. An example column would be 'R_DPB1_TYP1_TCE'.
        """
        # if isinstance(dataset, str) or isinstance(dataset, pd.core.frame.DataFrame):
        ds = Dataset(dataset, stack_loci=stack_loci)
        ds.annotate(self, seqs=seqs, annotations=annotations, nrows=nrows, n_print_out=n_print_out,
                        loci=loci, output=output, verbose=verbose)
        return ds
        # return dataset
    
    def parse_hml(self, filepath : Union[str, List[str]] = None, xmlstring :str = None, output : str = None, n_files_out : int = 0,
                verbose : bool = False, stratify_loci : bool = False,
                query : str = None, db : str = None, n_files : str = None,
                pwd : Union[str, dict] = None) -> pd.core.frame.DataFrame:
        return self.hml_parser.process_hml(filepath=filepath, xmlstring=xmlstring, outdir=output,
            verbose=verbose, n_files_out=n_files_out, query=query, db=db, n_files=n_files,
            stratify_loci=stratify_loci, pwd=pwd)
        
    def parse_db(self, db, query : str = None, fileout : str = None) -> pd.core.frame.DataFrame:
        return self.hml_parser.parse_db(db, query=query, fileout=fileout)
    
    def create_allotype(self, typing : str):
        return Allotype(typing, ref_data=self.ref_data, ard=self.ard)
    
    def create_genotype(self, genotype : str, id : str = None, sire : str = None):
        return Genotype(genotype, id=id, sire=sire, ref_data=self.ref_data, ard=self.ard)

    def create_match(self,  genotype_patient: Genotype, 
            genotype_donor: Genotype) -> None:
        return GenotypeMatch(genotype_patient, genotype_donor, ref_data=self.ref_data, ard=self.ard)

    def create_dataset(self,  filepath : str = None,
            headers_subj : List[str] = ['R', 'D', 'RECIP', 'DONOR', 'RID', 'DID'],
            headers_locus : List[str] = ['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPB1'],
            headers_allele : List[str] = ['TYP', 'ALLOTYPE', '1', '2'],
            headers_seq : List[str] = ['SEQ'],
            headers_gene_feats : List[str] = ['UTR5', 'EXON_1', 'INTRON_1', 'EXON_2', 'INTRON_2',
                                                 'EXON_3', 'INTRON_3', 'EXON_4', 'INTRON_4',
                                                 'EXON_5', 'INTRON_5', 'EXON_6', 'INTRON_6', 
                                                 'EXON_7', 'INTRON_7', 'EXON_8', 'UTR3'],
            stack_loci : bool = False,
            loci : List[str] = ['A', 'C', 'B', 'DRB1', 'DQB1', 'DPB1']) -> Dataset:
        return Dataset(filepath=filepath, headers_subj=headers_subj, hla_db=self,
            headers_locus=headers_locus, headers_allele=headers_allele, headers_gene_feats=headers_gene_feats,
            headers_seq=headers_seq, stack_loci=stack_loci, loci=loci)

    def annotate_ref_data(self, modes : List[str], loci : List[str] = None) -> None:
        """
        Passes on annotation types to the Annotator class.
        Possible values include:
            'gfe' - Gene Feature Enumeration
            'CIWD' - Common, Intermediate, Well-Documented alleles
            'dpb1_expr' - HLA-DPB1 expression motifs
            'b_leader' - HLA-B leader peptide P2 alleles
        """
        self.annotator.annotate_ref_data(modes, loci=loci)
