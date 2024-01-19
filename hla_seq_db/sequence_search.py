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
from Bio.Seq import Seq
from Bio import pairwise2
from hla_seq_db.snp import Snp
from typing import Tuple, Union, List, Dict
from hla_seq_db.sequence_features import SeqFeatures
from Bio.SeqFeature import SeqFeature, FeatureLocation
from hla_seq_db.sequence_match import SeqMatch
from hla_seq_db.sequence_feat_searched import SeqFeatSearched
import copy
import time
# import hla_seq_db
# from hla_seq_db.ref_data import RefData``
# from hla_seq_db.extractor import Extractor
import re
import regex

class SeqSearch(object):
    
    def __init__(self,
                 seq : Sequence, seq_features_ref : SeqFeatures,
                 seq_features : SeqFeatures = None,
                 feats : List[str] = None,
                 locus : str = None,
                 name : str = None, ref_data = None, exact : bool = True,
                 verbose : bool = False):
        self.exact = exact
        if seq and isinstance(seq, str):
            seq = Sequence(seq, name=name)
        if seq_features_ref and isinstance(seq_features_ref, dict):
            seq_features_ref = SeqFeatures(seq_features_ref)
        self.verbose = verbose
        self.locus = locus
        self.seq = seq
        self.seq_features_ref = seq_features_ref
        self.gene_feat_search_order = ['exon_2', 'exon_3', 'exon_1', 'exon_4', 'exon_5', 'exon_6', 'exon_7', 'intron_7', 'exon_8',
                                       'intron_2', 'intron_3', 'intron_1', 'intron_4', 'intron_5', 'intron_6', 
                                       'utr5', 'utr3']
        self.seqs = [seq]
        self.seqs_searched = []
        self.locations_searched = []
        self.ref_data = ref_data
        # self.extractor = extractor
        self.feats = feats
        self.seq_features = seq_features
        if not self.seq_features:
            self.seq_features = self._search_seq_features(self.seq, 
                self.seq_features_ref, feats)
        self.feat_names = self._calc_feat_names(self.seq_features['seqs'].feats)
        # self.gfe_nt_dist, self.gfe_aa_dist = self._calc_gfe(self.seq_features['gene_features'].feats)
        self.snps = self._parse_out_snps(self.seq_features)
    
    def _calc_feat_names(self, seq_features : Dict[str, SeqFeatSearched]) -> str:
        feat_names = []
        for seq_name, seq_feat in seq_features.items():
            prefix = 'p_' if seq_feat.partial else ''
            feat_names.append(prefix + seq_name)
        return '~'.join(feat_names)

    def _search_seq_features(self, seq : Sequence, seq_features : SeqFeatures,
            feats : List[str]) -> Union[Dict[str, SeqFeatures], SeqFeatures]:
        gene_feats = self._search_gene_features()
        results = {'seqs' : gene_feats}
        if feats:
            ann_feats = self._search_ann_features(gene_feats, feats)
            results['seq_anns'] = ann_feats
        return results
    
    def _search_ann_features(self, gene_feats : SeqFeatures, feats : List[str]) -> SeqFeatures:
        ann_feats = {}
        for feat in feats:
            if self.locus:
                motif = self.ref_data.alleles[self.locus].get_motif(feat)
                if motif and motif.type in gene_feats.feats:
                    seq_feat_searched = gene_feats.feats[motif.type]
                    ann_feat_searched = seq_feat_searched.extract(motif, name=feat)
                    ann_feats[feat] = ann_feat_searched
        return SeqFeatures(ann_feats)

    
    def _search_gene_features(self, 
            features : Dict[str, SeqFeatSearched] = {},
            wildcard : bool = False) -> SeqFeatures:

        skipped_feat_lists = [['utr5'], None]
        fuzzy = True
        finished = False
        if self.verbose:
            print('Searching sequence with {:,} nts'.format(len(self.seq)))
        for skipped_feat_list in skipped_feat_lists:
            if self.verbose:
                print('Skipping searches on these feats:', skipped_feat_list)
            # if not finished and not self.exact:
            #     if self.verbose:
            #         print('Searching features via wildcard.')
            #     finished = self._search_seq(skipped_feats=skipped_feat_list, wildcard=True)
            if not finished:
                if self.verbose:
                    print('Searching core features')
                finished = self._search_seq(skipped_feats=skipped_feat_list)
            if not finished:
                if self.verbose:
                    print('Searching partial flanking seq features.')
                    print('Unsearched seqs:', [seq for seq in self.seqs if not seq.searched])
                finished = self._search_seq(partial_seq=True, flanking_seq=True, fuzzy=fuzzy, skipped_feats=skipped_feat_list)
            if not finished:
                if self.verbose:
                    print('Searching flanking features')
                    print('Unsearched seqs:', [seq for seq in self.seqs if not seq.searched])
                finished = self._search_seq(flanking_seq=True, fuzzy=fuzzy, skipped_feats=skipped_feat_list)
            if not finished:
                if self.verbose:
                    print('Searching partial seq features.')
                    print('Unsearched seqs:', [seq for seq in self.seqs if not seq.searched])
                finished = self._search_seq(partial_seq=True, fuzzy=fuzzy, skipped_feats=skipped_feat_list)
            if self.verbose:
                print('Finished?', finished)
                print('Unsearched seqs:', [seq for seq in self.seqs if not seq.searched])
        
        # Clear search history
        for feat_seq_ref in self.seq_features_ref.feats.values():
            feat_seq_ref.searched = False
        # if feats_unsearched:
        #     if feats_searched:
        #         if self.verbose:
        #             print('Features unsearched:', feats_unsearched)
        #             print('Searching flanking features')
        #         features = self._search_flanking_seq_features(seq, seq_features,
        #                         features, feats_unsearched, feats_unsearched_ref)
        #         features, feats_searched, feats_unsearched, feats_unsearched_ref = self._order_and_search_unsearched(seq, seq_features, features=features)

        #         if self.verbose: 
        #             print('Features unsearched:', feats_unsearched)
        #             print('Imputing feature locations')
        #         features = self._impute_features(seq_features, features,
        #                          feats_unsearched)
        #         features, feats_searched, feats_unsearched, feats_unsearched_ref = self._order_and_search_unsearched(seq, seq_features, features=features)
                
        #     if feats_unsearched:
        #         if self.verbose:
        #             print('Features unsearched:', feats_unsearched)
        #             print('Searching partial seq features.')
        #         features = self._search_partial_seq_features(seq, seq_features,
        #                         features, feats_unsearched)
        #         features, feats_searched, feats_unsearched, feats_unsearched_ref = self._order_and_search_unsearched(seq, seq_features, features=features)
        #         if feats_unsearched and not wildcard:
        #             if self.verbose:
        #                 print('Features unsearched:', feats_unsearched)
        #                 print('Obtaining wildcard', feats_unsearched)
        #             seq_features = self.ref_data.alleles[self.locus].get_features('DPB1*XX:XX').seqs
        #             if self.verbose: print('Wildcard obtained')
        #             features = self._search_gene_features(seq, seq_features, features=features, wildcard=True).feats
        #             features, feats_searched, feats_unsearched, feats_unsearched_ref = self._order_and_search_unsearched(seq, seq_features, features=features)

        #         if feats_unsearched:
        #             raise Exception('Not completely parsed', [feat.location for feat in feats_unsearched])
        
        features = {seq_match.name : seq_match
               for seq_match in self.seqs 
               if seq_match.searched}
        return SeqFeatures(features)
        
    
    # def _search_flanking_seq_features(self,
    #         seq : Sequence, seq_features : SeqFeatures,
    #         features : Dict[str, SeqFeatSearched],
    #         feats_unsearched : List[Sequence], feats_unsearched_ref : List[Sequence]) -> Dict[str, SeqFeatSearched]:
    #     # Search flanking sequences
    #     indices_searched = [feat.index for feat in features.values() if feat.searched]
    #     for feat_unsearched in feats_unsearched:
    #         for feat_unsearched_ref in feats_unsearched_ref:
    #             index = feat_unsearched_ref.index
    #             # print(feat_unsearched_ref.name, index < min(indices_searched), index > max(indices_searched))
    #             feature = self._search_seq(feat_unsearched_ref, feat_unsearched,
    #                 flank_start=index < min(indices_searched),
    #                   flank_end=index > max(indices_searched))
    #             if feature:
    #                 # print(feature)
    #                 features[feature.name] = feature
    #                 features, feats_searched, feats_unsearched, feats_unsearched_ref = self._order_and_search_unsearched(seq, seq_features, features=features)
    #                 if not feats_unsearched:
    #                     break
    #     return features

    # def _search_partial_seq_features(self, 
    #         seq : Sequence, seq_features : SeqFeatures,
    #         features : Dict[str, SeqFeatSearched], feats_unsearched : List[Sequence]) -> Dict[str, SeqFeatSearched]:
    #     # Search core sequences
    #     feats_unsearched_ref = []
    #     for feat_unsearched in feats_unsearched:
    #         for feat_name, feat_seq in seq_features.feats.items():
    #             if feat_seq.gene_feat and not feat_seq.searched:
    #                 feature = self._search_seq(feat_seq, feat_unsearched, verbose=feat_name == 'exon_3',
    #                     partial_seq=True)
    #                 # if feat_name == 'exon_3':
    #                 #     print(feat_seq, feat_unsearched)
    #                 if feature:
    #                     if self.verbose:
    #                         print('Successfully found partial seq', feat_name)
    #                     features[feat_name] = feature
                # self._search_seq(feat_seq, )
                # if feat not in features:
                #     feature = None
                #     if len(seq) > len(feat_seq):
                #         feature = self._search_seq(seq, feat_seq)
                #     else:
                #         print('OTHER')
                #         feature = self._search_seq(feat_seq, seq)
                #     if feature:
                #         features[feat] = feature
                #     else:
                #         features[feat] = feat_seq
                        # feats_unsearched_ref.append(feat_seq)
        # return features
    
    def _parse_out_snps(self, seq_features : Dict[str, SeqFeatures]) -> Dict[str, List[Snp]]:
        return {feat_name : feature.snps
            for feat_type, seq_feature in seq_features.items()
            for feat_name, feature in seq_feature.feats.items() 
            if feature.snps}
    
    def _impute_features(self, seq_features : SeqFeatures,
            features : Dict[str, SeqFeatSearched],
            feats_unsearched : List[Sequence]) -> Dict[str, SeqFeatSearched]:
        locations = [feat.location for feat in features.values() if feat.searched]
        for feat_unsearched in feats_unsearched:
            loc_unsearched = feat_unsearched.location
            index_before = index_imputed = index_after = feat_imputed = None
            for feat_name, feat in features.items():
                if feat.searched and not index_imputed:
                    loc_searched = feat.location
                    if loc_searched.end == loc_unsearched.start:
                        index_before = feat.index
                    elif loc_searched.start == loc_unsearched.end:
                        index_after = feat.index
                    if index_before and index_after:
                        if index_after - index_before == 2:
                            index_imputed = index_after - 1
                if index_before and (feat.index == index_before + 1):
                    feat_imputed = feat_name
            if index_imputed:
                if self.verbose:
                    print(feat_imputed, 'was imputed.')
                seq_features = self.ref_data.alleles[self.locus].get_features('DPB1*XX:XX').seqs
                feature_ref_wildcard = seq_features.feats[feat_imputed]
                features[feat_imputed] = SeqFeatSearched(None,
                    self._align_seq_w_spacers(feature_ref_wildcard, feat_unsearched),
                                            location=feat_unsearched.location,
                                            index=index_imputed,
                                            searched=True, imputed=True)
        #     [feat.location for feat in features.values()
        #          if feat.searched and loc.start == feat.location.end]
        #     print(locations)
        return features
    
    # def _calc_unsearched_features(self, 
    #         seq : Sequence, 
    #         features : Dict[str, SeqFeatSearched]) -> Tuple[List[Sequence], List[Sequence]]:
    #     """
    #     Calculates uncovered features and filters down the possible features.
    #     """
    #     location_start = None
    #     locations_unsearched = []
    #     matched_feat = None
    #     prev_loc_end = 0
    #     prev_loc_start = 0
    #     n = 0
    #     feat_potential_matches = [] #feat_potential_matches
    #     feat_potential_matches_cache = []
    #     for i, (feat_name, seq_match) in enumerate(features.items()):
    #         last_feat = (i == len(features) - 1)
    #         if seq_match.searched:
    #             index = seq_match.seq_one.index
    #             loc = seq_match.location
    #             if prev_loc_end or (n == 0 and loc.start != 0):
    #                 if prev_loc_end > loc.start:
    #                     raise Exception('Overlapping features')
    #                     # print('Overlapping features')
    #                 elif prev_loc_end != loc.start:
    #                     feat_potential_matches += feat_potential_matches_cache
    #                     feat_potential_matches_cache = []
    #                     locations_unsearched.append(FeatureLocation(
    #                         prev_loc_end,
    #                         loc.start
    #                     ))
    #             prev_loc_end = loc.end
    #             n += 1
    #         else:
    #             feat_potential_matches_cache.append(seq_match)
    #             if last_feat:
    #                 feat_potential_matches += feat_potential_matches_cache
    #     if prev_loc_end < len(seq):
    #         locations_unsearched.append(FeatureLocation(prev_loc_end, len(seq)))
    #     feat_potential_matches.reverse()
    #     seqs_unsearched = [Sequence(location.extract(seq), location=location)
    #                         for location in locations_unsearched]
    #     return seqs_unsearched, feat_potential_matches

    def _search_seq(self, partial_seq : bool = False,
            flanking_seq : bool = False, fuzzy : bool = False, skipped_feats : bool = ['utr5', 'utr3'],
            verbose : bool = False, wildcard : bool = False) -> bool:
        # for index, (feat_name, feat_seq_ref) in enumerate(self.seq_features_ref.feats.items()):
        if wildcard:
            seq_features_ref = self.ref_data.alleles[self.locus].get_features('XX:XX').seqs
            for gene_feat in seq_features_ref.feats.values():
                gene_feat.searched = False
        else:
            seq_features_ref = self.seq_features_ref
        for feat_name in self.gene_feat_search_order:
            if feat_name in seq_features_ref.feats:
                feat_seq_ref = seq_features_ref.feats[feat_name]
                index = feat_seq_ref.index
                for i, seq in enumerate(self.seqs):
                    if not feat_seq_ref.searched and not seq.searched:
                        # print('seq', seq)
                        if skipped_feats and feat_seq_ref.name in skipped_feats:
                            continue
                        flank_start = flank_end = False
                        if flanking_seq:
                            indices_searched = [seq.index for seq in self.seqs if seq.index != None]
                            if indices_searched:
                                flank_start = index < min(indices_searched)
                                flank_end = index > max(indices_searched)
                            # print(feat_name, flank_start, flank_end)
                        # if (i == 0) or (index == 0):
                        index = feat_seq_ref.index
                        if self.verbose:
                            start_time = time.time()
                        # print(feat_name, i, len(seq.seq.replace('.', '')), len(feat_seq_ref.seq.replace('.', '')))
                        # if len(seq.seq.replace('.' ,'')) < len(feat_seq_ref.seq.replace('.', '')) and not partial_seq:
                        #     continue
                        # if len(seq) > len(feat_seq_ref):
                        #     continue asdf
                        seq_feat_searched = self._find_seq_feat_searched(seq, feat_seq_ref,
                                    partial_seq=partial_seq, 
                                    flank_start=flank_start, 
                                    flank_end=flank_end, 
                                    fuzzy=fuzzy,
                                    verbose=verbose)
                        if self.verbose:
                            dur = time.time() - start_time
                            if dur > 1:
                                print('{} took longer than expected: {:.2f} secs'.format(feat_seq_ref.name, dur))
                        if seq_feat_searched:
                            feat_seq_ref.searched = True
                            self.seqs_searched.append(seq_feat_searched)
                            # print('self.seqs', self.seqs)
                            # print('location', seq_feat_searched.location)
                            self.seqs = self._partition_seqs()
                            self.seqs = self._impute_seqs()
                            finished = not bool(len([seq for seq in self.seqs if not seq.searched]))
                            # print('self.seqs partiioned', self.seqs)
                            if finished:
                                return finished
                            continue
        return False

        # seq_feat_searched = self._find_seq_feat_searched(seq_feat, partial_seq=partial_seq, flank_start=flank_start, flank_end=flank_end, verbose=verbose)
        # asdf
        # if seq_feat_searched:
            

        #     print(self.seqs)
        #     print(seq_feat_searched.location.extract(seq_a))

    def _impute_seqs(self) -> List[Union[Sequence, SeqFeatSearched]]:
        feat_indices = {feat.index : feat for feat_name, feat in self.seq_features_ref.feats.items()}
        feat_names =   {feat.name : feat  for feat_name, feat in self.seq_features_ref.feats.items()}
        seqs = self.seqs
        for i in range(1, len(self.seqs) - 1):
            seq = self.seqs[i]
            if not seq.searched:
                seq_before = self.seqs[i - 1]
                seq_after = self.seqs[i + 1]
                if ((seq_before.index != None) and (seq_after.index != None) and
                    ((seq_before.index + 2) == seq_after.index)):
                        index = feat_names[seq_before.name].index + 1
                        feat_bridging = feat_indices[index]
                        location = seq.location
                        seq = self._align_seq_w_spacers(str(feat_bridging), seq)
                        seq_feat = Sequence(seq, name=feat_bridging.name)
                        seq_feat_searched = SeqFeatSearched(feat_bridging, seq_feat, index=index,
                            location=location, searched=True)
                        if self.verbose:
                            print('Successfully imputed {} at {} - {}'.format(
                                feat_bridging.name, location, seq_feat_searched.seq_two.__repr__()
                                ))
                        seqs[i] = seq_feat_searched
                        self.seqs_searched.append(seq_feat_searched)
        return seqs
        # prev_index = 0
        # for seq in self.seqs:
        #     if seq.index 
        #     prev_index = seq.index
        #     print(seq.index, seq.name, seq.searched)
        # return self.seqs

    def _partition_seqs(self) -> List[Union[Sequence, SeqFeatSearched]]:
        locations_unsearched = []
        prev_loc_end = 0
        prev_loc_start = 0
        i = 0
        seqs = []
        self.seqs_searched.sort(key=lambda x: x.index)
        # seq_searched = [seq_searched for seq_searched in self.seqs_searched if seq_searched.location.start == 0]
        # if len(seq_searched):
        #     seqs.append(seq_searched)
        for seq_searched in self.seqs_searched:
            # print(feat_name, index, loc)
            loc = seq_searched.location
            if prev_loc_end or (i == 0 and loc.start != 0):
                if prev_loc_end > loc.start:
                    raise Exception('Overlapping features', prev_loc_end, loc.start)
                    # print('Overlapping features', prev_loc_end, loc.start)
                elif prev_loc_end != loc.start:
                    location_unsearched = FeatureLocation(prev_loc_end, loc.start)
                    seqs.append(Sequence(location_unsearched.extract(self.seq), location=location_unsearched))
                    # locations_unsearched.append(FeatureLocation(
                    #     prev_loc_end,
                    #     loc.start
                    # ))
            seqs.append(seq_searched)
            prev_loc_end = loc.end
            i += 1
        if prev_loc_end < len(self.seq):
            location_unsearched = FeatureLocation(prev_loc_end, len(self.seq))
            seqs.append(Sequence(location_unsearched.extract(self.seq), location=location_unsearched))
            # locations_unsearched.append(FeatureLocation(prev_loc_end, len(self.seq)))
            # seqs.
        return seqs
        # return locations_unsearched

    def _find_seq_feat_searched(self, seq_a : Sequence, seq_b : Sequence, partial_seq : bool = False,
            flank_start : bool = False, flank_end : bool = False, fuzzy : bool = False, verbose : bool = False) -> Union[None, SeqFeatSearched]:
        """
        Search for sequence b (reference) within sequence a. Sequence b is assumed to be 
        an IMGT reference sequence unless this is a partial search (flank_start or flank_end)
        where the search is instead looking for sequence b within sequence a (reference).
        """
        if Sequence(str(seq_b).replace(' ', '')).empty: return None
        # print('START', seq_a, seq_b, seq_b.__dict__)

        if flank_start and flank_end:
            raise Exception('Cannot have two flanks (start and end) enabled.')
        # flipped = flank_start or flank_end
        flipped = partial_seq
        if flipped:
            seq = seq_a
            seq_a = seq_b
            seq_b = seq
        seq_b_regex = self.create_regex(seq_b, fuzzy=fuzzy, flank_start=flank_start, flank_end=flank_end, partial_seq=partial_seq)
        name = seq_b.name or seq_a.name
        if partial_seq:
            seq_b_regex = self.create_regex(seq_b, fuzzy=fuzzy, spacers=True, flank_start=flank_start, flank_end=flank_end)
            for seq in seq_a.seqs:
                m = regex.search(seq_b_regex, str(seq)) #raw_seq()
                # if (name == 'intron_2' and len(seq_b) > 3000):
                #     print(seq_b)
                #     print('seq_b_regex', seq_b_regex)
                #     print('seq', seq)
                #     print(m)
                if m:
                    # f = SeqFeature(location)
                    # seq_feat = Sequence(f.extract(seq))
                    location = seq_b.location or FeatureLocation(0, len(self.seq))
                    seq_feat = self._add_flanking_spacers(seq_a, str(seq_b), 
                            location=FeatureLocation(m.span()[0], m.span()[1]), 
                            flank_both=True)
                    verbose = False
                    # if ((name == 'intron_2')):
                    #     print('seq_b_regex', seq_b_regex)
                    #     print('seq', seq)
                    #     print('m', m)
                        # print('f', FeatureLocation(m.span()[0], m.span()[1]).extract(seq))
                    #     print('c', seq_feat)
                        # print('seq_a', seq_a)
                        # print('seq_feat', seq_feat)
                        # verbose = True
                    seq_feat = self._align_seq_w_spacers(seq_a, seq_feat, verbose=verbose)
                    # if (name == 'exon_3'):
                    #     print('d', seq_feat)
                    if self.verbose:
                        print('Successfully found partial seq {} at {}'.format(name, location))
                    return SeqFeatSearched(seq_a, seq_feat, index=seq_a.index, location=location, searched=True)
        else:
            # if (name == 'intron_2'):
            #     # print('a', seq_b)
            #     print('seq_b_regex', seq_b_regex)
            #     print('seq_a.raw_seq()', seq_a.raw_seq()) # asdf
            m = regex.search(seq_b_regex, seq_a.raw_seq()) #.replace('.', ''))
        if m:
            pos_one, pos_two = m.span()
            location = FeatureLocation(pos_one, pos_two)
            f = SeqFeature(location) #    type=seq_b.name)
            seq_feat_ref_full = seq_feat_ref = seq_b if not flipped else seq_a
            if flipped:
                seq = seq_b
                seq_feat_ref = seq_a #Sequence(f.extract(seq_a), name=name)
                # print(seq_a)
                # print(seq_b)
            else:
                seq = f.extract(seq_a)
                # print(str(seq_feat_ref))
                # print(seq)
            # if name == 'utr3':
            #     print('a', seq_feat_ref)
            #     print('b', seq)
            #     print(location)
            seq = self._align_seq_w_spacers(str(seq_feat_ref), seq)
            if flipped:
                seq = self._add_flanking_spacers(seq_feat_ref_full, seq, flank_start=flank_start)
            seq_feat = Sequence(seq, name=name, index=seq_b.index,
                location=location, searched=True)
            if flank_start:
                location = FeatureLocation(0, len(seq_b))
            elif flank_end:
                location = FeatureLocation(len(self.seq) - len(seq_b), len(self.seq))            
            if seq_a.location:
                location += seq_a.location.start
            # print(name, seq_b_regex)
            # print(location, pos_one, pos_two)
            # print('a', seq_feat_ref_full)
            # print('b', seq_feat)
            # print(seq_b_regex)
            # print(m)
            self.locations_searched.append(location)
            seq_feat_searched = SeqFeatSearched(seq_feat_ref_full, seq_feat,
                                   location=location, searched=True)
            # if (name == 'exon_4'):
            #     print('seq_b_regex', seq_b_regex)
            #     print('seq_b', seq_b)
            #     print('seq_a', seq_a.raw_seq())
            #     print(pos_one, pos_two)
            #     print(m)
            #     print('c', seq)
            #     print('c', seq_feat_searched.seq_one)
            #     print('d', str(seq_a.raw_seq()))
            if self.verbose:
                print('Successfully found {} at {} - {}'.format(name, location, str(seq_feat_searched.seq_two.__repr__())), )
            return seq_feat_searched
        return None

    def _add_flanking_spacers(self, seq_full : Sequence, seq_slice : str,
            flank_start : bool = False, flank_both : bool = False, location : FeatureLocation = None) -> str:
        """
        Adds flanking spacers ('.') to a seq_slice string to match the length
        of a full sequence (seq_full). Flank_end is assumed to be true if
        flank_start is false.
        """
        if location:
            seq = (location.start * '*') + seq_slice + ((len(seq_full) - location.end) * '*')
            return seq
        else:
            spacers = '*' * (len(seq_full) - len(seq_slice))
            if flank_start:
                seq_slice = spacers + seq_slice
            else:
                seq_slice = seq_slice + spacers
            return seq_slice
        
    def _get_n_nts_regex(self, seq : str, n : int, regex : bool = False) -> Union[str, None]:
        i = 0
        result = ""
        for nt in seq:
            if i < n:
                if nt not in ['.', ' ']:
                    i += 1
                result += nt
        if len(result) < n:
            return None
        if regex:
            return self._format_substring(result).replace(' ','')
        return result

    def _compare_lookahead(self, seq1 : str, seq2 : str, n_lookahead : int,
            shift1 : int = 0, shift2 : int = 0) -> bool:
        seq_ahead_ref = self._get_n_nts_regex(seq1[shift1:], n_lookahead, regex=True)
        seq_ahead_ref_regex = '({}){{s<=1}}'.format(seq_ahead_ref)
        seq_ahead =  self._get_n_nts_regex(seq2[shift2:], n_lookahead)
        # if len(seq_ahead_ref) == (n_lookahead - 1):
        if seq_ahead_ref and seq_ahead:
            # print(seq_ahead, seq_ahead_ref, bool(regex.match(seq_ahead_ref_regex, seq_ahead)), shift1, shift2)
            return bool(regex.match(seq_ahead_ref_regex, seq_ahead))
        return False

    def _align_seq_w_spacers(self, seq_ref : str, seq : str, reverse : bool = False, verbose : bool = False) -> str:
        """
        Aligns feat_seq with subseq to introduce spacers where required.
        """
        feat_seq = seq_ref
        subseq = seq
        seq_aligned = ""
        # print('a', feat_seq)
        # print('b', subseq)
        errata = []
        feat_seq = feat_seq[::(-1 if reverse else 1)]
        subseq = subseq[::(-1 if reverse else 1)]
        subseq_finished = False

        m = re.match('\*+', str(seq))
        if m:
            prefix_len = len(m.group())
            if (str(seq_ref[prefix_len:]).replace('.', '') == str(seq[prefix_len:])):
                return prefix_len * '*' + str(seq_ref[prefix_len:])

        # alignments = pairwise2.align.globalms(str(seq_ref), str(seq), 5, -4, -.5, -.1)
        # seq_ref = alignments[0][0]
        # seq = alignments[0][1]
        # seq_aligned = ""
        # spacers = ""
        # spaces = ""
        # # print(pairwise2.format_alignment(*alignments[0]))
        # for i, (nt_ref, nt) in enumerate(zip(seq_ref, seq)):
        #     if nt_ref != '-':
        #         if nt_ref == nt:
        #             seq_aligned += nt_ref
        #         else:
        #             if nt_ref in [' ', '.']:
        #                 seq_aligned += nt_ref
        #             else:
        #                 j = i
        #                 nt_next = seq[j]
        #                 # print(nt_ref, nt, seq_aligned)
        #                 while nt_next == '-' and j < len(seq):
        #                     # print(j)
        #                     # print(seq[j])
        #                     nt_next = seq[j]
        #                     j += 1
        #                 seq_aligned += nt_next
        #         # print(nt_ref, nt, seq_aligned)
        # return seq_aligned

        for i, nt in enumerate(feat_seq):
            # if verbose:
            #     print(i)
            if subseq:
                n_lookahead = 8
                matched = ((nt == subseq[0]) or (nt == 'X') or ('*' == subseq[0]))
                # print(nt, subseq[0], matched, seq_aligned)
                if nt == ' ':
                    seq_aligned += ' '
                elif nt == '*':
                    seq_aligned += '*'
                elif matched:
                    seq_aligned += subseq[0]
                    subseq = subseq[1:]
                else:
                    if self._compare_lookahead(feat_seq[i:], str(subseq), n_lookahead):
                        if nt == '.':
                            seq_aligned += '.'
                        else:
                            seq_aligned += subseq[0]
                            subseq = subseq[1:]
                    elif self._compare_lookahead(feat_seq[i:], str(subseq), n_lookahead, shift1=1):
                        # seq_aligned += '.'
                        # seq_aligned += subseq[0]
                        # subseq = subseq[1:]
                        # if nt == '.':
                        # seq_aligned += subseq[0]
                        # subseq = subseq[1:]
                        # else:
                        seq_aligned += '.'
                    elif self._compare_lookahead(feat_seq[i:], str(subseq), n_lookahead, shift2=1):
                        # print('D')
                        seq_aligned += subseq[0]
                        subseq = subseq[1:]
                    else:
                        if subseq[0] == '*':
                            seq_aligned += '*'
                            subseq = subseq[1:]
                        elif nt == '.':
                            seq_aligned += '.'
                        else:
                            seq_aligned += subseq[0]
                            subseq = subseq[1:]
                    # elsif:
                    #     seq_aligned += subseq[0]
                    #     subseq = subseq[1:]
            else:
                if nt == ' ':
                    seq_aligned += ' '
                elif nt == '*':
                    seq_aligned += '*'
                else:
                    seq_aligned += '.'
            # if (not matched and
            #     (str(subseq).replace('.', '')[:n_lookahead] == feat_seq.replace('.', '')[i:i + n_lookahead])):
            #         seq_aligned += subseq[0]
            #         subseq = subseq[1:]
            # elif nt in ['.', ' ']:
            #     if subseq:
            #         if nt == subseq[0]:
            #             seq_aligned += subseq[0]
            #             subseq = subseq[1:]
            #         else:
            #             seq_aligned += nt
            # elif subseq:
            #     # if 3932 <= i and i < 3940:
            #     #     print('nt', nt)
            #     #     print('nt', subseq[0])
            #     #     print('aligned', seq_aligned[-10:])
            #     #     print(feat_seq.replace('.', '')[i: i + n_lookahead])
            #     #     print(str(subseq).replace('.', '')[:n_lookahead])
            #     if (not matched and 
            #         (str(subseq).replace('.', '')[:n_lookahead] == feat_seq.replace('.', '')[i:i + n_lookahead])):
            #         # errata.append(i)
            #             print(i)
            #             seq_aligned += '.'
            #     # elif (not matched and
            #     #      (str(subseq).replace('.', '')[:n_lookahead] == feat_seq.replace('.', '')[i:i + n_lookahead])):
            #     #         seq_aligned += 
            #     else:
            #         seq_aligned += subseq[0]
            #         subseq = subseq[1:]
            # else:
            #     # if not subseq_finished and not reversed:
            #     #     seq_aligned += nt
            #     #     subseq_finished = True
            #     # else:
            #     seq_aligned += '.'
        seq_aligned = seq_aligned[::(-1 if reverse else 1)]
        # if not reverse:
        #     seq_aligned = seq_aligned[:-1]
        return seq_aligned

    def create_regex(self, seq : Sequence, spacers : bool = False, 
            fuzzy : bool = False, flank_start : bool = False,
            flank_end : bool = False, partial_seq : bool = False) -> str:

        regex = str(seq.expand().replace(' ', '').replace('*', ''))
        
        # n = 17
        n1 = 3
        n2 = 5
        # regex = '\.*'.join(list(regex)).replace('**', '*')
        pre_subs = errs = subs = ins = dels = 0

        if fuzzy:
            pre_subs = 0
            errs = 3
            subs = 3
            ins = 1
            dels = 1
        if len(seq) < (n1 + n2):
            return self._format_substring(regex)
        if spacers:
            regex = '({}){{s<={}}}({}){{e<={},s<={},i<={},d<={}}}({})'.format(
                '\.*'.join(regex[:n1].replace('.', '')), 
                pre_subs,
                '\.*'.join(regex[n1:-n2].replace('.', '')),
                errs, subs, ins, dels,
                '\.*'.join(regex[-n2:].replace('.', '')) + '.*?')
        else:
            regex = '({}){{s<={}}}({}){{e<={},s<={},i<={},d<={}}}({})'.format(
                self._format_substring(regex[:n1]),
                pre_subs,
                self._format_substring(regex[n1:-n2]),
                errs, subs, ins, dels,
                self._format_substring(regex[-n2:]))
        if flank_start:
            regex = regex + '$'
        elif flank_end:
            regex = '^' + regex
        # elif not partial_seq:
        #     regex = '^' + regex + '$'
        return regex

    def _format_substring(self, regex : str) -> str:
        regex = regex.replace('.', '.??').replace('X', '.')
        regex = regex.replace('\.\.+$', '')
        regex = re.sub('^\.\.+', '', regex)
        regex = re.sub('(\.\?)(\.\?)(\.\?)+', '.*?', regex)
        regex = re.sub('\.[*?]$', '', regex)
        regex = re.sub('[*.?]*$', '', regex)
        return regex

    def serialize(self):
        # results = {'snps' : self.snps}
        snps = {feat_name : str(snp) for feat_name, snp in self.snps.items()}
        return {'snps' : snps,
                'feat_names' : self.feat_names,
                'seq_features' : self.seq_features}

    def __add__(self, seq_search2 : any) -> any:
        if ('seqs' in self.seq_features) and ('seqs' in seq_search2.seq_features):
            feats1 = self.seq_features['seqs'].feats
            feats2 = seq_search2.seq_features['seqs'].feats
            for feat_name2 in feats2:
                if feat_name2 in feats1:
                    try:
                        feats1[feat_name2] = feats1[feat_name2] + feats2[feat_name2]
                    except Exception as e:
                        raise Exception(feat_name2, str(e))
                else:
                    feats1[feat_name2] = feats2[feat_name2]
        self.seq_features['seqs'].feats = feats1
        return SeqSearch(self.seq, self.seq_features_ref, seq_features=self.seq_features)
