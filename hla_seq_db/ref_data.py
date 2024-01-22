#
# Copyright (c) 2024 NMDP.
#
# This file is part of ExPAT 
# (see https://github.com/nmdp-bioinformatics/hla-db).
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
import os
from typing import Dict, List, Union
from .allotype import Allotype
from .sequence import Sequence
from .util import get_two_field_allele, calc_hla_class
from .ref_data_locus import RefDataLocus
from hla_seq_db.features import Features
from pyard.ard import ARD
import pandas as pd
import time

class RefData(object):
    def __init__(self, loci : List[str] = ['A', 'C', 'B', 'DRB1', 'DQB1', 'DPB1'],
                       db_version : str = '3460', ard : ARD = None, verbose :  bool = False):
        self.db_version = db_version.replace('.', '')
        self.ard = ard
        self.verbose = verbose
        self.loci = loci
        self.alleles = {}
        self.proteins_nt = {}
        self.proteins_aa = {}
        self.directory = os.path.dirname(__file__)
        self.ser_antigens = {}
        self.g_groups = self._load_groups('g')
        self.p_groups = self._load_groups('p')
        self.ser_antigens = self._load_serology()
        self.tce = None
        if 'DPB1' in loci:
            self.tce = self._get_tce_assignments()
        self.ciwd = self._download_ciwd()
        self._load_hla()
        self.cons_seq_whole_locus = {}

    def get_features(self, allotype : Allotype, feats: Union[str, List[str]] = None,
            field_level : int = 2) -> Features:
        if isinstance(allotype, str):
            allotype = Allotype(allotype, ref_data=self, ard=self.ard)
        return self.alleles[allotype.locus].get_features(allotype, feats, field_level=field_level)
    
    def _load_hla(self) -> None:
        """
        Loads gene_feature sequence data for HLA loci alleles and hi-res alleles/proteins.
        """
        motifs = {'DPB1' : 
                        {'exon_3_motif' : ('exon_3', [9, 16, 43, 80, 228, 236, 265]),
                            'rs9277534' : ('utr3', [791]),
                                 'ctcf' : ('intron_2', [1954, 1971, 1983]),
                                  'str' : ('intron_2', range(3893, 3960))}}
        if self.db_version >= '3500':
            motifs['DPB1']['exon_3_motif'] = ('exon_3', [9, 16, 43, 81, 229, 237, 266])
            motifs['DPB1']['str'] = ('intron_2', range(3894, 3961))
        for locus in self.loci:
            # self.alleles[locus] = self.load_seqs('allele', locus, db_type = 'gen')
            df = self.load_seqs('allele', locus, db_type = 'gen')
            self.alleles[locus] = RefDataLocus(df, locus, g_groups=self.g_groups[locus],
                proteins=self.load_seqs('prot', locus, db_type = 'prot'),
                tce=self.tce, ciwd=self.ciwd, ard=self.ard,
                motifs=motifs[locus] if locus in motifs else None)
            # self.proteins_nt[locus] = self.load_seqs('allele_hi_res', locus,
            #          df = self.alleles[locus])
            # self.proteins_aa[locus] = self.load_seqs('prot', locus, db_type = 'prot')
    
    # def _parse_gene_feats(self, df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    #     for col in df.columns:
    #         df[col].apply(lambda seq : Sequence(seq, name=col))
    #     return df

    def get_filepath(self, 
            locus : str = None, db_type : str = None,
            label : str = None) -> str:
        if locus and db_type:
            attribs = [locus, db_type, self.db_version]
            return '{}/data/alignments_{}.csv.gz'.format(self.directory,
                '_'.join([attrib for attrib in attribs if attrib]))
        elif label:
            return '{}/data/{}_{}.csv.gz'.format(self.directory,
                        label, self.db_version)

    def load_seqs(self, mode : str, locus : str, 
              db_type : str = None, df : pd.core.frame.DataFrame = None) -> pd.core.frame.DataFrame:
        """
        Loads sequences from IMGTHLA GitHub repository
        """
        db_version = self.db_version
        filepath = self.get_filepath(locus=locus, db_type=db_type)
        if os.path.exists(filepath):
            if self.verbose:
                print('Loaded', filepath)
            return pd.read_csv(filepath, index_col=0)
        if mode == 'allele':
            df_result = self._load_allele_seqs(locus)
        elif mode == 'allele_hi_res':
            df_result = self._collapse_to_prots(locus)
        elif mode == 'prot':
            df_result = self._create_df(self._load_seqs(locus, mode), locus,
                db_type=mode)
        else:
            raise Exception('Please input correct mode as either "allele" or "allele_hi_res"')
        if db_version != self.db_version:
            filepath = self.get_filepath(locus=locus, db_type=db_type)
        df_result.to_csv(filepath)
        return df_result
    
    def _collapse_to_prots(self, locus : str) -> pd.core.frame.DataFrame:
        filepath = 'hla_seq_db/data/proteins_{}_{}.csv'.format(locus, self.db_version)
        prot_con_seqs = {}
        for allele in self.alleles[locus].index:
            prot = get_two_field_allele(allele)
            if prot not in prot_con_seqs:
                prot_con_seqs[prot] = self.get_cons_seq(prot).split('|')
        return self._create_df(prot_con_seqs, locus)

    def load_allele_seqs(self, locus : str) -> pd.core.frame.DataFrame:
        df_gen = self._create_df(self._load_seqs(locus, 'gen'), locus, db_type='gen')
        df_nuc = self._create_df(self._load_seqs(locus, 'nuc'), locus, db_type='nuc')
        df_merged = df_nuc.join(df_gen, how="left", lsuffix='', rsuffix='_gen')
        return df_nuc

    def _load_allele_seqs(self, locus : str) -> pd.core.frame.DataFrame:
        df_gen = self._create_df(self._load_seqs(locus, 'gen'), locus, db_type='gen')
        df_nuc = self._create_df(self._load_seqs(locus, 'nuc'), locus, db_type='nuc')
        df_merged = df_nuc.join(df_gen, how="left", lsuffix='', rsuffix='_gen')
        # Fill in NA non-exon gene feats
        for col in df_merged.columns:
            seqs = df_merged[col]
            if len(seqs[seqs.isnull()]):
                seq = seqs.dropna().iloc[0]
                seq_empty = re.sub('[A-Z]', '*', seq)
                seqs = seqs.fillna(seq_empty)
                df_merged[col] = seqs
        df_merged = df_merged[df_gen.columns]
        return df_merged

    def _create_df(self, seqs : Dict[str, str], locus : str, db_type : str = 'gen') -> pd.core.frame.DataFrame:
        # print(locus, db_type, self._calc_features(locus, db_type))
        # print([(k, v) for (k, v) in seqs.items()][:3])
        return pd.DataFrame.from_dict(seqs, orient='index',
                        columns=self._calc_features(locus, db_type))
    
    def _calc_features(self, locus : str, db_type : str = 'gen') -> List[str]:
        """
        Calculates a list of feature names HLA loci.
        """
        if db_type in ['gen', 'nuc']:
            if locus == 'B':
                n_exons = 7
            elif locus in ['DRB1', 'DQB1']:
                n_exons = 6
            elif locus in ['DQA1']:
                n_exons = 4
            else:
                hla_class = calc_hla_class(locus)
                if hla_class == 1:
                    n_exons = 8
                elif hla_class == 2:
                    n_exons = 5
            feature_names = (['utr5'] + 
                            [feature + '_' + str(i)
                                for i in range(1, n_exons + 1) for feature in ['exon', 'intron'] 
                                if not (i == n_exons and feature == 'intron')] +
                            ['utr3'])
            if db_type=='nuc':
                feature_names = [feat for feat in feature_names if 'exon' in feat]
        elif db_type == 'prot':
            feature_names = ['aa']
        return feature_names

    def _load_seqs(self, locus : str, db_type : str) -> Dict[str, str]:
        """
        Loads sequences from IMGTHLA GitHub repository
        """
        if self.verbose:
            print('Loading seqs', locus, db_type)
        url = ("https://raw.githubusercontent.com"
                    "/ANHIG/IMGTHLA/%s/alignments/"
                    "%s_%s.txt" % (self.db_version, 
                        locus if not (db_type in ["prot", "nuc"] and locus == 'DRB1') else locus.replace('1', ''),
                        db_type))
        try:
            file = urlopen(url)
        except Exception as e:
            raise Exception('%s %s' % (url, str(e)))
        headers = {}
        pos = None
        alleles = {}
        ref_allele = None
        split_index = None
        for i, line in enumerate(file):
            decoded_line = line.decode("utf-8")
            if decoded_line[0] == '#':
                header_match = re.match("# ([a-z]+): (.+)\n", decoded_line)
                header, value = header_match.group(1,2)
                headers[header] = value
                if (self.db_version == 'Latest') and (header == 'version'):
                    self.db_version = value.replace('IPD-IMGT/HLA ', '').replace('.', '')
            else:
                if not split_index and ('|' in decoded_line) and (decoded_line.replace('|', '').strip() == ''):
                    split_index = decoded_line.index('|')
                else:
                    row = decoded_line.split()
                    if len(row) > 1:
                        if db_type == "gen" and not pos and re.match('\-?\d+', row[1]):
                            db_label, pos = row
                        elif split_index:
                            allele, row = decoded_line[:split_index].strip(), decoded_line[split_index:]
                            sequence = ''.join(row.split())
                            # allele, sequence = row[0], ''.join(row[1:])
                            if locus + '*' not in allele:
                                continue
                            if not ref_allele:
                                ref_allele = allele
                            if allele not in alleles:
                                alleles[allele] = ""
                            alleles[allele] += sequence
        ref_seq = alleles[ref_allele]
        for allele, seq_aligned in alleles.items():
            seq = ""
            for i, nt in enumerate(seq_aligned):
                if nt == '-': nt = ref_seq[i]
                seq += nt
            alleles[allele] = seq.split('|') 
        # print([(k, len(v)) for (k, v) in alleles.items() if len(v) > 10][:10])
        return alleles

    def _get_tce_assignments(self) -> Dict[str, str]:
        """
        Obtains TCE assignments from IMGT.
        :return: Dictionary of alleles to TCE groups (3, 2, 1, 0).
        :rtype: Dict[str, str]
        """
        allele_header = 'allele'
        tce_header = 'v2'
        try:
            url = 'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{}/tce/dpb_tce.csv'.format(self.db_version)
            tce_df = pd.read_csv(url, comment='#', 
                names=['allele', 'protein', 'v1', 'v2', 'comments'], header=0)
            tce_df_prot = tce_df[['protein', tce_header]].drop_duplicates().rename(columns={'protein' : 'allele'})
            tce_df = pd.concat([tce_df, tce_df_prot])
            tce_df = tce_df[[allele_header, tce_header]]
            tce_df.replace({'\$' : '', 'a' : ''}, regex=True, inplace=True)
            tce_df = tce_df[~tce_df[tce_header].isnull()]
            tce_df.set_index(allele_header, inplace=True)
            tce_map = tce_df.to_dict()[tce_header]
            return tce_map
        except Exception as e:
            if self.verbose:
                print('Error when loading TCE:', e)
            return None

    def _download_ciwd(self) -> pd.DataFrame:
        """
        Downloads the CIWD 3.0.0 table from the NMDP Bioinformatics repository
        as a Pandas DataFrame.
        :return: CIWD table as Pandas Dataframe with 'allele' and 'CIWD' status columns.
        :rtype: pd.DataFrame
        """
        url = 'https://raw.githubusercontent.com/nmdp-bioinformatics/CIWD/master/ciwd.3.0.0.csv'
        ciwd_df = pd.read_csv(url).add_prefix('CIWD_')\
                    .rename(columns={"CIWD_ALLELE" : "allele"})\
                    .set_index('allele')
        return ciwd_df

    # def _load_csv(self) -> pd.core.frame.DataFrame:
        

    def _load_groups(self, group_type : str) -> Dict[str, List[str]]:
        """
        Obtains the P group allele lists from IMGT and
        returns them in a list of lists for available parsing.
        """
        try:
            file = urlopen("https://raw.githubusercontent.com/ANHIG/IMGTHLA/{}/wmda/hla_nom_{}.txt"\
                .format(self.db_version, group_type))
        except Exception as e:
            # try:
            #     file = pd.read_csv(self.get_filepath(label=group_type))
            # except Exception as e:
            #     print(self.get_filepath(label=group_type + '_groups'))
            #     raise e
            raise e
        groups = {}
        for i, line in enumerate(file):
            decoded_line = line.decode("utf-8").replace('\n', '')
            if decoded_line[0] != '#':
                locus, alleles, group = decoded_line.split(';')
                locus = locus.replace('*', '')
                if locus not in groups:
                    groups[locus] = {}
                if group and locus in self.loci:
                    groups[locus][locus + '*' + group] = [locus + '*' + allele for allele in alleles.split('/')]
        return groups

    def _load_serology(self) -> Dict[str, List[str]]:
        """
        Obtains the different serology typings (keys) and their lists of possible alleles (values).
        """
        file = urlopen("https://raw.githubusercontent.com/ANHIG/IMGTHLA/{}/wmda/rel_dna_ser.txt"\
            .format(self.db_version))
        ser_typings = {}
        for i, line in enumerate(file):
            decoded_line = line.decode("utf-8").replace('\n', '')
            if decoded_line[0] != '#':
                # if self.locus + '*' in decoded_line:
                locus, allele, ser_unambig, ser_pos, ser_assum, ser_excep = decoded_line.split(';')
                sers = [ser for ser in [ser_unambig, ser_pos, ser_assum] if ser]
                ser = locus[:-1] + sers[0]
                if len(sers) > 1:
                    raise Exception('More serology assignments than expected for this allele', locus + allele)
                label = ''
                if ser_unambig:
                    label = 'unambiguous'
                elif ser_pos:
                    label = 'possible'
                elif ser_assum:
                    label = 'assumed'
                if not ser_excep and label == 'unambiguous':
                    if ser not in ser_typings:
                        ser_typings[ser] = []
                    ser_typings[ser].append(locus + allele)
        return ser_typings

    def get_map(self, locus : str,
            group_col : str, feature : str,
            unk_val : str = 'unknown') -> Dict[str, str]:
        return self.alleles[locus].df[[group_col, feature]]\
            .dropna().groupby(group_col).agg(lambda x: '/'.join(set([el for el in x if el != unk_val]))
                                                if set(x) != {unk_val} else unk_val).to_dict()[feature]


    def get_seq(self, allele : str, gene_feat : str = None) -> str:
        allele = Allotype(allele, ref_data=self)
        locus = allele.locus
        alleles = self.alleles[locus]
        if allele.typing in alleles:
            gene_feats = alleles.loc[allele.typing]
            if gene_feat:
                return gene_feats[gene_feat]
            return '|'.join(gene_feats)
        else:
            return self.get_cons_seq(allele.typing, gene_feat=gene_feat)
    
    def get_cons_seq(self, typing : str = None, gene_feat : str = None,
            df : pd.core.frame.DataFrame = None,
            alleles : List[str] = None, spacers : bool = True,
            most_freq_nt : bool = False, cds_only : bool = False) -> str:
        """
        Get consensus sequence of the provided HLA typing.
        Any ambiguity is denoted with an 'X'. You may also
        give a locus as a parameter i.e., 'B' to obtain 
        the consensus sequence for the whole locus.
        """
        # TODO: Refactor this
        whole_locus = False
        if isinstance(df, pd.core.frame.DataFrame) and not df.empty:
            df = df
        elif typing in self.loci and not alleles:
            df = self.alleles[typing]
        else:
            if typing in self.loci:
                locus = typing
                if locus in self.cons_seq_whole_locus:
                    return self.cons_seq_whole_locus[locus]
                whole_locus = True
            elif typing:
                if not isinstance(typing, Allotype):
                    allele = Allotype(typing, ref_data=self)
                else:
                    allele = typing
                locus = allele.locus
            if not alleles:
                alleles = [str(allele) for allele in allele.get_potential_alleles()]
            df = self.alleles[locus]
            df = df[df.index.isin(alleles)]
        seqs_feat_out = []
        for col in df.columns:
            if not ('utr' in col or 'exon' in col or 'intron' in col or 'str' in col):
                continue
            if cds_only and 'exon' not in col:
                continue
            if gene_feat and gene_feat != col:
                continue
            seqs_feat = df[col]
            seq_feat_out = ""
            seqs_feat = seqs_feat.dropna()
            if seqs_feat.empty:
                break
            if col == 'expression_experimental':
                continue #TODO: Remove this
            for i in range(len(seqs_feat.iloc[0])):
                nts = []
                nt = None
                for seq in seqs_feat:
                    if isinstance(seq, str):
                        if not most_freq_nt:
                            if seq[i] != '*':
                                nts.append(seq[i])
                            # nts.add(seq[i])
                            # if len(set(nts)) >= 2:
                            #     nt = 'X'
                            #     break
                        else:
                            if seq[i] != '*':
                                nts.append(seq[i])
                if not nt:
                    if not len(nts):
                        nt = '*'
                    else:
                        if most_freq_nt:
                            nt = max(set(nts), key = nts.count)
                        else:
                            nt_counts = pd.Series(nts).value_counts(normalize=True).to_dict()
                            # if col == 'exon_2':
                            #     print(nt_counts)
                            if (len(nt_counts) == 1) or list(nt_counts.values())[0] >= .99:
                                nt = list(nt_counts.keys())[0]
                            else:
                                nt = 'X'
                            # if len(set(nts)) == 1:
                            #     nt = nts[0]
                            # else:
                            #     nt = 'X'
                seq_feat_out += nt
            seqs_feat_out.append(seq_feat_out)
        seq_joined = '|'.join(seqs_feat_out)
        if not spacers:
            seq_joined = seq_joined.replace('.', '')
        if whole_locus:
            self.cons_seq_whole_locus[locus] = seq_joined
        return seq_joined
    
    def __repr__(self):
        return str(list(self.__dict__.keys()))