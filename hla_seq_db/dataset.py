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
import pandas as pd
import numpy as np
from typing import List, Tuple, Dict, Union
from os.path import exists
import re
import time
import signal

class Dataset(object):
    def __init__(self, filepath : str = None,
                headers_subj : List[str] = ['R', 'D', 'RECIP', 'DONOR'],
                headers_ids : List[str] = ['ID', 'RID', 'DID'],
                headers_locus : List[str] = ['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPB1'],
                headers_allele : List[str] = ['TYP', 'ALLOTYPE', '1', '2'],
                headers_num : List[str] = ['1', '2'],
                headers_gene_feats : List[str] = ['UTR5', 'EXON_1', 'INTRON_1', 'EXON_2', 'INTRON_2',
                                                 'EXON_3', 'INTRON_3', 'EXON_4', 'INTRON_4',
                                                 'EXON_5', 'INTRON_5', 'EXON_6', 'INTRON_6', 
                                                 'EXON_7', 'INTRON_7', 'EXON_8', 'UTR3'],
                headers_gene_feat_label : str = 'GENE_FEAT',
                headers_ref_snp_label : str = 'REF_SNPS',
                header_glstring : str = 'GL_STRING',
                headers_seq : List[str] = ['SEQ'],
                separator : str = '_',
                stack_loci : bool = False,
                hla_db = None,
                loci : List[str] = ['A', 'C', 'B', 'DRB1', 'DQB1', 'DPB1']):
        self.hla_db = hla_db
        self.loci = loci
        self.stack_loci = stack_loci
        self.filepath = filepath
        self.cols_stacked = None
        self.headers_subj = headers_subj
        self.headers_ids = headers_ids
        self.headers_locus = headers_locus
        self.headers_allele = headers_allele
        self.header_glstring = header_glstring
        self.headers_num = headers_num
        self.headers_gene_feats = headers_gene_feats
        self.headers_gene_feat_label = headers_gene_feat_label
        self.headers_ref_snp_label = headers_ref_snp_label
        self.headers_seq = headers_seq
        self.header_subj_index = 0
        self.header_locus_index = 1
        self.header_allele_index = 2
        self.headers_genotypes = {}
        self.separator = separator
        self.headers_metadata = None
        self.df = self._parse_dataframe(filepath)
        self.columns_orig = self.df.columns
        self.headers_data = None
        self._detect_headers()
        self.headers = {'subjects' : self.headers_subj,
                        'ids' : self.headers_ids,
                        'loci' : self.headers_locus,
                        'allotype' : self.headers_allele,
                        'gene features' : self.headers_gene_feats}
        self.paired = self._is_paired()
        # self.headers = self._parse_headers()
    
    def _is_paired(self) -> bool:
        return (len(self.headers_subj) >= 2)

    def _parse_dataframe(self, filepath : str) -> pd.core.frame.DataFrame:
        # ext = filepath.split('.')[1]
        dtype = {'sample-center-code' : 'string'}
        dtype = 'str'
        if isinstance(filepath, pd.core.frame.DataFrame):
            df = filepath
            filepath = None
        elif not filepath:
            self.df = pd.core.frame.DataFrame()
            headers = self.get_headers(headers_subj=self.headers_subj,
                    headers_locus=self.headers_locus,
                    headers_allele=self.headers_allele,
                    headers_type=[],
                    headers_gene_feats=self.headers_gene_feats)
            return pd.core.frame.DataFrame({header : [''] for header in headers})
        else:
            try:
                df = pd.read_csv(filepath, dtype=dtype)
            except:
                df = pd.read_excel(filepath, dtype=dtype)
        first_val = str(df.iloc[0][0])
        if first_val and (first_val[0] == '=') and (first_val[-1] == '"'):
            df = df.replace('(^=")|("$)', '', regex=True) # Trim '="example"' characters to 'example'.
        tce_cols = [col for col in df.columns if 'TCE' in col]
        df[tce_cols] = df[tce_cols].astype(str)
        # Rename duplicate columns
        suffix_def = '.1'
        remapper = {}
        for col in df.columns:
            col_bare = col.replace(suffix_def, '')
            if suffix_def in col and col_bare in df.columns:
                remapper.update({col_bare : col_bare + '_1',
                                col : col_bare + '_2'})
        df = df.rename(remapper, axis=1)

        if self.stack_loci:
            df = self._stack_df(df)
        return df             

    def _stack_df(self, df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
        col_stacked = [col for col in df.columns if re.search('_\d$', col)]
        cols_stacked = []
        for col in df.columns:
            if col not in col_stacked:
                df.rename({col : col + '_1'}, axis=1, inplace=True)
                df[col + '_2'] = df[col + '_1']
            else:
                cols_stacked.append(col.split('_')[0])
        self.cols_stacked = list(set(cols_stacked))
        cols = [tuple(col.split('_')) for col in df.columns]

        df.columns = pd.MultiIndex.from_tuples(cols)
        df = df.stack()
        col_order = ['Age', 'Ethnicity', 'GRID/ID', 'DRS %', 'Avail Date', 
                    'CCR5', 'CMV', 'Contact Type', 'Country',
                    'Pr (10/10) = %', 'Pr (9/10) = %', 'Pr (8/10) = %', 'Pr (7/10) = %', 'Pr (6/10) = %', 'Pr (5/10) = %', 
                    'Pr (8/8) = %', 'Pr (7/8) = %', 'Pr (6/8) = %', 'Pr (5/8) = %', 'Pr (4/8) = %', 'Pr (3/8) = %',
                    'A', 'A Match Grade', 'A%', 'ABO', 
                    'B', 'B Leader Match Status', 'B Match Grade', 'B%', 
                    'C', 'C Match Grade', 'C%', 'Antibody Indicator',
                    'DRB1', 'DRB1 Match Grade', 'DRB1%', 
                    'DRB3', 'DRB4', 'DRB5',
                    'DQB1', 'DQB1 Match Grade', 'DQB1%', 
                    'DPB1', 'DPB1 Match Grade', 'DPB1 TCE', 'TCE Grp',
                    'DQA1', 'DQA1 Match Grade',
                    'DPA1', 'DPA1 Match Grade',
                    'ION/DC ID', 'ION/DC Name', 'KIR CenB',
                    'Last Contact Date', 'Likely DRB3/4/5 Mismatches', 'List',
                    'Match Category (out of 10)', 'Match Category (out of 8)', 'Num Preg',
                    'Population', 'Prev. Donations', 'Previous CT',
                    'Previously released', 'Race', 'Reg Date', 'Registry Donor ID',
                    'Repository Sample', 'RhD Type', 'Sex', 'Status', 
                    'Weight (kg)']
        if all([col in col_order for col in df.columns]):
            df = df[col_order]
        return df

    def _detect_headers(self) -> None:
        """Filters out headers that are not present in file"""
        self.df.columns = [col.upper() for col in self.df.columns]
        headers_subj = set()
        headers_locus = set()
        headers_allele = set()
        headers_num = set()
        headers_gene_feat = set()
        headers_gene_feat_label = set()
        headers_ref_snp_label = set()
        headers_data = []
        for col in self.df.columns:
            header_num = col[-1] 
            if header_num in self.headers_num:
                col = col[:-1]
                headers_num.add(header_num)
            else:
                header_num = None
            subheaders = col.split(self.separator)
            if len(subheaders) == 3:
                headers_subj.add(subheaders[self.header_subj_index])
                headers_locus.add(subheaders[self.header_locus_index])
                headers_allele.add(subheaders[self.header_allele_index])
                if header_num:
                    col += header_num
                headers_data.append(col)
            elif len(subheaders) > 3:
                indices_joined = []
                for i in range(len(subheaders) - 1):
                    subheaders_adj = subheaders[i] + '_' + subheaders[i + 1]
                    if subheaders_adj in [self.headers_gene_feats, self.headers_ref_snp_label, self.headers_gene_feat_label]:
                        indices_joined.append(i)
                subheaders_upd = []
                for i in range(len(subheaders)):
                    if (i in indices_joined):
                        subheaders_upd.append(subheaders[i] + '_' + subheaders[i + 1])
                    elif (i - 1) not in indices_joined:
                        subheaders_upd.append(subheaders[i])
                subheaders = subheaders_upd
                if len(subheaders) == 5:
                    headers_subj.add(subheaders[0])
                    headers_locus.add(subheaders[1])
                    headers_gene_feat.add(subheaders[2])
                    if subheaders[3] in self.headers_ref_snp_label:
                        headers_ref_snp_label.add(subheaders[3])
                    elif subheaders[3] in self.headers_gene_feat_label:
                        headers_gene_feat_label.add(subheaders[3])
                if header_num:
                    col += header_num
                headers_data.append(col)
        headers_subj = [header for header in self.headers_subj if header in headers_subj]
        self.headers_locus = [header for header in self.headers_locus if header in headers_locus]
        self.headers_num = [header for header in self.headers_num if header in headers_num]
        headers_allele = [header for header in self.headers_allele if header in headers_allele]
        if not len(headers_allele):
            headers_allele = [header for header in self.headers_allele if header in self.df.columns]
        if not len(headers_subj):
            headers_subj = [header for header in self.df.columns if header in self.headers_subj]
        self.headers_allele = headers_allele
        self.headers_seq = [header for header in self.headers_seq if header in self.df.columns]
        if self.header_glstring not in self.df.columns:
            self.header_glstring = None
        self.headers_subj = headers_subj
        self.headers_data = headers_data
        self.headers_ids = [header for header in self.df.columns if header in self.headers_ids]
        self.headers_metadata = [col for col in self.df.columns if col not in headers_data]

    def get_dataframe(self, headers_subj : List[str] = None,
                          headers_locus : List[str] = None,
                          headers_allele : List[str] = None,
                          headers_gene_feats : List[str] = None,
                          headers_type : List[str] = ['REF_SNPS', 'GENE_FEAT'],
                          headers_num : List[str] = ['1', '2'],
                          metadata_only : bool = False,
                          metadata : bool = False,
                          suffix : str = '') -> pd.core.frame.DataFrame:
        headers = []
        if not metadata_only:
            headers = self.get_headers(headers_subj=headers_subj,
                headers_locus=headers_locus,
                headers_allele=headers_allele,
                headers_gene_feats=headers_gene_feats, metadata_only=metadata_only,
                metadata=metadata, headers_num=headers_num,
                headers_type=headers_type, suffix=suffix)
        if not self.df.empty:
            headers = [header for header in headers if header in self.df and
                                                     ~self.df[header].isnull().all()]
        return self.df[headers]

    def get_headers(self, headers_subj : List[str] = None,
                          headers_locus : List[str] = None,
                          headers_allele : List[str] = None,
                          headers_gene_feats : List[str] = None,
                          headers_type : List[str] = ['REF_SNPS', 'GENE_FEAT', 'SNPS'],
                          headers_num : List[str] = ['1', '2'],
                          metadata_only : bool = False,
                          metadata : bool = False,
                          suffix : str = '', new_headers : bool = False) -> List[str]:
        headers = []
        if headers_subj == None: headers_subj = self.headers_subj
        if headers_locus == None: headers_locus = self.headers_locus
        if headers_allele == None: headers_allele = self.headers_allele
        if headers_gene_feats == None: headers_gene_feats = self.headers_gene_feats
        if isinstance(headers_allele, list):
            headers_allele = sorted(headers_allele)
        if isinstance(headers_subj, str):
            headers_subj = [headers_subj]
        if isinstance(headers_locus, str):
            headers_locus = [headers_locus]
        if isinstance(headers_allele, str):
            headers_allele = [headers_allele]
        if isinstance(headers_gene_feats, str):
            headers_gene_feats = [headers_gene_feats]
        if isinstance(headers_type, str):
            headers_type = [headers_type]
        # if headers_gene_feat:
        for i, header_subj in enumerate(headers_subj):
            for header_locus in headers_locus:
                for header_num in headers_num:
                    for header_allele in headers_allele:
                    # if header_subj not in self.headers_subj: continue
                    # if header_locus not in self.headers_locus: continue
                    # if header_allele not in self.headers_allele: continue
                        row = [None] * 3
                        row[self.header_subj_index] = header_subj
                        row[self.header_locus_index] = header_locus
                        row[self.header_allele_index] = header_allele
                        headers.append(self.separator.join(row) + header_num)
                    for header_type in headers_type:
                        for header_gene_feat in headers_gene_feats:
                            if (header_type == 'SNPS') and (i == 0):
                                headers.append(self.separator.join([header_locus, header_gene_feat, 'SNPS']))
                                headers.append(self.separator.join([header_locus, header_gene_feat, 'NUM_SNPS_TOTAL']))
                                headers.append(self.separator.join([header_locus, header_gene_feat, 'NUM_SNPS_MISSENSE']))
                            elif 'GENE_FEAT' == header_type:
                                header = self.separator.join([header_subj, header_locus, header_gene_feat, 'GENE_FEAT', header_num])
                                header = header.replace('GENE_FEAT_', 'GENE_FEAT')
                                headers.append(header)
                            else:
                                header = self.separator.join([header_subj, header_locus, header_gene_feat, header_type])
                                headers.append(header)
        if 'SNPS' in headers_type:
            headers += ['NUM_SNPS_TOTAL', 'NUM_SNPS_MISSENSE']
        if suffix:
            headers = [header + self.separator + suffix for header in headers]
        if not new_headers and not self.df.empty:
            headers = [header for header in headers if header in self.df]
        if (metadata_only or metadata) and self.headers_metadata:
            headers = self.headers_metadata + headers
        return headers

    def get_genotypes(self, header_subj : str, locus : str, suffix : str = '') -> Tuple[List[str], pd.core.series.Series]:
        headers = self.get_headers(headers_locus=[locus], headers_subj=[header_subj], suffix=suffix)
        fields = self.df[headers].replace('^(\d):', '0\\1:', regex=True)
        df_alleles = locus + '*' + fields
        df_alleles = df_alleles.replace("^\s*(\w+)\*(\d+)+$", "\\1\\2", regex=True)

        # TEMPORARY:
        def assign_loci(el : str) -> str:
            if isinstance(el, str):
                alleles = el.split('/');
                first_allele = alleles[0]
                alleles_out = []
                if '*' in first_allele:
                    locus = alleles[0].split('*')[0]
                    for allele in alleles:
                        if '*' not in allele:
                            allele = locus + '*' + allele
                        alleles_out.append(allele)
                    return '/'.join(alleles_out)
            return el

        # df_out = df.copy()
        for col in df_alleles.columns:
            series = df_alleles[col]
            df_alleles[col] = series.apply(assign_loci)

        df_alleles = df_alleles.replace(locus + '*', np.NaN)
        for col in df_alleles.columns:
            df_alleles[col] = df_alleles[col].str.replace(locus + '*' + locus, locus, regex=False)

        return headers, df_alleles.apply(lambda alleles: '+'.join(([allele 
                                            for allele in alleles
                                            if isinstance(allele, str)] * 2)[:2]), axis=1)

    def annotate(self,
            seqs : pd.core.frame.DataFrame = None,
            annotations : List[str] = [], 
            nrows : int = None,
            loci : List[str] = None,
            index : int = None,
            output : str = None,
            overwrite : bool = False,
            n_print_out : int = 50,
            verbose : bool = False) -> pd.core.frame.DataFrame:
        """
        Specific annotations can be provided in the output such as:
            'genotypes' : Columns of formatted genotypes such as 'R_DPB1_geno'.
            'tce_match' : Includes TCE matching in the 'tce_match' column.
            'tce' : TCE group assignments. An example column would be 'R_DPB1_TYP1_TCE'.
        """
        if output and exists(output):
            return pd.read_csv(output)
        annotations = [ann.upper() for ann in annotations]
        if isinstance(loci, str):
            loci = [loci]
        loci = loci or self.headers_locus
        if verbose:
            print('Annotating a file with {:,} rows'.format(len(self.df)))
        if self.paired:
            if isinstance(seqs, pd.core.frame.DataFrame):
                self.annotate_individual_matches_hi_res(seqs=seqs, nrows=nrows, n_print_out=n_print_out, verbose=verbose)
            else:
                if 'DQA1' in loci and 'DQB1' in loci:
                    self.annotate_individual_matches(loci=loci, nrows=nrows, n_print_out=n_print_out, verbose=verbose)
                else:
                    for locus in loci:
                        self.annotate_matches(locus, annotations=annotations, nrows=nrows, n_print_out=n_print_out, verbose=verbose)
        elif self.header_glstring:
            self.annotate_individuals(nrows=nrows, n_print_out=n_print_out, verbose=verbose)
        else:
            self.annotate_allotypes(index=index, annotations=annotations, nrows=nrows, n_print_out=n_print_out, verbose=verbose, overwrite=overwrite)
            #     print('annotating')
        
                # if verbose:
                #     # if i % 10 == 0:
                #     print(i, geno_pair, len(cache))
        if 'genotypes' not in annotations:
            self.df = self.df[[col for col in self.df.columns]] # if col not in headers_subj_locus
        if output:
            try:
                self.df.to_csv(output, index=False)
            except:
                self.df.to_excel(output, index=False)
        return self.df

    def annotate_allotypes(self, annotations : List[str] = [], nrows : int = None, index : int = None,
        annotate_seqs : bool = True,
        n_print_out : int = None, verbose : bool = False, overwrite : bool = False) -> pd.core.frame.DataFrame:

        def handler(signum, frame):
            raise Exception("Ran too long!")

        def flatten_dict(annotation : Union[str, Dict[str, any]], prefix : str = "") -> Dict[str, str]:
            if isinstance(annotation, str):
                return {prefix : annotation}
            elif isinstance(annotation, dict):
                ann_items = annotation.items()  
            elif isinstance(annotation, list):
                print(prefix, annotation)
                return None
            elif isinstance(annotation, bool):
                return {prefix : annotation}
            elif isinstance(annotation, int):
                return {prefix : annotation}
            else:
                ann_items = annotation.serialize().items()
            annotation_flat = {}
            for k, v in ann_items:
                if k not in ['alleles_hi_res', 'seq_diffs']:
                    if prefix:
                        prefix = (prefix + '_').replace('__', '_')
                    # print(annotation_flat, flatten_dict(v, prefix=k), type(v))
                    annotation_flat.update(flatten_dict(v, prefix=prefix + k))
            return annotation_flat

        n = 0
        gene_feat = None
        headers = ['ALLOTYPE']
        annotate_seqs = annotate_seqs and ('SEQ' in self.df)
        if annotate_seqs:
            headers += ['SEQ']
        df_unique = self.df[headers].drop_duplicates()
        if verbose:
            print('Parsing through {:,} unique seq-allotypes.'.format(len(df_unique)))
        start_time_abs = time.time()
        # if verbose:
        for i, row in df_unique.iterrows():
            if nrows == None or i < nrows:
                if (index == None) or (index == i):
                    if annotate_seqs:
                        seq = row['SEQ']
                    else:
                        seq = None
                    allotype = row['ALLOTYPE']

                    seq_locator = (self.df['ALLOTYPE'] == allotype)
                    if annotate_seqs:
                        seq_locator = (seq_locator & self.df['SEQ'])
                    rows = self.df.loc[seq_locator]
                    rows_ran = 'RAN' in rows and not rows.empty and (rows.iloc[0]['RAN'] == True)
                    if not overwrite and rows_ran:
                        continue
                    if 'gene_feat' in row:
                        gene_feat = row['gene_feat']
                    allotype_seq = (allotype, seq)
                    if (index != None) and verbose:
                        print(allotype_seq, len(seq), gene_feat)

                    signal.signal(signal.SIGALRM, handler)
                    signal.alarm(30)
                    try:
                    # if True:
                        if allotype and (not annotate_seqs or seq):
                            start_time = time.time()
                            # if allotype_seq in cache:
                            #     anno_allo = cache[allotype_seq]
                            # else:
                            anno_allo = self.hla_db.annotate_allotype(allotype, feats=annotations, seq=seq, verbose=verbose)
                            # gene_feat_names = '-'.join(anno_allo.gene_feats.keys())
                            anno_allo = anno_allo.serialize()
                            # print(anno_allo)
                            anno_allo = flatten_dict(anno_allo)
                                # cache[allotype_seq] = anno_allo
                            for k, v in anno_allo.items():
                                # if k not in ['gene_feats_features']:
                                    # print(k, v)
                                    # self.df.at[i, k] = v
                                self.df.loc[seq_locator, k.upper()] = v
                            # self.df.at[i, 'gene_feat_names'] = gene_feat_names
                            # self.df.at[i, 'seq_len'] = len(seq)
                            # self.df.loc[seq_locator, 'GENE_FEAT_NAMES'] = gene_feat_names
                            if annotate_seqs:
                                self.df.loc[seq_locator, 'SEQ_LEN'] = len(seq)
                            self.df.loc[seq_locator, 'RAN'] = True
                            n += 1
                            if (n_print_out and i and i % n_print_out == 0):
                                dur = time.time() - start_time
                                dur_abs = time.time() - start_time_abs
                                spf = dur_abs / n
                                print('Ran through {} rows. {:.1f} secs have elapsed. {:.2f} secs per result.'\
                                    .format(n, dur_abs, spf))
                    except Exception as e:
                        print(i, allotype_seq, 'error', e)
                        self.df.loc[seq_locator, 'error'] = str(e)
                    signal.alarm(0)
        return self.df

    def annotate_individuals(self, annotations : List[str] = [], nrows : int = None, index : int = None,
            n_print_out : int = None, verbose : bool = False, overwrite : bool = False) -> pd.core.frame.DataFrame:
        
        hla_db = self.hla_db

        def handler(signum, frame):
            raise Exception("Ran too long!")

        def flatten_dict(annotation : Union[str, Dict[str, any]], prefix : str = "") -> Dict[str, str]:
            if isinstance(annotation, str):
                return {prefix : annotation}
            elif isinstance(annotation, dict):
                ann_items = annotation.items()  
            elif isinstance(annotation, list):
                print(prefix, annotation)
                return None
            elif isinstance(annotation, bool):
                return {prefix : annotation}
            elif isinstance(annotation, int):
                return {prefix : annotation}
            else:
                ann_items = annotation.serialize().items()
            annotation_flat = {}
            for k, v in ann_items:
                if k not in ['alleles_hi_res', 'seq_diffs']:
                    if prefix:
                        prefix = (prefix + '_').replace('__', '_')
                    # print(annotation_flat, flatten_dict(v, prefix=k), type(v))
                    annotation_flat.update(flatten_dict(v, prefix=prefix + k))
            return annotation_flat

        n = 0
        gl_strings = self.df[self.header_glstring].drop_duplicates()
        if verbose:
            print('Parsing through {:,} unique seq-allotypes.'.format(len(gl_strings)))
        start_time_abs = time.time()
        # if verbose:
        for i, gl_string in enumerate(gl_strings):
            if nrows == None or i < nrows:
                if (index == None) or (index == i):
                    locator = self.df[self.header_glstring] == gl_string
                    rows = self.df.loc[locator]
                    rows_ran = 'RAN' in rows and not rows.empty and (rows.iloc[0]['RAN'] == True)
                    if not overwrite and rows_ran:
                        continue

                    signal.signal(signal.SIGALRM, handler)
                    signal.alarm(30)
                    try:
                        start_time = time.time()
                        annotation = hla_db.annotate_individual(gl_string).annotation
                        annotation = flatten_dict(annotation)
                        for k, v in annotation.items():
                            self.df.loc[locator, k.upper()] = v
                        self.df.loc[locator, 'RAN'] = True
                        n += 1
                        if (n_print_out and i and i % n_print_out == 0):
                            dur = time.time() - start_time
                            dur_abs = time.time() - start_time_abs
                            spf = dur_abs / n
                            print('Ran through {} rows. {:.1f} secs have elapsed. {:.2f} secs per result.'\
                                .format(n, dur_abs, spf))
                    except Exception as e:
                        print(i, gl_string, str(e))
                        self.df.loc[locator, 'error'] = str(e)
                    signal.alarm(0)
        return self.df


    def annotate_matches(self, locus : str,
            headers_subj_locus : List[str] = [],
            annotations : List[str] = [], nrows : int = None, 
            prefix : str = '', suffix : str = '',
            n_print_out : int = 10, verbose : bool = False) -> pd.core.frame.DataFrame:

        def handler(signum, frame):
            raise Exception("Ran too long!")

        headers_genotypes = {}
        # Genotype assignment
        cache = {}
        start_time_abs = time.time()
        if not headers_subj_locus:
            for header_subj in self.headers_subj:
                if header_subj not in headers_genotypes:
                    headers_genotypes[header_subj] = []
                header_subj_locus = self.get_headers(headers_subj=header_subj, headers_locus=locus, 
                                headers_gene_feats=[], headers_num=[''], headers_type=[], 
                                headers_allele=['GENO'], new_headers=True)[0]
                headers_subj_locus.append(header_subj_locus)
                if header_subj_locus not in self.df:
                    header_geno, genos = self.get_genotypes(header_subj, locus)
                    self.df[header_subj_locus] = genos
                    headers_genotypes[header_subj].append(header_subj_locus)
        n = 0
        if verbose:
            print('Parsing through {:,} rows.'.format(len(self.df)))
        for i, row in self.df.iterrows():
            if nrows == None or i < nrows:
                geno_recip = row[headers_subj_locus[0]]
                geno_donor = row[headers_subj_locus[1]]
                geno_pair = (geno_recip, geno_donor)
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(60)
                # try:
                if True:
                    if geno_recip and geno_donor:
                        start_time = time.time()
                        if geno_pair in cache:
                            match = cache[geno_pair]
                        else:
                            match = hla_db.annotate_match(geno_recip, geno_donor, feats=annotations)
                            cache[geno_pair] = match
                        if len(annotations):
                            if 'EXPR_GENO' in annotations:
                                self.df.at[i, prefix + 'R_EXPR_GENO' + suffix] = ', '.join(match.genotype_patient.annotation['expression'])
                                self.df.at[i, prefix + 'D_EXPR_GENO' + suffix] = ', '.join(match.genotype_donor.annotation['expression'])
                            if 'HAP_GRADE' in annotations:
                                self.df.at[i, prefix + 'HAP_GRADE' + suffix] = match.grade
                            if locus == 'B' and 'b_leader' in annotations:
                                self.df.at[i, prefix + 'B_LEADER' + suffix] = match.annotation['b_leader']
                            if locus == 'DPB1':
                                if 'TCE_MATCH' in annotations:
                                    self.df.at[i, prefix + 'TCE_MATCH' + suffix] = match.annotation['tce_match']
                                if 'EXPR_MATCH' in annotations:
                                    self.df.at[i, prefix + 'EXPR_MATCH' + suffix] = match.annotation['expr_match']
                            # if ('diffs' in annotations) and (locus == 'B'):
                            #     self.df.at[i, locus + '_diffs'] = match.annotation['b_leader']
                            if 'TCE' in annotations:
                                tces_recip = match.genotype_patient.annotation['tce']
                                self.df.at[i, prefix + 'R_DPB1_TYP1_TCE' + suffix] = tces_recip[0]
                                self.df.at[i, prefix + 'R_DPB1_TYP2_TCE' + suffix] = tces_recip[1]
                                tces_donor = match.genotype_donor.annotation['tce']
                                self.df.at[i, prefix + 'D_DPB1_TYP1_TCE' + suffix] = tces_donor[0]
                                self.df.at[i, prefix + 'D_DPB1_TYP2_TCE' + suffix] = tces_donor[1]
                        n += 1
                        if (i and n_print_out and i % n_print_out == 0) and verbose:
                            dur = time.time() - start_time
                            dur_abs = time.time() - start_time_abs
                            spf = dur_abs / n
                            print('Ran through {} rows. {:.1f} secs have elapsed. {:.2f} secs per row.'\
                                .format(n, dur_abs, spf))
                # except Exception as e:
                #     print(i, geno_pair, 'error', e)
                #     self.df.at[i, 'error'] = str(e)
                signal.alarm(0)
        self.headers_genotypes = headers_genotypes
        return self.df

    def annotate_individual_matches(self, loci : List[str], index : int = None,
            nrows : int = None, 
            prefix : str = '', suffix : str = '', overwrite : bool = False,
            n_print_out : int = 10, verbose : bool = False) -> pd.core.frame.DataFrame:

        hla_db = self.hla_db
        def handler(signum, frame):
            raise Exception("Ran too long!")

        def flatten_dict(annotation : Union[str, Dict[str, any]], prefix : str = "") -> Dict[str, str]:
            if isinstance(annotation, str):
                return {prefix : annotation}
            elif isinstance(annotation, dict):
                ann_items = annotation.items()  
            elif isinstance(annotation, list):
                print(prefix, annotation)
                return None
            elif isinstance(annotation, bool):
                return {prefix : annotation}
            elif isinstance(annotation, int):
                return {prefix : annotation}
            else:
                ann_items = annotation.serialize().items()
            annotation_flat = {}
            for k, v in ann_items:
                if k not in ['alleles_hi_res', 'seq_diffs']:
                    if prefix:
                        prefix = (prefix + '_').replace('__', '_')
                    # print(annotation_flat, flatten_dict(v, prefix=k), type(v))
                    annotation_flat.update(flatten_dict(v, prefix=prefix + k))
            return annotation_flat

        n = 0
        # gl_strings = self.df[self.header_glstring].drop_duplicates()
        typing = self.get_dataframe(headers_locus=loci).drop_duplicates()
        if verbose:
            print('Parsing through {:,} unique rows.'.format(len(typing)))
        start_time_abs = time.time()
        # if verbose:
        for i, row in typing.iterrows():
            if nrows == None or i < nrows:
                if (index == None) or (index == i):
                    row_typing = row.to_dict()
                    query = ' & '.join(["({} == '{}')".format(locus, hla)
                                 for locus, hla in row_typing.items()])
                    rows = self.df.query(query)
                    locator = rows.index
                    
                    rows_ran = 'RAN' in rows and not rows.empty and (rows.iloc[0]['RAN'] == True)
                    if not overwrite and rows_ran:
                        continue

                    signal.signal(signal.SIGALRM, handler)
                    signal.alarm(30)
                    try:
                        start_time = time.time()
                        annotation = hla_db.annotate_individual_match(row).annotation
                        annotation = flatten_dict(annotation)
                        for k, v in annotation.items():
                            self.df.loc[locator, k.upper()] = v
                        self.df.loc[locator, 'RAN'] = True
                        n += 1
                        if (n_print_out and i and i % n_print_out == 0):
                            dur = time.time() - start_time
                            dur_abs = time.time() - start_time_abs
                            spf = dur_abs / n
                            print('Ran through {} rows. {:.1f} secs have elapsed. {:.2f} secs per result.'\
                                .format(n, dur_abs, spf))
                    except Exception as e:
                        print(i, row, str(e))
                        self.df.loc[locator, 'error'] = str(e)
                    signal.alarm(0)
        return self.df

    def annotate_individual_matches_hi_res(self, seqs : pd.core.frame.DataFrame = None,
            annotations : List[str] = [], nrows : int = None, loci : List[str] = ['A', 'C', 'B', 'DQB1', 'DRB1' 'DPB1'],
            n_print_out : int = None, verbose : bool = False) -> pd.core.frame.DataFrame:
        
        hla_db = self.hla_db
        def handler(signum, frame):
            raise Exception("Ran too long!")
        seqs = Dataset(seqs).df

        n = 0
        start_time_abs = time.time()
        for i, row in self.df.iterrows():
            if nrows == None or i < nrows:
                id1 = row[self.headers_ids[0]]
                id2 = row[self.headers_ids[1]]
                df_seqs = pd.concat([seqs[seqs['ID'] == id1], seqs[seqs['ID'] == id2]])
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(180)
                try:
                # if True:
                    if not df_seqs.empty:
                        start_time = time.time()
                        match = hla_db.annotate_individual_match(df_seqs, loci=loci, align_mismatches_only=False)

                        missense_mms = 0
                        all_mms = 0
                        for geno_match in match.geno_matches:
                            for allele_match in geno_match.allele_matches:
                                # for allotype in allele_match.allotypes:
                                #     if allotype.features_searched and allotype.feats_searched.snps:
                                for k, individual in enumerate(match.individuals):
                                    label = ['R_', 'D_'][k]
                                    for genotype in individual.genos:
                                        for j, allotype in enumerate(genotype.allotypes):
                                            gene_feats = allotype.feats_searched.seqs.feats.values()
                                            self.df.at[i, label + str(allotype.locus) + '_TYP' + str(j + 1)] = allotype.typing
                                            for gene_feat in gene_feats:
                                                self.df.at[i, label + str(allotype.locus) + '_' + gene_feat.name.upper() + '_GENE_FEAT' + str(j + 1)] = str(gene_feat.seq_two)
                                            if allotype.feats_searched.snps:
                                                for feat_name, snps in allotype.feats_searched.snps.items():
                                                    locus_feat = allele_match.locus + '_' + feat_name.upper()
                                                    self.df.at[i, label + locus_feat + '_REF_SNPS_' + str(j + 1)] = ', '.join([str(snp) for snp in snps])
                                if allele_match.gene_feat_snps:
                                    for feat_name, snps in allele_match.gene_feat_snps.items():
                                        locus_feat = allele_match.locus + '_' + feat_name.upper()
                                        missense_mms_local = len([snp for snp in snps if snp.missense])
                                        all_mms_local = len(snps)
                                        snps = ', '.join([str(snp) for snp in snps])
                                        print(i, locus_feat, snps)
                                        self.df.at[i, locus_feat + '_SNPs'] = snps
                                        self.df.at[i, locus_feat + '_NUM_SNPS_TOTAL'] = all_mms_local
                                        self.df.at[i, locus_feat + '_NUM_SNPS_MISSENSE'] = missense_mms_local
                                        missense_mms += missense_mms_local
                                        all_mms += all_mms_local
                        self.df.at[i, 'NUM_SNPS_TOTAL'] = all_mms
                        self.df.at[i, 'NUM_SNPS_MISSENSE'] = missense_mms
                        self.df.at[i, 'RAN'] = True
                        n += 1
                        if (i and n_print_out and i % n_print_out == 0) and verbose:
                            dur = time.time() - start_time
                            dur_abs = time.time() - start_time_abs
                            spf = dur_abs / n
                            print('Ran through {} rows. {:.1f} secs have elapsed. {:.2f} secs per row.'\
                                .format(n, dur_abs, spf))
                except Exception as e:
                    print(i, id1, id2, 'error', e)
                    self.df.at[i, 'error'] = str(e)
                    self.df.at[i, 'RAN'] = True
                signal.alarm(0)
        # self.headers_genotypes = headers_genotypes
        return self.df

    def format_counts(self, header : str, replacements : Dict[str, str] = None,
                 label : str = None, limit : int = None, indent : int = 0) -> None:
        series = self.df[header].astype(str)
        output = ((series.value_counts(normalize=False, dropna=False)).apply('{:,} '.format).str.replace(',', ',') +
                (series.value_counts(normalize=True, dropna=False) * 100).apply('({:.1f}%)'.format))
        output = output.reset_index()
        if replacements:
            output = output.replace({'index' : replacements})
        output = output.agg(lambda x : ' with '.join([y for y in x if isinstance(y, str)]), axis=1)
        output = output.str.replace('^0.0%', '<0.1%', regex=True)
        if label:
            indent += 1
        for n, line in enumerate(output):
            indent_str = ('\t' * indent)
            if label and (n == 0):
                print(('\t' * (indent - 1)) + '- ' + "{} ({:,})".format(label, len(series)))
            if (limit != None) and (n == limit):
                print(indent_str + ' ...')
                break
            print(indent_str + '- ' + line)

    def __repr__(self):
        return str(self.df)

    # def _parse_headers(self):
    #     subject_header_index = set()
    #     locus_header_index = set()
    #     allele_header_index = set()
    #     for col in self.df.columns:
    #         if self.separator in col:
    #             headers = col.split(self.separator)
    #             print(col)
    #             print(self.subject_headers)
    #             print(headers)
    #             subject_header_index.add(self.subject_headers.index(headers))
    #             locus_header_index.add(self.locus_headers.index(headers))
    #             allele_header_index.add(self.allele_headers.index(headers))
    #     print(subject_header_index,
    #             locus_header_index,
    #             allele_header_index)

    # def get_genotype(loci : List[str]):
        # for locus in loci:
# 