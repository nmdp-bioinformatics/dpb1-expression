#
# Copyright (c) 2024 NMDP.
#
# This file is part of 
# (see https://github.com/nmdp-bioinformatics/).
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
import xml.etree.ElementTree as ET
import os
from typing import List, Dict, Union
import pandas as pd
import time
from zipfile import ZipFile, is_zipfile
import copy
import sqlalchemy
from sqlalchemy_utils import database_exists
import sqlite3
import signal

class HmlParser(object):

    def __init__(self) -> None:
        self.name = "hml_database"
        self.db = self._create_db(self.name)
        # self.db = create_engine('sqlite:///{}.db'.format(self.name))
    
    def _create_db(self, name : str) -> sqlalchemy.engine.base.Engine:
        try:
            db = sqlalchemy.create_engine("sqlite:///{}.db".format(name))
        except Exception as e:
            raise e
        self.db = db
        self.name = name
        return db
    
    def _process_el(self, el : ET.Element) -> str:
        return el.tag.split('}')[1], el.attrib, (el.text and el.text.strip()), bool(len(el))

    def _prepend_attributes(self, tag : str, attrib : Dict[str, str]) -> Dict[str, str]:
        return {tag + '-' + k : v for k, v in attrib.items()}
    
    def recurse_element(self, el : ET.Element):
        tag, attrib, value, has_children = self._process_el(el)
        contents = self._prepend_attributes(tag, attrib)
        if has_children:
            for child in el:
                contents.update(self.recurse_element(child))
        elif attrib:
            return self._prepend_attributes(tag, attrib)
        else:
            return {tag : value}
        return contents

    def process_sequences(self, el : ET.Element, metadata : Dict = {}, filename : str = None) -> List[Dict[str, str]]:
        tag, attrib, value, has_children = self._process_el(el)
        sequences = []
        # headers = ["hml", "sample", "sample-id", "glstring", "filename", "reference-sequence-name", "consensus-sequence-block-description",
        #     "hml-project-name", "sbt-ngs-test-id", "sbt-ngs-test-id-source", "typing-date", "allele-assignment-allele-version",
        #     "sequence"
        # ]
        # if tag not in headers:
        #     return []
        if filename:
            metadata['filename'] = filename
        if attrib and not value:
            if tag == 'haploid':
                # if not ((attrib['locus'] in ['HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']) and (attrib['type'] == 'NNNN')):
                    # print(tag, attrib, value)
                #TODO: Figure out DPB1
                pass
            else:
                metadata.update(self._prepend_attributes(tag, attrib))
        else:
            metadata.update({tag : value})
        if tag == 'consensus-sequence':
            sequences = []
            seq_id_header = 'consensus-sequence-block-reference-sequence-id'
            ids_added = []
            for child in el:
                tag, attrib, value, has_children = self._process_el(child)
                if has_children:
                    consensus_block = self.recurse_element(child)
                    consensus_block.update(metadata)
                    if tag == 'consensus-sequence-block':
                        seq_id = consensus_block[seq_id_header]
                        sequences_new = []
                        for seq in sequences:
                            if (seq['reference-sequence-id'] == seq_id):
                                if seq_id in ids_added:
                                    sequences_new.append(seq)
                                    seq = copy.copy(seq)
                                seq.update(consensus_block)
                                # if filename:
                                #     seq['filename'] = filename
                                seq.pop(seq_id_header)
                            sequences_new.append(seq)
                        sequences = sequences_new
                        # print(len(sequences_new))
                        # if seq_id in ids_added:
                        #     sequences.append(seq_dup)
                        ids_added.append(seq_id)
                    else:
                        sequences.append(consensus_block)
                else:
                    print(tag, attrib, value)
            return sequences
        elif has_children:
            for child in el:
                sequences += self.process_sequences(child, metadata=metadata)
            # df = pd.DataFrame(metadata)
            # df.to_sql('metadata', self.db, if_exists='append', index=True)
            # print('meta', metadata)
            return sequences
        else:
            return []

    def process_hml_file(self, 
            filepath : str = None, 
            xmlstring : str = None,
            outdir : str = None) -> List[Dict[str, str]]:
        if isinstance(xmlstring, bytes):
            if len(xmlstring):
                return self.process_sequences(ET.fromstring(xmlstring), filename=filepath)
        elif filepath:
            # print('c', filepath)
            # ET.parse(filepath).getroot()
            return self.process_sequences(ET.parse(filepath).getroot(), filename=filepath)
        

    def process_hml(self, filepath : Union[str, List[str]] = None, 
                xmlstring : str = None,
                outdir : str = None, 
                n_files_out : int = 0,
                n_print_out : int = 100,
                n_files : int = None,
                query : str = None, db : str = None,
                verbose : bool = False, stratify_loci : bool = False, pwd : Union[str, dict] = None) -> pd.core.frame.DataFrame:
        if db:
            self._create_db(db)
        start_time = time.time()
        last_file = None
        started = True

        if isinstance(filepath, str):
            filepath = [filepath]
        dfs = []
        filepaths = filepath
        for filepath in filepaths:
            if isinstance(pwd, dict):
                filepath_pwd = filepath
                if filepath_pwd not in pwd:
                    filepath_pwd = filepath_pwd.split('/')[-1]
                pwd = pwd[filepath_pwd]
            
            def handler(signum, frame):
                raise Exception("Ran too long!")

            if filepath:
                rows = []
                # print('Is {} zipped?'.format(filepath))
                is_zipped = '.zip' in filepath or is_zipfile(filepath)
                if os.path.isdir(filepath) or is_zipped:
                    if is_zipped:
                        # print('Obtaining files from zip:', filepath)
                        # start_time = time.time()
                        zf = ZipFile(filepath)
                        files = zf.namelist()
                        if database_exists(self.db.url) and (self.name not in ['hml_database', 'test']):
                            last_file = pd.read_sql_query("SELECT filename from 'table'", self.db).iloc[-1][0]
                            print('Database exists. Last file is ', last_file)
                        # last_file = 'HML201603-20200729-to-CIBMTR/1352.hml101.xml'
                        # print('Took {:,} minutes to obtain files.'.format((time.time() - start_time) / 60))
                    else:
                        files = os.listdir(filepath)
                    if verbose:
                        print('There are {:,} files'.format(len(files)))
                    n = 0
                    columns = []
                    if last_file:
                        started = False
                    for i, file in enumerate(files):
                        if n_files and (n_files < n):
                            break
                        # if i < 15190:
                        #     continue
                        if started and not file.startswith('.'):
                            filepath_input = xmlstring = None
                            if is_zipped:
                                if pwd and not isinstance(pwd, bytes):
                                    pwd = pwd.encode('utf-8')
                                xmlstring = zf.open(file, pwd=pwd).read()
                                filepath_input = file
                            else:
                                filepath_input = os.path.join(filepath, file)
                            signal.signal(signal.SIGALRM, handler)
                            signal.alarm(60)
                            try:
                                output = self.process_hml_file(filepath=filepath_input, 
                                                            xmlstring=xmlstring,
                                                                outdir=outdir)
                                if output:
                                    rows +=  output
                            except Exception as e:
                                print(i, file, e)
                                rows += [{'filename' : file, 'error' : str(e)}]
                            signal.alarm(0)
                            n = n + 1
                            # n += 1 THIS IS TRUE
                            if (n and n % n_files_out == 0) and outdir:
                                # print('PRINTED OUT', len(rows), i, n, n_files_out, file)
                                df = pd.DataFrame(rows)
                                if not columns:
                                    columns = list(df.columns)
                                for col in df.columns:
                                    if col not in columns:
                                        # print(col)
                                        self.add_column(col)
                                        # df[col] = ""
                                        columns.append(col)
                                # df = df[columns]
                                df.to_sql('table', self.db, if_exists='append' if n != n_files_out else 'replace',
                                    index=False)

                                # df.to_csv(outdir, mode='a' if n != n_files_out else 'w', index=False)
                                rows = []
                            if (n and n % n_print_out == 0) and verbose:
                                dur = time.time() - start_time
                                spf = dur / n
                                print('{}: Ran through {:,} files. {:,} left. {:.1f} secs have elapsed. {:.2f} secs per file. {:.1f} mins left.'\
                                    .format(filepath, n, len(files) - n, dur, spf, (len(files) - n) * spf / 60))
                        if last_file == file:
                            started = True
                elif os.path.isfile(filepath):
                    rows = self.process_hml_file(filepath, outdir=outdir)
                else:
                    raise Exception('This file cannot be processed.')
            elif xmlstring:
                rows = self.process_hml_file(xmlstring=xmlstring, outdir=outdir)
            output_df = pd.DataFrame(rows)
            dfs.append(output_df)
        output_df = pd.concat(dfs)
        if outdir:
            if n_files_out == 0:
                if stratify_loci:
                    for locus in set(output_df['sbt-ngs-locus']):
                        if locus in ['HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5']:
                            locus = 'HLA-DRB[345]'
                        if '.' in outdir:
                            filename, ext = outdir.split('.', 1)
                        else:
                            filename, ext = outdir, 'csv'
                        output_df[output_df['sbt-ngs-locus'].str.contains(locus)].\
                            to_csv('{}_{}.{}'.format(filename, locus, ext), index=False)
                else:
                    output_df.to_csv(outdir, mode='w', index=False)
            else:
                # output_df = pd.read_csv(outdir)
                output_df = self.parse_db( self.db, query, fileout=outdir)
        if verbose:
            print('Finished all files. {:.2f} secs have elapsed.'.format(time.time() - start_time))
        return output_df

    def parse_db(self, db : Union[str, sqlalchemy.engine.base.Engine], query : str = None, fileout : str = None) -> pd.core.frame.DataFrame:
        if isinstance(db, str):
            db = self._create_db(db)
        if not query:
            query = "SELECT * FROM 'table'"
        if db:
            self._create_db(db)
        df = pd.read_sql_query(query, db)
        print(type(db))
        if fileout:
            df.to_csv(fileout, index=False)
        return df

    def add_column(self, col : str) -> None:
        """
        Adds column to sqlite databse
        """
        # conn = sqlite3.connect(self.name)
        conn = self.db.connect()
        # c = conn.cursor()
        cmd = "ALTER TABLE 'table' ADD column '{}' 'VARCHAR(100)'".format(col)
        conn.execute(cmd)
        # c.execute(cmd)
        # conn.commit()
        # conn.close()

