#
# Copyright (c) 2024 NMDP.
#
# This file is part of DPB1 Expression 
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
import sys
import os
import argparse
from py2neo import Graph
import pandas as pd
from difflib import SequenceMatcher


# high expression pattern
high_pat = "ACCACTC"
# low expression pattern
low_pat = "GTTGTCT"

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

dpb1_expression_cypher = """
MATCH (n:IMGT_HLA)-
[hf:HAS_FEATURE]-(f:FEATURE)
WHERE hf.imgt_release = "{imgt}"
    and n.locus = "HLA-DPB1" 
    and f.name = "EXON"
    and f.rank = "3"
RETURN n.name as allele,
substring(f.sequence, 374-365,1)+
substring(f.sequence, 381-365,1)+
substring(f.sequence, 406-365,1)+
substring(f.sequence, 441-365,1)+
substring(f.sequence, 588-365,1)+
substring(f.sequence, 596-365,1)+
substring(f.sequence, 624-365,1) as motif
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--address",
        required=False, default="localhost",
        help="address of the GFE neo4j database server. eg: neo4j.b12x.org",
        type=str)
    parser.add_argument("-u", "--user",
        required=False, default="neo4j",
        help="username on neo4j db",
        type=str)
    parser.add_argument("-i", "--imgt",
        required=False, default="3.42.0",
        help="IMGT/HLADB version (e.g. 3.42.0)", 
        type=str)
    args = parser.parse_args()
    dpb1_expression_cypher = dpb1_expression_cypher.format(imgt=args.imgt)
    print ("dpb1_expression_cypher: ", dpb1_expression_cypher)
    # neo4j
    neo4j_server = args.address
    neo4j_user = args.user

    # neo4j password
    neo4j_password = os.getenv('NEO4J_PASSWORD')

    if neo4j_password is None or neo4j_password == "password":
        print("Exiting. NEO4J_PASSWORD not supplied")
        sys.exit(1)

    outfile = "dpb1_expression.{imgt}.csv".format(imgt=args.imgt)
    neo4j_server = "http://{us}:{pw}@neo4j.b12x.org:80".format(us=neo4j_user, pw=neo4j_password)
    
    print ("connecting using host= ", neo4j_server)
    graph = Graph(neo4j_server)
    df = graph.run(dpb1_expression_cypher).to_data_frame()

    # match 
    df['high'] = df.apply(lambda row: similar(high_pat, row['motif']), axis=1)
    df['low'] = df.apply(lambda row: similar(low_pat, row['motif']), axis=1)
    print ("output to: ", outfile)
    df.to_csv(outfile, index=False) 


