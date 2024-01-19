# dpb1-expression
DPB1 expression


This program creates DPB1 expression table.

## Requirements

1. An accessible instance of GFEDB loaded with the target HLADB version 
2. py2neo

## instructions
Setup:
```
virtualenv -p /usr/local/homebrew/bin/python3 venv
source venv/bin/activate
pip install py2neo
pip install pandas
```
Run:
python extract.py 


Output:
dpb1_expression.3.42.0.csv

Sample:
```
allele,motif,high,low
HLA-DPB1*994:01,ACCACTC,1.0,0.2857142857142857
...
```

# Interpretation

The columns "high" and "low" provide the sequence similarity to the exon 3 positions identified in [this paper](https://doi.org/10.1016/j.humimm.2017.11.001)

These values can form the basis for a heuristic that classifies the allele as being high or low expression when the similarity to the high or low motif is above a threshold (e.g. 0.8).

Note: some alleles do not have a similarity score of 0.8 to either motif.  These include ```HLA-DPB1*10:01:04``` which is equally distant to both patterns and these 5 alleles (as of HLADB version 3.42.0) all of which happen to have null or questionable expression:
```
HLA-DPB1*986:01N
HLA-DPB1*876:01N
HLA-DPB1*866:01N
HLA-DPB1*878:01N
HLA-DPB1*1148:01Q
```

# TODO

This method only works for HLA-DPB1 alleles that are defined in exon 3.
It remains an open issue to make this assignment for alleles that are not defined in exon 3 (e.g. exon2 only).
