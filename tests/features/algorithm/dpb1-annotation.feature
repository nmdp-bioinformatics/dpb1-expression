Feature: DPB1 Annotation (Expression + TCE)
    HLA-DP expression levels are highly correlated with the rs9277534 single nucleotide polymorphism (snp)
    in the 3' untranslated region (UTR). Since the 3' UTR is not routinely genotyped, this expression marker
    is not always accessible. Instead, standard genotyping of exons 2 and exons 3 can be used to impute the expression marker.
    More specifically, seven nucleotide positions (9, 16, 41, 76, 223, 231, and 259) on exon 3 were found to have 100% linkage
    to the rs9277534 expression marker by Schöne et al.(https://doi.org/10.1016/j.humimm.2017.11.001). At these seven positions,
    a motif of ACCACTC correlates with a G-rs9277534 expression marker and high expression. Conversely, a motif of GTTGTCT 
    correlates with an A-rs9277534 expression marker and low expression. Additionally, Schöne et al
    found that specific amino acid residues on exon 2 had 99% linkage with the expression marker. The following scenarios
    go over unambiguous and ambiguous DPB1 expression assignment based on these motifs on exon 2 and exon 3 as well as rs9277534 itself.

    T-Cell Epitope (TCE) groups are important for identifying permissive vs. non-permissive
    DPB1 mismatches.
    
    Scenario Outline: Allelic resolution
        
        Given that the allele typing is "<Allele typing>"
        When annotating genomic features (TCE, expression)
        Then the expression is found to be "<Expression>"
        And the expression confirmation is found to be "<Experimentally Confirmed?>"

        Examples: If the rs9277534 expression marker is known, then the expression level is known.
            | Allele typing    | exon3 motif | rs9277534 | Expression | Experimentally Confirmed? |
            | DPB1*01:01:01:01 | ACCACTC     | G         | high       | yes |
            | DPB1*02:01:02:01 | GTTGTCT     | A         | low        | yes |

        Examples: If the rs9277534 expression marker is unknown, but exon 3 sequence feature is known, then expression level is still known.
            | Allele typing    | exon3 motif | rs9277534 | Expression | Experimentally Confirmed? |
            | DPB1*01:01:01:05 | ACCACTC     | Unknown   | high       | yes |
            | DPB1*02:01:02:41 | GTTGTCT     | Unknown   | low        | yes |

        Examples: If the rs9277534 expression marker and exon 3 sequence feature are unknown, but exon 2 sequence feature is known
            | Allele typing    | exon3 motif | rs9277534 | Expression | Experimentally Confirmed? |
            | DPB1*03:01:02    | Unknown     | Unknown   | unknown    | yes |
            | DPB1*100:01      | Unknown     | Unknown   | unknown    | no |

    Scenario: Ambiguous Multiple Allele Code (MAC) with ambiguous expression marker

        Given that the allele typing is "DPB1*01:AETTG"
        And the possible hi-res alleles and expression motifs.
            | Allele typing    | CIWD    | exon3 motif | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
            | DPB1*01:01       | C/R/WD  | ACCACTC     | G         | 3   | high       | high       | yes         |
            | DPB1*162:01      | WD      | GTTGTCT     | A         | 3   | high       | low        | no        |
            | DPB1*417:01      | WD      | ACCACTC     | unknown   | 3   | high       | high       | no        |
        When annotating genomic features (TCE, expression)
        Then the expression is found to be "~high"
        And the expression confirmation is found to be "~yes"
        And the TCE group is found to be "3"
        And the possible hi-res alleles and expression motifs are as expected.

    # Scenario: Ambiguous G group with ambiguous expression marker

    #     Given that the allele typing is "DPB1*02:02:01G"
    #     And the possible hi-res alleles and expression motifs.
    #         | Allele typing    | CIWD    | exon3 motif | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
    #         | DPB1*02:02       | C/R     | GTTGTCT     | A         | 3   | high       | low        | yes         |
    #         | DPB1*547:01      | WD      | ACCACTC     | unknown   | 3   | allelic    | high       | no        |
    #         | DPB1*721:01      | unknown | GTTGTCT     | unknown   | 3   | allelic    | low        | no        |
    #         | DPB1*766:01      | unknown | GTTGTCT     | unknown   | 3   | allelic    | low        | no        |
    #         | DPB1*1188:01     | unknown | GTTGTCT     | A        | unknown | allelic | low       | no        |
    #     When annotating genomic features (TCE, expression)
    #     Then the expression is found to be "~low"
    #     And the expression confirmation is found to be "~yes"
    #     And the TCE group is found to be "3"
    #     And the possible hi-res alleles and expression motifs are as expected.

    # Scenario: Ambiguous G group with ambiguous expression marker

    #     Given that the allele typing is "DPB1*01:01:01G"
    #     And the possible hi-res alleles and expression motifs.
    #         | Allele typing    | CIWD    | exon3 motif | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
    #         | DPB1*01:01  | C | ACCACTC | G | 3 | high | high | yes |
    #         | DPB1*462:01 | WD | GTTGTCT | unknown | 3 | allelic | low | no |
    #         | DPB1*162:01 | WD | GTTGTCT | A | 3 | high | low | no |
    #         | DPB1*417:01 | WD | ACCACTC | unknown | 3 | high | high | no |
    #         | DPB1*616:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*733:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*800:01N | unknown | ACCACTC | unknown | 0 | allelic | high | no |
    #         | DPB1*807:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*810:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*862:01N | unknown | ACCACTC | unknown | 0 | allelic | high | no |
    #         | DPB1*953:01 | unknown | ACCACTC | G | 3 | allelic | high | no |
    #         | DPB1*979:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*998:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*999:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*1024:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*1038:01Q | unknown | ACCACTC | unknown | unknown | allelic | high | no |
    #         | DPB1*1068:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*1076:01 | unknown | ACCACTC | unknown | 3 | allelic | high | no |
    #         | DPB1*1151:01 | unknown | ACCACTC | unknown | unknown | allelic | high | no |
    #         | DPB1*1162:01 | unknown | ACCACTC | unknown | unknown | allelic | high | no |
    #         | DPB1*1183:01 | unknown | ACCACTC | unknown | unknown | allelic | high | no |
    #         | DPB1*1256:01N | unknown | ACCACTC | G | 0 | allelic | high | no |
    #         | DPB1*1050:01 | unknown | ACCACTC | unknown | 3 | high | high | no |
    #     When annotating genomic features (TCE, expression)
    #     Then the expression is found to be "~high"
    #     And the expression confirmation is found to be "~yes"
    #     And the TCE group is found to be "3"
    #     And the possible hi-res alleles and expression motifs are as expected.

    # Scenario Outline: Null alleles? How about null alleles?

    Scenario: Ambiguous Multiple Allele Code (MAC) with highly ambiguous expression marker

        Given that the allele typing is "DPB1*03:ACMGK"
        And the possible hi-res alleles and expression motifs.
            | Allele typing    | CIWD | exon3 motif          | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
            | DPB1*03:01       | C/I/R/WD | ACCACTC/ACCATTC/ACCACTT | G         | 2   | high       | high      | yes |
            | DPB1*104:01      | C/I | ACCACTC/ATCACTC         | G         | 2   | high       | high      | no |
            | DPB1*124:01      | C/I | GTTGTCT                 | A         | 2   | high       | low       | no |
            | DPB1*351:01      | WD | GTTGTCT                | unknown   | 2   | allelic    | low       | no |
        When annotating genomic features (TCE, expression)
        Then the expression is found to be "?high"
        And the expression confirmation is found to be "?yes"
        And the TCE group is found to be "2"
        And the possible hi-res alleles and expression motifs are as expected.

    # Scenario: Ambiguous Multiple Allele Code (MAC) with highly ambiguous expression marker. SIRE (self-identified race & ethnicity) provided.

    #     Given that the allele typing is "DPB1*03:ACMGK"
    #     And the associated SIRE is "API"
    #     And the possible hi-res alleles and expression motifs.
    #         | Allele typing    | CIWD | exon3 motif          | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
    #         | DPB1*03:01       | C | ACCACTC/ACCATTC/ACCACTT | G         | 2   | high       | high      | yes |
    #         | DPB1*104:01      | C | ACCACTC/ATCACTC         | G         | 2   | high       | high      | no |
    #         | DPB1*124:01      | I | GTTGTCT                 | A         | 2   | high       | low       | no |
    #         | DPB1*351:01      | R | GTTGTCT                | unknown   | 2   | allelic    | low       | no |
    #     When annotating genomic features (TCE, expression)
    #     Then the expression is found to be "~high"
    #     And the expression confirmation is found to be "~yes"
    #     And the TCE group is found to be "2"
    #     And the possible hi-res alleles and expression motifs are as expected.

    Scenario: Ambiguous Multiple Allele Code (MAC) with highly ambiguous TCE group

        Given that the allele typing is "DPB1*09:TH"
        And the possible hi-res alleles and expression motifs.
            | Allele typing    | CIWD     | exon3 motif  | rs9277534 | TCE | Resolution | Expression | Experimentally Confirmed? |
            | DPB1*09:01       | C/I/R/WD | ACCACTC/ACCATTC | G         | 1   | high       | high       | yes |
            | DPB1*13:01       | C/I/WD   | ACCACTC/ACCATTC | G         | 3   | high       | high       | yes |
        When annotating genomic features (TCE, expression)
        Then the expression is found to be "high"
        And the expression confirmation is found to be "yes"
        And the TCE group is found to be "?3"
        And the possible hi-res alleles and expression motifs are as expected.