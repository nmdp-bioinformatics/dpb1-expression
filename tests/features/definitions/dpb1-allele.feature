Feature: Define HLA-DPB1 Allele
    HLA-DPB1 is a key gene in the major histocompability complex that codes
    for the β subunit in the class II DP protein. It is key for 
    hematopoietic cell transplantation via antigen recognition function in antigen-presenting cells.
    The following scenarios established the correct nomenclature for representing HLA-DPB1 alleles.
    Alleles represent alternate versions of a specific gene.
    
    Scenario Outline: Allele name and fields
        Valid alleles names start with `DPB1*` or `HLA-DPB1*`, followed by a list of
        numeric field IDs separated by colons. MACs (Multiple Allele Codes) such as
        DPB1*03:CXD are valid as well as G groups such as DPB1*01:01:02G.

        Given that the allele typing is "<Allele typing>"
        When evaluating the validity of the allele
        Then the allele typing is found to be "<Validity>"

        Examples: Valid examples of HLA-DPB1 alleles
            | Allele typing  | Validity |
            | DPB1*02:01     |  valid   |
            | HLA-DPB1*02:01 |  valid   |
            | DPB1*03:CXD    |  valid   |
            | DPB1*01:01:02G |  valid   |

        Examples: Invalid examples of HLA-DPB1 alleles
            | Allele typing | Validity |
            | 😄            | invalid  |
            | DPB1*02:01+DPB1*02:01     | invalid  |
    
    # Scenario Outline: Allele within G group

    #     Given that the allele typing is "<Allele typing>"
    #     When obtaining the G group it might be within
    #     Then the G group is found to be "<G group>"

    #     Examples: Examples of alleles within G groups
    #         | Allele typing    | G group        |
    #         | DPB1*01:01:02:01 | DPB1*01:01:02G |

    Scenario Outline: Ambiguous Typing

        Given that the allele typing is "<Allele typing>"
        When extracting the possible alleles
        Then the first three allelic-res alleles are found to be "<Allele list>"
        And the first three hi-res alleles are found to be "<Hi-res allele list>"

        Examples: Examples of ambiguous typing
            | Allele typing  | Allele list                                        | Hi-res allele list    |
            | DPB1*02:01     | DPB1*02:01:02:01,DPB1*02:01:02:02,DPB1*02:01:02:03 | DPB1*02:01            |
            | HLA-DPB1*02:01 | DPB1*02:01:02:01,DPB1*02:01:02:02,DPB1*02:01:02:03 | DPB1*02:01            | 
            | DPB1*03:CXD    | DPB1*03:01:01:01,DPB1*03:01:01:02,DPB1*03:01:01:03 | DPB1*03:01,DPB1*78:01 |
            | DPB1*01:01:02G | DPB1*01:01:02:01,DPB1*01:01:02:02,DPB1*01:01:02:03 | DPB1*01:01,DPB1*1050:01,DPB1*162:01 |