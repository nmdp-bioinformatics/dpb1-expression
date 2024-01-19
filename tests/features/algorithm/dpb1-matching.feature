Feature: Define HLA-DPB1 Genotype Match
    
    A genotype match is defined by its HLA allotype match grades as folfavorable:
        AA - A genotype match occurs if both allotype pairs are allele matches.
        MA - A single mismatch constitutes one allele match and one mismatch.
        MM - A double mismatch occurs when both allotype pairs are mismatches.
        PA/PP/LP/MP - A potential single/double genotype match occurs when at least one pair has a potential match.
        LA/LP/LL/ML - An allele single/double genotype mismatch occurs when at least one pair has an allele mismatch.
    
    Scenario Outline: HLA-DPB1 Genotype Matches

        Given two HLA-DPB1 genotype names as <Genotype Pat> and <Genotype Don> for patient and donor, respectively
        When evaluating the TCE match, DPB1 match grades, directionality, and matched/mismatch alleles between the two genotypes
        Then the donor's genotype was <Flipped>
        And the match grades are found to be <Match Grades>
        And the expression match is found to be <Expr Level>
        And the directionality is found to be <Directionality>
        And the matched alleles are <Matched Alleles Patient> and <Matched Alleles Donor> for patient and donor, respectively
        And the mismatched alleles are <Mismatched Alleles Patient> and <Mismatched Alleles Donor> for patient and donor, respectively

        Examples: Hi-res genotype Match
            |  Genotype Pat               |  Genotype Don               | Flipped  | Match Grades | Expr Level | Directionality | Matched Alleles Patient | Matched Alleles Donor | Mismatched Alleles Patient | Mismatched Alleles Donor |
            | DPB1*09:01+DPB1*40:01       | DPB1*09:01+DPB1*40:01       | unflipped |      AA      | None       | None           | DPB1*09:01,DPB1*40:01   | DPB1*09:01,DPB1*40:01 | None                       | None                     |
            | DPB1*09:01:01+DPB1*40:01:01 | DPB1*09:01+DPB1*40:01       | unflipped |      AA      | None       | None           | DPB1*09:01:01,DPB1*40:01:01 | DPB1*09:01,DPB1*40:01 | None                       | None                     |
            | DPB1*09:01+DPB1*40:01       | DPB1*40:01+DPB1*09:01       | flipped   |      AA      | None       | None           | DPB1*09:01,DPB1*40:01   | DPB1*09:01,DPB1*40:01 | None                       | None                     |
            | DPB1*40:01+DPB1*09:01       | DPB1*40:01+DPB1*09:01       | unflipped |      AA      | None       | None           | DPB1*40:01,DPB1*09:01   | DPB1*40:01,DPB1*09:01 | None                       | None                     |
            | DPB1*04:AETTB+DPB1*13:KHMN  | DPB1*04:AETTB+DPB1*13:KHMN  | unflipped |      PP      | None       | None           | DPB1*04:AETTB,DPB1*13:KHMN  | DPB1*04:AETTB,DPB1*13:KHMN | None                       | None                     |
            | DPB1*02:AEMGJ+DPB1*09:01:01 | DPB1*02:AEMGJ+DPB1*09:01:01 | unflipped |      PA      | None       | None           | DPB1*02:AEMGJ,DPB1*09:01:01 | DPB1*02:AEMGJ,DPB1*09:01:01 | None                       | None                     |
        
        Examples: Single Mismatch
            |  Genotype Pat               |  Genotype Don              | Flipped  | Match Grades | Expr Level | Directionality | Matched Alleles Patient | Matched Alleles Donor | Mismatched Alleles Patient | Mismatched Alleles Donor |
            | DPB1*09:01+DPB1*40:01       | DPB1*08:01+DPB1*40:01      | unflipped |      MA      | unfavorable       | bidirectional  | DPB1*40:01              | DPB1*40:01            | DPB1*09:01                | DPB1*08:01              |
            | DPB1*09:01+DPB1*40:01       | DPB1*02:01+DPB1*09:01      | flipped   |      AM      | favorable       | bidirectional  | DPB1*09:01              | DPB1*09:01            | DPB1*40:01                | DPB1*02:01              |
            # | DPB1*40:01+DPB1*40:01       | DPB1*09:01+DPB1*40:01      | unflipped |      MA      | favorable       | HvG            | DPB1*40:01,DPB1*40:01   | DPB1*40:01            | None                      | DPB1*09:01              |
            # | DPB1*09:01+DPB1*40:01       | DPB1*40:01+DPB1*40:01      | unflipped |      MA      | unfavorable       | GvH            | DPB1*40:01              | DPB1*40:01,DPB1*40:01 | DPB1*09:01                | None                    |
            | DPB1*47:01+DPB1*09:01       | DPB1*132:01+DPB1*09:01     | unflipped |      MA      | favorable       | bidirectional  | DPB1*09:01              | DPB1*09:01            | DPB1*47:01                | DPB1*132:01              |
            | DPB1*40:01+DPB1*132:01      | DPB1*35:01+DPB1*132:01     | unflipped |      MA      | favorable       | bidirectional  | DPB1*132:01             | DPB1*132:01            | DPB1*40:01                | DPB1*35:01              |
            | DPB1*04:01+DPB1*04:FNVS     | DPB1*02:01+DPB1*04:02      | unflipped |      MP      | favorable       | bidirectional  | DPB1*04:FNVS            | DPB1*04:02            | DPB1*04:01                | DPB1*02:01              |
            | DPB1*01:AETTA+DPB1*04:AETTB | DPB1*04:AETTB+DPB1*04:FNVS | flipped   |      MP      | unfavorable       | bidirectional  | DPB1*04:AETTB           | DPB1*04:AETTB          | DPB1*01:AETTA             | DPB1*04:FNVS           |
            # | DPB1*02:01+DPB1*04:VM       | DPB1*04:VR+DPB1*04:VR      | unflipped |      MP      | favorable       | GvH            | DPB1*04:VM              | DPB1*04:VR,DPB1*04:VR  | DPB1*02:01                | None                   |
            | DPB1*04:VZS+DPB1*11:01      | DPB1*04:01+DPB1*11:01      | unflipped |      LA      | favorable       | bidirectional  | DPB1*11:01              | DPB1*11:01             | DPB1*04:VZS               | DPB1*04:01             | 
            # AMBIGUOUS TCE | DPB1*01:01+DPB1*02:01       | DPB1*01:ANZW+DPB1*17:ANZX  | unflipped |      PM      | bidirectional  | DPB1*01:01              | DPB1*01:ANZW           | DPB1*02:01                | DPB1*17:ANZX           |

        Examples: Double Mismatch
            |  Genotype Pat         |  Genotype Don         | Flipped  | Match Grades | Expr Level  | Directionality | Matched Alleles Patient | Matched Alleles Donor | Mismatched Alleles Patient | Mismatched Alleles Donor |
            | DPB1*142:01+DPB1*05:01 | DPB1*09:01+DPB1*06:01 | unflipped |      MM      | None       | bidirectional | None                    | None                  | DPB1*142:01,DPB1*05:01     | DPB1*09:01,DPB1*06:01   |
            | DPB1*05:01+DPB1*142:01 | DPB1*09:01+DPB1*06:01 | unflipped |      MM      | None       | bidirectional | None                    | None                  | DPB1*05:01,DPB1*142:01     | DPB1*09:01,DPB1*06:01   |

        Examples: Potential Single/Double Genotype Match
            | Genotype Pat                   |  Genotype Don             | Flipped  | Match Grades | Expr Level | Directionality | Matched Alleles Patient | Matched Alleles Donor | Mismatched Alleles Patient | Mismatched Alleles Donor |
            | DPB1*01:01:01G+DPB1*51:01:01G  | DPB1*01:01:01G+DPB1*09:01 | unflipped |      PM      | favorable       | bidirectional  | DPB1*01:01:01G          | DPB1*01:01:01G        | DPB1*51:01:01G            | DPB1*09:01              |
            # | DPB1*01:01:01G+DPB1*01:01:01G  | DPB1*01:01:01G+DPB1*09:01 | unflipped |      PM      | ~unfavorable      | HvG            | DPB1*01:01:01G,DPB1*01:01:01G | DPB1*01:01:01G  | None                      | DPB1*09:01              |
            # | DPB1*44:01+DPB1*15             | DPB1*44:01+DPB1*15:01     | unflipped |      AP      | None       | None           | DPB1*44:01              | DPB1*44:01            | DPB1*15                   | DPB1*15:01              |
            | DPB1*09+DPB1*09                | DPB1*09+DPB1*09           | unflipped |      AA      | None       | None           | DPB1*09,DPB1*09         | DPB1*09,DPB1*09       | None                      | None                    |
            | DPB1*04:02+DPB1*09             | DPB1*09:01+DPB1*04:01     | flipped   |      LA      | favorable       | bidirectional  | DPB1*09                 | DPB1*09:01            | DPB1*04:02                | DPB1*04:01              |
            | DPB1*15+DPB1*09                | DPB1*42+DPB1*15           | flipped |      PM        | unfavorable       | bidirectional  | DPB1*15                 | DPB1*15               | DPB1*09                   | DPB1*42                 |

        # Examples: Single/Double Allele Genotype Mismatch
        #     |  Genotype Pat   |  Genotype Don   | Flipped  | Match Grades | Directionality |
        #     | DPB1*09:01+DPB1*05:01 | DPB1*09:04+DPB1*05:01 | unflipped |      LA      | None |
        #     | DPB1*14:01+DPB1*09:01 | DPB1*14:03+DPB1*09    | unflipped |      LP      | None |
        #     | DPB1*14:01+DPB1*09:01 | DPB1*14:03+DPB1*09:12 | unflipped |      LL      | None |
        #     | DPB1*42:01+DPB1*09:01 | DPB1*40:01+DPB1*09:03 | unflipped |      ML      | None |