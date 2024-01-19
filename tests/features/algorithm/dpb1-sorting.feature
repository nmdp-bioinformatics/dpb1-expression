Feature: Sorting DPB1 mismatches based on DPB1 expression and T-cell epitope (TCE) groups
    

    Scenario: Simple list with low/high expression and permissive/nonpermissive TCE matching

        Given that the genotype is "DPB1*01:01+DPB1*02:01" for a patient
        And the genotypes and ranks from the expected donor list
            | rank | Genotype              | TCE genotype | TCE match | Pat. mismatched allele | Expr Match   | Direction     |
            | 1    | DPB1*01:01+DPB1*04:01 | 3+3          | P         | DPB1*02:01 (L)         | Favorable    | bidirectional |
            | 2    | DPB1*01:01+DPB1*01:01 | 3+3          | P         | DPB1*02:01 (L)         | Favorable    | GvH           |
            | 3    | DPB1*01:01+DPB1*17:01 | 3+1          | NP (HvG)  | DPB1*02:01 (L)         | Favorable    | bidirectional |
            | 4    | DPB1*04:01+DPB1*02:01 | 3+3          | P         | DPB1*01:01 (H)         | Unfavorable  | bidirectional |
            | 5    | DPB1*02:01+DPB1*02:01 | 3+3          | P         | DPB1*01:01 (H)         | Unfavorable  | GvH           |
            | 6    | DPB1*17:01+DPB1*02:01 | 3+1          | NP (HvG)  | DPB1*01:01 (H)         | Unfavorable  | bidirectional |

        When evaluating and sorting the donor list
        Then the patient's 'expression' genotype is found to be 'high+~low'
        And the patient's 'tce' genotype is found to be '3+3'
        And the expected and observed donor lists are the same

    # Scenario: Simple list with low/high expression and permissive/nonpermissive TCE matching

    #     Given that the genotype is "DPB1*02:01+DPB1*03:01" for a patient
    #     And the genotypes and ranks from the expected donor list
    #         | rank | Genotype              | TCE genotype | TCE match | Pat. mismatched allele | Expr Match   | Direction     |
    #         | 1    | DPB1*04:01+DPB1*03:01 | 2+3          | P         | DPB1*02:01 (L)         | Favorable    | bidirectional |
    #         | 2    | DPB1*03:01+DPB1*03:01 | 2+2          | P         | DPB1*02:01 (L)         | Favorable    | GvH           |
    #         | 3    | DPB1*17:01+DPB1*03:01 | 2+1          | NP (HvG)  | DPB1*02:01 (L)         | Favorable    | bidirectional |
    #         | 4    | DPB1*02:01+DPB1*14:01 | 3+2          | P         | DPB1*03:01 (H)         | Unfavorable  | bidirectional |
    #         | 5    | DPB1*02:01+DPB1*09:01 | 3+1          | NP (HvG)  | DPB1*03:01 (H)         | Unfavorable  | bidirectional |
    #         | 6    | DPB1*02:01+DPB1*04:01 | 3+3          | NP (GvH)  | DPB1*03:01 (H)         | Unfavorable  | bidirectional |
    #         | 7    | DPB1*02:01+DPB1*02:01 | 3+3          | NP (GvH)  | DPB1*03:01 (H)         | Unfavorable  | GvH           |
    #     When evaluating and sorting the donor list
    #     Then the patient's 'expression' genotype is found to be 'L+H'
    #     And the patient's 'tce' genotype is found to be '3+2'
    #     And the expected and observed donor lists are the same