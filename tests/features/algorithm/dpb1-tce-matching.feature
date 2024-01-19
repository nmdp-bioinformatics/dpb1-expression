Feature: Define HLA-B Genotype Match
    
    Scenario Outline: HLA-B Genotype Matches

        Given two HLA-DPB1 genotype names as <Genotype Pat> and <Genotype Don> for patient and donor, respectively
        When evaluating the TCE match, DPB1 match grades, directionality, and matched/mismatch alleles between the two genotypes
        Then the TCE match category is found to be <Match Category>

        Examples: Allele
            | Genotype Pat          |  Genotype Don         | TCE Genotype Pat | TCE Genotype Don | Match Category    |
            | DPB1*03:01+DPB1*40:01 | DPB1*03:01+DPB1*40:01 | 3+3              | 3+3              | Allele            |

        Examples: Permissive
            | Genotype Pat          |  Genotype Don         | TCE Genotype Pat | TCE Genotype Don | Match Category    |
            | DPB1*03:01+DPB1*40:01 | DPB1*08:01+DPB1*40:01 | 2+3              | 2+3              | Permissive |
            
        Examples: Nonpermissive (graft-versus-host)
            | Genotype Pat          |  Genotype Don         | TCE Genotype Pat | TCE Genotype Don | Match Category |
            | DPB1*08:01+DPB1*40:01 | DPB1*40:01+DPB1*40:01 | 2+3              | 3+3              | GvH_nonpermissive     |
            
        
        Examples: Nonpermissive (host-versus-graft)
            | Genotype Pat          |  Genotype Don         | TCE Genotype Pat | TCE Genotype Don | Match Category    |
            | DPB1*04:01+DPB1*40:01 | DPB1*03:01+DPB1*40:01 | 3+3              | 2+3              | HvG_nonpermissive |