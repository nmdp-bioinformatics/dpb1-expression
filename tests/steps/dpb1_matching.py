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
from behave import *
from hamcrest import assert_that, is_
from typing import List
from hla_seq_db.allele import Allele

@given('two HLA-DPB1 genotype names as {genotype_name_one} and {genotype_name_two} for patient and donor, respectively')
def step_impl(context, genotype_name_one, genotype_name_two):
    context.genotype_one = genotype_name_one
    context.genotype_two = genotype_name_two

def concatenate(alleles : List[Allele]) -> str:
    """
    Contatenates a list the string names of Allele objects using a comma.
    Returns 'None' if the list is empty.
    """
    return ','.join([str(allele) for allele in alleles]) or 'None'

@when('evaluating the TCE match, DPB1 match grades, directionality, and matched/mismatch alleles between the two genotypes')
def step_impl(context):
    context.match = context.hla_db.annotate_match(context.genotype_one, context.genotype_two).serialize()
    context.matched_alleles_pat = concatenate(context.match['matched_alleles_pat'])
    context.matched_alleles_don = concatenate(context.match['matched_alleles_don'])
    context.mismatched_alleles_pat = concatenate(context.match['mismatched_alleles_pat'])
    context.mismatched_alleles_don = concatenate(context.match['mismatched_alleles_don'])

@then("the donor's genotype was {flipped}")
def step_impl(context, flipped):
    assert_that(str(context.match['genotype_donor_flipped'] and 'flipped' or 'unflipped'), 
        is_(flipped))

@then('the match grades are found to be {match_grades}')
def step_impl(context, match_grades):
    assert_that(str(context.match['grade']), is_(match_grades))

@then("the expression match is found to be {expr_level}")
def step_impl(context, expr_level):
    print(context.match)
    assert_that(str(context.match['annotation']['expr_match']), is_(expr_level))


@then('the directionality is found to be {directionality}')
def step_impl(context, directionality):
    assert_that(str(context.match['directionality']), is_(directionality))

@then('the matched alleles are {matched_alleles_pat} and {matched_alleles_don} for patient and donor, respectively')
def step_impl(context, matched_alleles_pat, matched_alleles_don):
    assert_that(str(context.matched_alleles_pat), is_(matched_alleles_pat))
    assert_that(str(context.matched_alleles_don), is_(matched_alleles_don))

@then('the mismatched alleles are {mismatched_alleles_pat} and {mismatched_alleles_don} for patient and donor, respectively')
def step_impl(context, mismatched_alleles_pat, mismatched_alleles_don):
    assert_that(str(context.mismatched_alleles_pat), is_(mismatched_alleles_pat))
    assert_that(str(context.mismatched_alleles_don), is_(mismatched_alleles_don))