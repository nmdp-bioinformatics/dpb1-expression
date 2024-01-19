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
import behave
# from behave import *
from hamcrest import assert_that, is_
from hla_seq_db.genotype import Genotype
from hla_seq_db.match_list import MatchList
from typing import Union, Tuple, List

@given('that the genotype is "{genotype}" for a patient')
def step_impl(context, genotype):
    context.patient = context.hla_db.create_genotype(genotype)

def format_ranks(rows : Union[behave.model.Table, MatchList]) -> List[Tuple[int, str]]:
    if isinstance(rows, MatchList):
        return([(match.rank, str(match.genotype_donor)) for match in rows.matches])
    elif isinstance(rows, behave.model.Table):
        return [(int(row['rank']), 
                    row['genotype_donor']['genotype']
                        if 'genotype_donor' in row 
                        else row['Genotype'])
                        for row in rows]

@given('the genotypes and ranks from the expected donor list')
def step_impl(context):
    context.ranks = format_ranks(context.table)
    context.donors = [context.hla_db.create_genotype(row["Genotype"]) for row in context.table]

@when('evaluating and sorting the donor list')
def step_impl(context):
    context.match_list = context.hla_db.annotate_match_list(context.patient,
                                          context.donors)
    context.match_list.sort()

@then("the patient's '{result}' genotype is found to be '{genotype}'")
def step_impl(context, result, genotype):
    geno_result = '+'.join(context.match_list.geno_recip.annotation[result])
    assert_that(geno_result, is_(genotype))


@then('the expected and observed donor lists are the same')
def step_impl(context):
    assert_that(format_ranks(context.match_list), is_(context.ranks))