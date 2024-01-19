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
from hamcrest import assert_that, is_

@given('that the allele typing is "{allele_name}"')
def step_impl(context, allele_name):
    try:
        context.allele = context.hla_db.create_allele(allele_name)
    except:
        context.allele = None

@given('the associated SIRE is "{sire}"')
def step_impl(context, sire):
    context.allele.sire = sire

@when('evaluating the validity of the allele')
def step_impl(context):
    context.validity = 'valid' if context.allele else 'invalid'

@when('extracting the possible alleles')
def step_impl(context):
    allele_dict = {'allele_list' : context.allele.alleles or context.allele.get_potential_alleles(),
                   'allele_hi_res_list' : context.allele.alleles_hi_res}
    if context.allele:
        for allele_type, alleles in allele_dict.items():
            alleles = sorted([str(allele) for allele in alleles])
            allele_list = [str(allele) for allele in alleles[:3]]
            allele_dict[allele_type] = allele_list
    context.allele_dict = allele_dict

@when('annotating genomic features (TCE, expression)')
def step_impl(context):
    context.annotation = context.hla_db.annotate_allotype(context.allele).annotation
    context.alleles_obs = []
    for allele in context.allele.alleles_hi_res:
        allele_dict = {'allele_two_field' : allele.typing}
        allele_dict.update(allele.annotation)
        context.alleles_obs.append(allele_dict)

@then('the allele typing is found to be "{validity}"')
def step_impl(context, validity):
    assert_that(context.validity, is_(validity))

@then('the first three {allele_type} alleles are found to be "{allele_list}"')
def step_impl(context, allele_type, allele_list):
    allele_type = 'allele_hi_res_list' if allele_type == 'hi-res' else 'allele_list'
    assert_that(context.allele_dict[allele_type], is_(allele_list.split(',')))