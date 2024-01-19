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

@then('the {output} is found to be "{value}"')
def step_impl(context, output : str, value: str):
    map = {'expression' : 'expression',
           'expression confirmation' : 'expression_experimental',
           'TCE group' : 'tce',
           }
    output = context.annotation[map[output]]
    if output == True:
        output = 'True'
    elif output == False:
        output = 'False'
    assert_that(output or 'None', is_(value))

@given('the possible hi-res alleles and expression motifs.')
def step_impl(context):
    results = []
    for row in context.table:
        result = {'allele_two_field' : row['Allele typing'],
                  'CIWD_{}_combined'.format('TOTAL') : row['CIWD'],
                  'expression' : row['Expression'],
                  'expression_experimental' : row['Experimentally Confirmed?'],
                  'tce' : row['TCE']
                  }
        results.append(result)
    context.alleles = results

@then('the possible hi-res alleles and expression motifs are as expected.')
def step_impl(context):
    assert_that(context.alleles_obs, is_(context.alleles))