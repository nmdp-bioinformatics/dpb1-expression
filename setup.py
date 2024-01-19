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

from setuptools import setup, find_packages
from setuptools_behave import behave_test

requirements = ['allure-behave==2.8.5',
                'behave==1.2.6',
                'biopython==1.75',
                'Flask==1.1.1',
                'flask-restx==0.5.0',
                'regex==2021.7.6',
                'requests==2.26.0',
                'py-ard==0.6.6',
                'py2neo==2021.0.1',
                'PyHamcrest==1.9.0']

test_requirements = [
]

setup(name='dpb1_expr',
      version='0.1',
      description='Package for annotating the expression levels and TCE groups of HLA-DPB1 typing.',
      url='https://github.com/nmdp-bioinformatics/dpb1-expression',
      author='NMDP Bioinformatics, CIBMTR',
      author_email='bioinformatics-web@nmdp.org',
      license='MIT',
      packages=['dpb1_expr', 'hla_seq_handler.hla_seq_db'],
      install_requires=requirements,
      tests_require=["behave>=1.2.6"],
      keywords='dpb1_expr')