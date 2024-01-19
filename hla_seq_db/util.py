#
# Copyright (c) 2024 NMDP.
#
# This file is part of HLA-DB 
# (see https://github.com/nmdp-bioinformatics/hla-db).
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
import pandas as pd
import re

def get_two_field_allele(allele_name : str) -> str:
    """
    Obtains allele name with only first two fields
    """
    return ':'.join(str(allele_name).split(':')[:2])

def calc_hla_class(locus : str) -> int:
    """
    Calculate the HLA class of this particular allele typing.
    """
    if locus in ['A', 'C','B']:
        return 1
    elif 'D' in locus:
        return 2
    else:
        raise InvalidAllotypeError(locus,
        "Check this locus. Cannot determine its HLA class")

def sort_df(df : pd.core.frame.DataFrame) -> pd.core.frame.DataFrame:
    df_headers = pd.core.frame.DataFrame([[float(val) for val in 
                                            re.sub(r'[A-Z]$', '.5', re.sub(r'.+\*', '', name)).split(':')]
                                 for name in df.index], index=df.index)
    df_headers = df_headers.sort_values(list(df_headers.columns))
    return df.loc[df_headers.index]