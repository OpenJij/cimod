# Copyright 2020 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import cxxcimod
import cimod
from cimod.vartype import to_cxxcimod
from cimod.utils.decolator import recalc
import dimod
import numpy as np

def make_BinaryPolynomialModel(polynomial):
    """BinaryPolynomialModel factory.
       Generate BinaryPolynomialModel class with the base class specified by the arguments linear and quadratic
    Args:
        polynomial (dict): polynomial bias including linear bias
    Returns:
        generated BinaryPolynomialModel class
    """
    # select base class
    index = set()
    base = None

    if polynomial != {}:
        first_tuple_size = len(next(iter(polynomial)))
        for i in range (first_tuple_size):
            index.add(next(iter(polynomial))[i])

    if len(set(type(i) for i in index)) != 1:
        raise TypeError("invalid types of linear and quadratic")
    else:

        ind = next(iter(index))

        if isinstance(ind, int):
            base = cxxcimod.BinaryPolynomialModel
        elif isinstance(ind, str):
            base = cxxcimod.BinaryPolynomialModel_str
        elif isinstance(ind, tuple):
            if len(ind) == 2:
                base = cxxcimod.BinaryPolynomialModel_tuple2
            elif len(ind) == 3:
                base = cxxcimod.BinaryPolynomialModel_tuple3
            elif len(ind) == 4:
                base = cxxcimod.BinaryPolynomialModel_tuple4
            else:
                raise TypeError("invalid length of tuple")
        else:
            raise TypeError("invalid types of linear and quadratic")

    # now define class
    class BinaryPolynomialModel(base):
        """Represents Binary polynomial model.
           Note that the indices are converted to the integers internally. 
           The dictionaries between indices and integers are self.ind_to_num (indices -> integers) and self.num_to_ind (integers -> indices).
           Indices are listed in self._indices.
        Attributes:
            var_type (cimod.VariableType): variable type SPIN or BINARY
            polynomial (dict): represents polynomial term including linear term
            adj (dict): represents adjacency
            indices (list): labels of each variables sorted by results variables
            ind_to_num (list): map which specifies where the index is in self._indices
        """
        def __init__(self, polynomial, var_type=dimod.SPIN, **kwargs):
            super().__init__(polynomial, to_cxxcimod(var_type))
        
        @property
        def polynomial(self):
            return self.get_polynomial()
    
        @property
        def adj(self):
            return self.get_adjacency()
            
        @property
        def variables(self):
            return self.get_variables()
      
        @property
        def vartype(self):
            vartype = super().get_vartype()
            if vartype == cxxcimod.Vartype.SPIN:
                return dimod.SPIN
            else:
                return dimod.BINARY
            
    return BinaryPolynomialModel

def BinaryPolynomialModel(polynomial, var_type=dimod.SPIN, **kwargs):
    Model = make_BinaryPolynomialModel(polynomial)
    return Model(polynomial, var_type, **kwargs)


