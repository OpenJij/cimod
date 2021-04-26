# Copyright 2021 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
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
        raise TypeError("invalid types of polynomial")
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
            raise TypeError("invalid types of polynomial")

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
            self._init_process()

        def _init_process(self):
            # indices
            self._re_calculate_indices = True
            self._indices    = None
            self._ind_to_num = None

        
        @property
        def indices(self):
            ind, _ = self.update_indices()
            return ind

        def update_indices(self):
            """calculate self._indices and self.ind_to_num
            Returns:
                self._indices and self._ind_to_num
            """
            if self._re_calculate_indices is True:
                self._indices = self._generate_indices()
                self._ind_to_num = {ind:num for num,ind in enumerate(self._indices)}

                self._re_calculate_indices = False

            return self._indices, self._ind_to_num
        
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

        def energy(self, sample, convert_sample=False):
            # convert samples to SPIN or BINARY
            if convert_sample:
                for k in sample.keys():
                    if sample[k] == -1 and self.vartype == dimod.BINARY:
                        sample[k] = 0
                    if sample[k] == 0  and self.vartype == dimod.SPIN:
                        sample[k] = -1
            return super().energy(sample)

        def energies(self, samples_like, convert_sample=False, **kwargs):
            en_vec = []
            for elem in samples_like:
                en_vec.append(self.energy(elem, convert_sample, **kwargs))
    
            return en_vec

        @classmethod
        def from_ising(cls, polynomial, **kwargs):
            return cls(polynomial, var_type=dimod.SPIN, **kwargs)
            
        @classmethod
        def from_pubo(cls, polynomial, **kwargs):
            return cls(polynomial, var_type=dimod.BINARY, **kwargs)
        
        @classmethod
        def from_serializable(cls, obj):
    
            variable_labels = [tuple(elem) if type(elem) == list else elem for elem in obj['variable_labels']]

            #convert to polynomial biases
            zipped_obj = zip(obj["polynomial_interactions"], obj["polynomial_biases"])

            polynomial = {}
            for elem in zipped_obj:
                temp = ()
                if (type(elem[0][0]) == int) or (type(elem[0][0]) == str):
                    for i in elem[0]:
                        temp += (variable_labels[i],)
                else:
                    for i in elem[0]:
                        temp += (tuple(variable_labels[i]),)
                polynomial.update({(temp):elem[1]})
                #polynomial |= {(temp):elem[1]}
            
            # set vartype
            vartype = cimod.SPIN if obj['variable_type'] == 'SPIN' else cimod.BINARY

            return cls(polynomial, vartype)


    return BinaryPolynomialModel

# for JSON
def make_BinaryPolynomialModel_from_JSON(obj):
    label = obj['variable_labels'][0]
    if isinstance(label, list):
        #convert to tuple
        label = tuple(label)
    mock_linear = {(label,):1.0}
    return make_BinaryPolynomialModel(mock_linear)

def BinaryPolynomialModel(polynomial, var_type=dimod.SPIN, **kwargs):
    Model = make_BinaryPolynomialModel(polynomial)
    return Model(polynomial, var_type, **kwargs)

#classmethods
BinaryPolynomialModel.from_serializable = lambda obj: make_BinaryPolynomialModel_from_JSON(obj).from_serializable(obj)
BinaryPolynomialModel.from_ising = \
lambda polynomial, **kwargs: make_BinaryPolynomialModel(polynomial).from_ising(polynomial, **kwargs)
BinaryPolynomialModel.from_pubo = \
lambda polynomial, **kwargs: make_BinaryPolynomialModel(polynomial).from_pubo(polynomial, **kwargs)
'''