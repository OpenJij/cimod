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

import cxxcimod
import cimod
from cimod.vartype import to_cxxcimod
import dimod
import numpy as np

def get_cxxcimod_class(linear, quadratic, sparse):
    # select base class
    index = set()
    base = None

    if linear != {}:
        index.add(next(iter(linear)))
    if quadratic != {}:
        index.add(next(iter(quadratic))[0])
        index.add(next(iter(quadratic))[1])

    if len(set(type(i) for i in index)) != 1:
        raise TypeError("invalid types of linear and quadratic")
    else:

        ind = next(iter(index))

        if sparse:
            if isinstance(ind, int):
                base = cxxcimod.BinaryQuadraticModel_Sparse
            elif isinstance(ind, str):
                base = cxxcimod.BinaryQuadraticModel_str_Sparse
            elif isinstance(ind, tuple):
                if len(ind) == 2:
                    base = cxxcimod.BinaryQuadraticModel_tuple2_Sparse
                elif len(ind) == 3:
                    base = cxxcimod.BinaryQuadraticModel_tuple3_Sparse
                elif len(ind) == 4:
                    base = cxxcimod.BinaryQuadraticModel_tuple4_Sparse
                else:
                    raise TypeError("invalid length of tuple")
            else:
                raise TypeError("invalid types of linear and quadratic")
        else:
            if isinstance(ind, int):
                base = cxxcimod.BinaryQuadraticModel_Dense
            elif isinstance(ind, str):
                base = cxxcimod.BinaryQuadraticModel_str_Dense
            elif isinstance(ind, tuple):
                if len(ind) == 2:
                    base = cxxcimod.BinaryQuadraticModel_tuple2_Dense
                elif len(ind) == 3:                             
                    base = cxxcimod.BinaryQuadraticModel_tuple3_Dense
                elif len(ind) == 4:
                    base = cxxcimod.BinaryQuadraticModel_tuple4_Dense
                else:
                    raise TypeError("invalid length of tuple")
            else:
                raise TypeError("invalid types of linear and quadratic")

    return base

def make_BinaryQuadraticModel(linear, quadratic, sparse):
    """BinaryQuadraticModel factory. 
       Generate BinaryQuadraticModel class with the base class specified by the arguments linear and quadratic
    Args:
        linear (dict): linear bias
        quadratic (dict): quadratic bias
        sparse (bool): if true, the inner data will be a sparse matrix, otherwise the data will be a dense matrix
    Returns:
        generated BinaryQuadraticModel class
    """

    Base = get_cxxcimod_class(linear, quadratic, sparse)

    # now define class
    class BinaryQuadraticModel(Base):
        """Represents Binary quadratic model. 
           Note that the indices are converted to the integers internally. 
           The dictionaries between indices and integers are self.ind_to_num (indices -> integers) and self.num_to_ind (integers -> indices).
           Indices are listed in self._indices.
        Attributes:
            vartype (cimod.VariableType): variable type SPIN or BINARY
            linear (dict): represents linear term
            quadratic (dict): represents quadratic term
            offset (float): represents constant energy term when convert to SPIN from BINARY
            num_variables (int): represents constant energy term when convert to SPIN from BINARY
            variables (list): represents constant energy term when convert to SPIN from BINARY
        """
        
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        @property
        def vartype(self):
            vartype = super().get_vartype()
            if vartype == cxxcimod.Vartype.SPIN:
                return dimod.SPIN
            else:
                #BINARY
                return dimod.BINARY

        @property
        def linear(self):
            return self.get_linear()
    
        @property
        def quadratic(self):
            return self.get_quadratic()
    
        @property
        def offset(self):
            return self.get_offset()

        @property
        def num_variables(self):
            return self.get_num_variables()
    
        @property
        def variables(self):
            return self.get_variables()

        def empty(self, vartype):
            return self.__class__(super().empty(to_cxxcimod(vartype)))

        def change_vartype(self, vartype, inplace=None):
            if inplace is None:
                super().change_vartype(to_cxxcimod(vartype))
            else:
                return self.__class__(super().change_vartype(to_cxxcimod(vartype), inplace))

        def __str__(self):
            return f"BinaryQuadraticModel({self.linear}, {self.quadratic}, {self.offset}, {self.vartype})"

        def __repr__(self):
            return f"BinaryQuadraticModel({self.linear}, {self.quadratic}, {self.offset}, {self.vartype})"

        @classmethod
        def from_numpy_matrix(cls, mat, variables: list, offset=0.0, vartype='BINARY', fix_format=True, **kwargs):
            shape = np.shape(mat)
            if not (len(shape) == 2 and shape[0] == shape[1]):
                raise TypeError("numpy matrix has to be a square matrix")

            num_variables = shape[0] - 1

            return cls(mat, variables, offset, to_cxxcimod(vartype), fix_format, **kwargs)

        @classmethod
        def from_qubo(cls, Q, offset=0.0, **kwargs):
            cxxbqm = Base.from_qubo(Q, offset)
            return cls(cxxbqm, **kwargs)

        @classmethod
        def from_ising(cls, linear, quadratic, offset=0.0, **kwargs):
            cxxbqm = Base.from_ising(linear, quadratic, offset)
            return cls(cxxbqm, **kwargs)

        @classmethod
        def from_serializable(cls, obj, **kwargs):
            cxxbqm = Base.from_serializable(obj)
            return cls(cxxbqm, **kwargs)



    return BinaryQuadraticModel

# for JSON
def make_BinaryQuadraticModel_from_JSON(obj):
    label = obj['variable_labels'][0]
    if isinstance(label, list):
        #convert to tuple
        label = tuple(label)

    mock_linear = {label:1.0}

    if obj['version']['bqm_schema'] == '3.0.0-dense':
        sparse = False
    elif obj['version']['bqm_schema'] == '3.0.0':
        sparse = True 
    else:
        raise TypeError("Invalid bqm_schema")

    return make_BinaryQuadraticModel(mock_linear, {}, sparse)


def BinaryQuadraticModel(linear, quadratic, *args, **kwargs):
    Model = make_BinaryQuadraticModel(linear, quadratic, kwargs.pop('sparse', False))

    # offset and vartype
    if len(args) == 2:
        [offset, vartype] = args
        return Model(linear, quadratic, offset, to_cxxcimod(vartype))
    elif len(args) == 1 and 'vartype' in kwargs:
        [offset] = args
        vartype = kwargs['vartype']
        return Model(linear, quadratic, offset, to_cxxcimod(vartype))
    elif len(args) == 1:
        [vartype] = args
        return Model(linear, quadratic, 0.0, to_cxxcimod(vartype))
    elif len(args) == 0 and 'vartype' in kwargs:
        vartype = kwargs['vartype']
        return Model(linear, quadratic, 0.0, to_cxxcimod(vartype))
    else:
        raise TypeError("invalid args for BinaryQuadraticModel. please check arguments")

def bqm_from_numpy_matrix(mat, variables: list=None, offset=0.0, vartype='BINARY', **kwargs):
    if variables is None:
        # generate array
        num_variables = mat.shape[0]
        variables = list(range(num_variables))

    return make_BinaryQuadraticModel({variables[0]: 1.0}, {}, kwargs.pop('sparse', False)).from_numpy_matrix(mat, variables, offset, to_cxxcimod(vartype), True, **kwargs)


BinaryQuadraticModel.from_numpy_matrix = bqm_from_numpy_matrix

BinaryQuadraticModel.from_qubo = \
lambda Q, offset=0.0, **kwargs: make_BinaryQuadraticModel({}, Q, kwargs.pop('sparse', False)).from_qubo(Q, offset, **kwargs)

BinaryQuadraticModel.from_qubo = \
lambda Q, offset=0.0, **kwargs: make_BinaryQuadraticModel({}, Q, kwargs.pop('sparse', False)).from_qubo(Q, offset, **kwargs)

BinaryQuadraticModel.from_ising = \
lambda linear, quadratic, offset=0.0, **kwargs: make_BinaryQuadraticModel(linear, quadratic, kwargs.pop('sparse', False)).from_ising(linear, quadratic, offset, **kwargs)

BinaryQuadraticModel.from_serializable = \
lambda obj, **kwargs: make_BinaryQuadraticModel_from_JSON(obj).from_serializable(obj, **kwargs)



