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
import dimod
import numpy as np
from cimod.vartype import to_cxxcimod
from cimod.utils.decolator import methoddispatch
from functools import singledispatch
from enum import Enum

class Polynomial:
    def __init__(self, bpm):
        self._bpm = bpm

    def __repr__(self):
        return  str(self._bpm.get_polynomial())

class Variables:
    def __init__(self, bpm):
        self._bpm = bpm

    def __repr__(self):
        return  str(self._bpm.get_variables())

class IndexType(Enum):
    INT = 1
    STRING = 2
    INT_TUPLE_2 = 3
    INT_TUPLE_3 = 4
    INT_TUPLE_4 = 5

def make_BinaryPolynomialModel(polynomial, index_type = None, tuple_size = 0):
    """BinaryPolynomialModel factory.
       Generate BinaryPolynomialModel class with the base class specified by the arguments linear and quadratic
    Args:
        polynomial (dict): polynomial bias including linear bias
    Returns:
        generated BinaryPolynomialModel class
    """

    def base_selector(index_type, index):
        if index_type == int or index_type == None:
            return cxxcimod.BinaryPolynomialModel, IndexType.INT
        elif index_type == str:
            return cxxcimod.BinaryPolynomialModel_str, IndexType.STRING
        elif index_type == tuple:
            if len(index) == 2:
                return cxxcimod.BinaryPolynomialModel_tuple2, IndexType.INT_TUPLE_2
            elif len(index) == 3:
                return cxxcimod.BinaryPolynomialModel_tuple3, IndexType.INT_TUPLE_3
            elif len(index) == 4:
                return cxxcimod.BinaryPolynomialModel_tuple4, IndexType.INT_TUPLE_4
            raise TypeError("Invalid length of tuple")
        else:
            raise TypeError("Invalid types of polynomial")

    index = set()
    base  = None

    if polynomial != {}:
        if len(polynomial) == 1 and tuple() in polynomial:
            base, base_type = base_selector(index_type, [1 for _ in range(min(tuple_size, 4))])
        elif len(polynomial) > 1 and next(iter(polynomial)) == tuple():
            iter_poly = iter(polynomial)
            next(iter_poly)
            second_tuple = next(iter_poly)
            if len(set(type(i) for i in second_tuple)) != 1:
                raise TypeError("Invalid types of polynomial")
            base, base_type = base_selector(type(second_tuple[0]), second_tuple[0])
        else:
            iter_poly = iter(polynomial)
            first_tuple = next(iter_poly)
            if len(set(type(i) for i in first_tuple)) != 1:
                raise TypeError("Invalid types of polynomial")
            base, base_type = base_selector(type(first_tuple[0]), first_tuple[0])
    else:
        base, base_type = base_selector(index_type, [1 for _ in range(min(tuple_size, 4))])

    # now define class
    class BinaryPolynomialModel(base):
        """Represents Binary polynomial model.
           Note that the indices are converted to the integers internally. 
           The dictionaries between indices and integers are self.ind_to_num (indices -> integers) and self.num_to_ind (integers -> indices).
           Indices are listed in self._indices.
        Attributes:
            vartype (cimod.VariableType): variable type SPIN or BINARY
            polynomial (dict): represents polynomial term including linear term
            variables (list): labels of each variables sorted by results variables
        """
        @methoddispatch
        def __init__(self, polynomial: dict, vartype):
            self.index_type = base_type
            super().__init__(polynomial, to_cxxcimod(vartype))
                
        @__init__.register
        def __init__from_list(self, keys: list, values: list, vartype):
            self.index_type = base_type
            super().__init__(keys, values, to_cxxcimod(vartype))

        def _model_selector(self):
            if self.index_type == IndexType.INT:
                return make_BinaryPolynomialModel({}, int)
            elif self.index_type == IndexType.STRING:
                return make_BinaryPolynomialModel({}, str)
            elif self.index_type == IndexType.INT_TUPLE_2:
                return make_BinaryPolynomialModel({}, tuple, 2)
            elif self.index_type == IndexType.INT_TUPLE_3:
                return make_BinaryPolynomialModel({}, tuple, 3)
            elif self.index_type == IndexType.INT_TUPLE_4:
                return make_BinaryPolynomialModel({}, tuple, 4)
            else:
                raise TypeError("invalid types of polynomial")

        @property
        def polynomial(self):
            return Polynomial(self)
            
        @property
        def variables(self):
            return Variables(self)

        @property
        def degree(self):
            return super().get_degree()

        @property
        def num_interactions(self):
            return super().get_num_interactions()

        @property
        def num_variables(self):
            return super().get_num_variables()

        @property
        def vartype(self):
            vartype = super().get_vartype()
            if vartype == cxxcimod.Vartype.SPIN:
                return dimod.SPIN
            elif vartype == cxxcimod.Vartype.BINARY:
                return dimod.BINARY
            else:
                raise Exception("Unknown vartype detected")

        def print(self):
            print("[BinaryPolynomialModel]")
            print("index_type =", self.index_type)
            print("Variables =")
            print(self.get_variables())
            print("polynomial =")
            print(self.get_polynomial())
            print("vartype =", self.get_vartype())
            print("num_variables =", self.get_num_variables())
            print("num_interactions =", self.get_num_interactions())

        def empty(self, vartype):
            Model = self._model_selector()
            return Model({}, to_cxxcimod(vartype))

        def add_interaction(self, key: list, value, vartype = cxxcimod.Vartype.NONE):
            return super().add_interaction(key, value, to_cxxcimod(vartype))

        def get_polynomial(self, *args, **kwargs):
            if kwargs != {}:
                return super().get_polynomial(*args, **kwargs)

            if args == tuple():
                return super().get_polynomial()
            elif args[0] == tuple() or args[0] == []:
                return super().get_polynomial(())
            elif self.index_type == IndexType.INT or self.index_type == IndexType.STRING:
                if type(args[0]) == int or type(args[0]) == str:
                    return super().get_polynomial(args)
                else:
                    return super().get_polynomial(*args, **kwargs)
            else:
                if type(args[0][0]) == int or type(args[0][0]) == str:
                    return super().get_polynomial(args)
                else:
                    return super().get_polynomial(*args, **kwargs)

        def remove_interaction(self, *args, **kwargs):
            if kwargs != {}:
                return super().remove_interaction(*args, **kwargs)

            if args == tuple():
                raise TypeError("Enter an argument.")

            elif args[0] == tuple() or args[0] == []:
                return super().remove_interaction(())

            elif self.index_type == IndexType.INT or self.index_type == IndexType.STRING:
                if type(args[0]) == int or type(args[0]) == str:
                    return super().remove_interaction(args)
                else:
                    return super().remove_interaction(*args, **kwargs)
            else:
                if type(args[0][0]) == int or type(args[0][0]) == str:
                    return super().remove_interaction(args)
                else:
                    return super().remove_interaction(*args, **kwargs)

        def remove_interactions_from(self, *args, **kwargs):
            if kwargs != {}:
                return super().remove_interactions_from(*args, **kwargs)

            if self.index_type == IndexType.INT or self.index_type == IndexType.STRING:
                if type(args[0][0]) == int or type(args[0][0]) == str:
                    return super().remove_interactions_from(args)
                else:
                    return super().remove_interactions_from(*args, **kwargs)
            else:
                if type(args[0][0][0]) == int or type(args[0][0][0]) == str:
                    return super().remove_interactions_from(args)
                else:
                    return super().remove_interactions_from(*args, **kwargs)

        @methoddispatch
        def add_interactions_from(self, polynomial: dict, vartype = cxxcimod.Vartype.NONE):
            return super().add_interactions_from(polynomial, to_cxxcimod(vartype))

        @add_interactions_from.register
        def _add_interactions_from_from_list(self, keys: list, values: list, vartype = cxxcimod.Vartype.NONE):
            return super().add_interactions_from(keys, values, to_cxxcimod(vartype))

        def change_vartype(self, vartype, inplace = None):
            vartype = to_cxxcimod(vartype)
            if inplace == None or inplace == True:
                return super().change_vartype(vartype)
            elif inplace == False:
                Model = self._model_selector()
                if to_cxxcimod(self.vartype) == vartype:
                    return Model(self._get_keys(), self._get_values(), vartype)
                else:
                    if vartype == cxxcimod.SPIN:
                        return Model(self.to_hising(), vartype)
                    elif vartype == cxxcimod.BINARY:
                        return Model(self.to_hubo(), vartype)
                    else:
                        raise Exception("Unknown vartype error")
            else:
                raise TypeError("Invalid inplace value")

        @classmethod
        def from_hising(cls, *args, **kwargs):
            return cls(*args, **kwargs, vartype = dimod.SPIN)
            
        @classmethod
        def from_hubo(cls, *args, **kwargs):
            return cls(*args, **kwargs, vartype = dimod.BINARY)
        
        @classmethod
        def from_serializable(cls, obj):
            if(obj["type"] != "BinaryPolynomialModel"):
                raise Exception("Type must be \"BinaryPolynomialModel\".\n")
            return cls(obj['poly_key_list'], obj['poly_value_list'], obj['vartype'])

        def __repr__(self):
            ss = "BinaryPolynomialModel(" + str(self.get_polynomial()) + ", " + str(self.get_vartype()) + ")"
            return ss

    return BinaryPolynomialModel

def make_BinaryPolynomialModel_from_JSON(obj):
    label = obj["poly_key_list"][0][0]
    if isinstance(label, list):
        label = tuple(label)
    mock_polynomial = {(label,):1.0}
    return make_BinaryPolynomialModel(mock_polynomial)

@singledispatch
def BinaryPolynomialModel(polynomial: dict, vartype):
    Model = make_BinaryPolynomialModel(polynomial)
    return Model(polynomial, vartype)

@BinaryPolynomialModel.register
def _BinaryPolynomialModel_from_list(keys: list, values: list, vartype):
    label = keys[0][0]
    if isinstance(label, list):
        label = tuple(label)
    mock_polynomial = {(label,):1.0}
    Model = make_BinaryPolynomialModel(mock_polynomial)
    return Model(keys, values, vartype)

@singledispatch
def make_BinaryPolynomialModel_from_hising(polynomial: dict):
    return make_BinaryPolynomialModel(polynomial).from_hising(polynomial)

@make_BinaryPolynomialModel_from_hising.register
def _make_BinaryPolynomialModel_from_hising_from_list(keys: list, values: list):
    label = keys[0][0]
    if isinstance(label, list):
        label = tuple(label)
    mock_polynomial = {(label,):1.0}
    return make_BinaryPolynomialModel(mock_polynomial).from_hising(keys, values)

@singledispatch
def make_BinaryPolynomialModel_from_hubo(polynomial: dict):
    return make_BinaryPolynomialModel(polynomial).from_hubo(polynomial)

def _make_BinaryPolynomialModel_from_hubo_from_list(keys: list, values: list):
    label = keys[0][0]
    if isinstance(label, list):
        label = tuple(label)
    mock_polynomial = {(label,):1.0}
    return make_BinaryPolynomialModel(mock_polynomial).from_hubo(keys, values)

#classmethods
BinaryPolynomialModel.from_serializable = lambda obj: make_BinaryPolynomialModel_from_JSON(obj).from_serializable(obj)
BinaryPolynomialModel.from_hising       = lambda *args, **kwargs: make_BinaryPolynomialModel_from_hising(*args, **kwargs)
BinaryPolynomialModel.from_hubo         = lambda *args, **kwargs: make_BinaryPolynomialModel_from_hubo(*args, **kwargs)

