# Copyright 2022 Jij Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

import cimod.cxxcimod as cxxcimod
import dimod

from cimod.vartype import to_cxxcimod


class Polynomial:
    def __init__(self, bpm):
        self.__bpm = bpm

    def items(self, *args, **kwargs):
        return self.__bpm.get_polynomial().items(*args, **kwargs)

    def keys(self, *args, **kwargs):
        return self.__bpm.get_polynomial().keys(*args, **kwargs)

    def values(self, *args, **kwargs):
        return self.__bpm.get_polynomial().values(*args, **kwargs)

    def copy(self):
        return self.__bpm.get_polynomial()

    def fromkeys(self, keys=None):
        return dict.fromkeys(self.__bpm.get_polynomial(), keys)

    def get(self, arg1, arg2=None):
        val = self.__bpm.get_polynomial(arg1)
        if val == 0.0:
            if arg2 is not None:
                return arg2
        else:
            return val

    def __len__(self):
        return self.__bpm.num_interactions

    def __repr__(self):
        return str(self.__bpm.get_polynomial())

    """
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
    """


def make_BinaryPolynomialModel(polynomial, index_type=None, tuple_size=0):
    """BinaryPolynomialModel factory.
       Generate BinaryPolynomialModel class with the base class specified by the arguments linear and quadratic
    Args:
        polynomial (dict): polynomial bias including linear bias
    Returns:
        generated BinaryPolynomialModel class
    """

    def base_selector(index_type, index):
        if index_type == int or index_type is None:
            return cxxcimod.BinaryPolynomialModel, "IndexType.INT"
        elif index_type == str:
            return cxxcimod.BinaryPolynomialModel_str, "IndexType.STRING"
        elif index_type == tuple:
            if len(index) == 2:
                return cxxcimod.BinaryPolynomialModel_tuple2, "IndexType.INT_TUPLE_2"
            elif len(index) == 3:
                return cxxcimod.BinaryPolynomialModel_tuple3, "IndexType.INT_TUPLE_3"
            elif len(index) == 4:
                return cxxcimod.BinaryPolynomialModel_tuple4, "IndexType.INT_TUPLE_4"
            raise TypeError("Invalid length of tuple")
        else:
            raise TypeError("Invalid types of polynomial")

    base = None

    if polynomial != {}:
        if len(polynomial) == 1 and tuple() in polynomial:
            base, base_type = base_selector(
                index_type, [1 for _ in range(min(tuple_size, 4))]
            )
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
        base, base_type = base_selector(
            index_type, [1 for _ in range(min(tuple_size, 4))]
        )

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

        def __init__(self, *args, **kwargs):
            self.index_type = base_type
            self.model_type = "cimod.BinaryPolynomialModel"
            super().__init__(*args, **kwargs)

        def _model_selector(self):
            if self.index_type == "IndexType.INT":
                return make_BinaryPolynomialModel({}, int)
            elif self.index_type == "IndexType.STRING":
                return make_BinaryPolynomialModel({}, str)
            elif self.index_type == "IndexType.INT_TUPLE_2":
                return make_BinaryPolynomialModel({}, tuple, 2)
            elif self.index_type == "IndexType.INT_TUPLE_3":
                return make_BinaryPolynomialModel({}, tuple, 3)
            elif self.index_type == "IndexType.INT_TUPLE_4":
                return make_BinaryPolynomialModel({}, tuple, 4)
            else:
                raise TypeError("invalid types of polynomial")

        @property
        def polynomial(self):
            return Polynomial(self)

        @property
        def variables(self):
            return super().get_variables()

        # This function is depricated
        @property
        def indices(self):
            return super().indices()

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

        def add_interaction(self, key: list, value, vartype=cxxcimod.Vartype.NONE):
            return super().add_interaction(key, value, to_cxxcimod(vartype))

        def get_polynomial(self, *args, **kwargs):
            if kwargs != {}:
                return super().get_polynomial(*args, **kwargs)

            if args == tuple():
                return super().get_polynomial()
            elif args[0] == tuple() or args[0] == []:
                return super().get_polynomial(())
            elif (
                self.index_type == "IndexType.INT"
                or self.index_type == "IndexType.STRING"
            ):
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

            elif (
                self.index_type == "IndexType.INT"
                or self.index_type == "IndexType.STRING"
            ):
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

            if (
                self.index_type == "IndexType.INT"
                or self.index_type == "IndexType.STRING"
            ):
                if type(args[0][0]) == int or type(args[0][0]) == str:
                    return super().remove_interactions_from(args)
                else:
                    return super().remove_interactions_from(*args, **kwargs)
            else:
                if type(args[0][0][0]) == int or type(args[0][0][0]) == str:
                    return super().remove_interactions_from(args)
                else:
                    return super().remove_interactions_from(*args, **kwargs)

        def get_variables_to_integers(self, *args, **kwargs):
            obj = super().get_variables_to_integers(*args, **kwargs)
            if obj != -1:
                return obj

        def to_serializable(self):
            obj = super().to_serializable()
            obj["index_type"] = self.index_type
            return obj

        def add_interactions_from(self, *args, **kwargs):
            if kwargs == {}:
                if len(args) == 0:
                    raise TypeError("Invalid argument for this function")
                if len(args) == 1:
                    if isinstance(args[0], dict):
                        super().add_interactions_from(args[0])
                    else:
                        raise TypeError("Invalid argument for this function")
                elif len(args) == 2:
                    key_condition = isinstance(args[0], list) or isinstance(
                        args[0], tuple
                    )
                    val_condition = isinstance(args[1], list) or isinstance(
                        args[1], tuple
                    )
                    if isinstance(args[0], dict):
                        super().add_interactions_from(args[0], to_cxxcimod(args[1]))
                    elif key_condition and val_condition:
                        super().add_interactions_from(args[0], args[1])
                    else:
                        raise TypeError("Invalid argument for this function")
                elif len(args) == 3:
                    key_condition = isinstance(args[0], list) or isinstance(
                        args[0], tuple
                    )
                    val_condition = isinstance(args[1], list) or isinstance(
                        args[1], tuple
                    )
                    if key_condition and val_condition:
                        super().add_interactions_from(
                            args[0], args[1], to_cxxcimod(args[2])
                        )
                    else:
                        raise TypeError("Invalid argument for this function")
                else:
                    raise TypeError("Invalid argument for this function")
            else:
                if "keys" in kwargs and "values" in kwargs and "vartype" in kwargs:
                    if len(args) != 0:
                        raise TypeError("Invalid argument for this function")
                    super().add_interactions_from(
                        kwargs["keys"], kwargs["values"], to_cxxcimod(kwargs["vartype"])
                    )
                elif "values" in kwargs and "vartype" in kwargs:
                    if len(args) != 1:
                        raise TypeError("Invalid argument for this function")
                    key_condition = isinstance(args[0], list) or isinstance(
                        args[0], tuple
                    )
                    val_condition = isinstance(kwargs["values"], list) or isinstance(
                        kwargs["values"], tuple
                    )
                    if key_condition and val_condition:
                        super().add_interactions_from(
                            args[0], kwargs["values"], to_cxxcimod(kwargs["vartype"])
                        )
                    else:
                        raise TypeError("Invalid argument for this function")
                elif "polynomial" in kwargs and "vartype" in kwargs:
                    if len(args) != 0:
                        raise TypeError("Invalid argument for this function")
                    if isinstance(kwargs["polynomial"], dict):
                        super().add_interactions_from(
                            kwargs["polynomial"], to_cxxcimod(kwargs["vartype"])
                        )
                    else:
                        raise TypeError("Invalid argument for this function")
                elif "keys" in kwargs and "values" in kwargs:
                    if len(args) != 0:
                        raise TypeError("Invalid argument for this function")
                    key_condition = isinstance(kwargs["keys"], list) or isinstance(
                        kwargs["keys"], tuple
                    )
                    val_condition = isinstance(kwargs["values"], list) or isinstance(
                        kwargs["values"], tuple
                    )
                    if key_condition and val_condition:
                        super().add_interactions_from(kwargs["keys"], kwargs["values"])
                    else:
                        raise TypeError("Invalid argument for this function")
                elif "vartype" in kwargs:
                    if len(args) == 1:
                        if isinstance(args[0], dict):
                            super().add_interactions_from(
                                args[0], to_cxxcimod(kwargs["vartype"])
                            )
                        else:
                            raise TypeError("Invalid argument for this function")
                    elif len(args) == 2:
                        key_condition = isinstance(args[0], list) or isinstance(
                            args[0], tuple
                        )
                        val_condition = isinstance(args[1], list) or isinstance(
                            args[1], tuple
                        )
                        if key_condition and val_condition:
                            super().add_interactions_from(
                                args[0], args[1], to_cxxcimod(kwargs["vartype"])
                            )
                        else:
                            raise TypeError("Invalid argument for this function")
                elif "values" in kwargs:
                    if len(args) != 1:
                        raise TypeError("Invalid argument for this function")
                    key_condition = isinstance(args[0], list) or isinstance(
                        args[0], tuple
                    )
                    val_condition = isinstance(kwargs["values"], list) or isinstance(
                        kwargs["values"], tuple
                    )
                    if key_condition and val_condition:
                        super().add_interactions_from(args[0], kwargs["values"])
                    else:
                        raise TypeError("Invalid argument for this function")
                elif "polynomial" in kwargs:
                    if len(args) != 0:
                        raise TypeError("Invalid argument for this function")
                    if isinstance(kwargs["polynomial"], dict):
                        super().add_interactions_from(kwargs["polynomial"])
                    else:
                        raise TypeError("Invalid argument for this function")
                else:
                    raise TypeError("Invalid argument for this function")

        def change_vartype(self, vartype, inplace=None):
            vartype = to_cxxcimod(vartype)
            if inplace is None or inplace:
                return super().change_vartype(vartype)
            elif inplace == False:
                Model = self._model_selector()
                if to_cxxcimod(self.vartype) == vartype:
                    return Model(self.get_key_list(), self.get_value_list(), vartype)
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
            return cls(*args, **kwargs, vartype=cxxcimod.SPIN)

        @classmethod
        def from_hubo(cls, *args, **kwargs):
            return cls(*args, **kwargs, vartype=cxxcimod.BINARY)

        @classmethod
        def from_serializable(cls, obj):
            if obj["type"] != "BinaryPolynomialModel":
                raise Exception('Type must be "BinaryPolynomialModel"')
            return cls(
                obj["variables"],
                obj["poly_key_distance_list"],
                obj["poly_value_list"],
                to_cxxcimod(obj["vartype"]),
            )

        def __repr__(self):
            ss = (
                "BinaryPolynomialModel("
                + str(self.get_polynomial())
                + ", "
                + str(self.get_vartype())
                + ")"
            )
            return ss

    return BinaryPolynomialModel


def make_BinaryPolynomialModel_from_JSON(obj):
    if obj["type"] != "BinaryPolynomialModel":
        raise Exception('Type must be "BinaryPolynomialModel"')
    mock_polynomial = {}
    if obj["index_type"] == "IndexType.INT":
        mock_polynomial = {(0, 1): 1}
    elif obj["index_type"] == "IndexType.STRING":
        mock_polynomial = {("a", "b"): 1}
    elif obj["index_type"] == "IndexType.INT_TUPLE_2":
        mock_polynomial = {((0, 1), (1, 2)): 1}
    elif obj["index_type"] == "IndexType.INT_TUPLE_3":
        mock_polynomial = {((0, 1, 2), (1, 2, 3)): 1}
    elif obj["index_type"] == "IndexType.INT_TUPLE_4":
        mock_polynomial = {((0, 1, 2, 3), (1, 2, 3, 4)): 1}
    else:
        raise TypeError("Invalid types of polynomial")
    return make_BinaryPolynomialModel(mock_polynomial)


def BinaryPolynomialModel(*args, **kwargs):
    if kwargs == {}:
        if len(args) <= 1:
            raise TypeError("Invalid argument for this function")
        elif len(args) == 2:
            if isinstance(args[0], dict):
                return _BinaryPolynomialModel_from_dict(args[0], to_cxxcimod(args[1]))
            else:
                raise TypeError("Invalid argument for this function")
        elif len(args) == 3:
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(args[1], list) or isinstance(args[1], tuple)
            if key_condition and val_condition:
                return _BinaryPolynomialModel_from_list(
                    args[0], args[1], to_cxxcimod(args[2])
                )
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")
    else:
        if "keys" in kwargs and "values" in kwargs and "vartype" in kwargs:
            key_condition = isinstance(kwargs["keys"], list) or isinstance(
                kwargs["keys"], tuple
            )
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _BinaryPolynomialModel_from_list(
                    kwargs["keys"], kwargs["values"], to_cxxcimod(kwargs["vartype"])
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "polynomial" in kwargs and "vartype" in kwargs:
            if isinstance(kwargs["polynomial"], dict):
                return _BinaryPolynomialModel_from_dict(
                    kwargs["polynomial"], to_cxxcimod(kwargs["vartype"])
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "values" in kwargs and "vartype" in kwargs:
            if len(args) != 1:
                raise TypeError("Invalid argument for this function")
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _BinaryPolynomialModel_from_list(
                    args[0], kwargs["values"], to_cxxcimod(kwargs["vartype"])
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "vartype" in kwargs:
            if len(args) == 1:
                if isinstance(args[0], dict):
                    return _BinaryPolynomialModel_from_dict(
                        args[0], to_cxxcimod(kwargs["vartype"])
                    )
                else:
                    raise TypeError("Invalid argument for this function")
            elif len(args) == 2:
                key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
                val_condition = isinstance(args[1], list) or isinstance(args[1], tuple)
                if key_condition and val_condition:
                    return _BinaryPolynomialModel_from_list(
                        args[0], args[1], to_cxxcimod(kwargs["vartype"])
                    )
                else:
                    raise TypeError("Invalid argument for this function")
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")


def _BinaryPolynomialModel_from_dict(polynomial: dict, vartype):
    Model = make_BinaryPolynomialModel(polynomial)
    return Model(polynomial, to_cxxcimod(vartype))


def _BinaryPolynomialModel_from_list(keys: list, values: list, vartype):
    if len(keys) == 0:
        Model = make_BinaryPolynomialModel({})
        return Model(keys, values, to_cxxcimod(vartype))
    i = 0
    label = None
    while i < len(keys):
        if len(keys[i]) > 0:
            label = keys[i][0]
            break
        i += 1
    if label is None:
        Model = make_BinaryPolynomialModel({(): 1.0})
        return Model(keys, values, to_cxxcimod(vartype))
    else:
        if isinstance(label, list):
            label = tuple(label)
        mock_polynomial = {(label,): 1.0}
        Model = make_BinaryPolynomialModel(mock_polynomial)
        return Model(keys, values, to_cxxcimod(vartype))


def make_BinaryPolynomialModel_from_hising(*args, **kwargs):
    if kwargs == {}:
        if len(args) == 0:
            raise TypeError("Invalid argument for this function")
        elif len(args) == 1:
            if isinstance(args[0], dict):
                return _make_BinaryPolynomialModel_from_hising_from_dict(args[0])
            else:
                raise TypeError("Invalid argument for this function")
        elif len(args) == 2:
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(args[1], list) or isinstance(args[1], tuple)
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hising_from_list(
                    args[0], args[1]
                )
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")
    else:
        if "keys" in kwargs and "values" in kwargs:
            key_condition = isinstance(kwargs["keys"], list) or isinstance(
                kwargs["keys"], tuple
            )
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hising_from_list(
                    kwargs["keys"], kwargs["values"]
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "values" in kwargs:
            if len(args) != 1:
                raise TypeError("Invalid argument for this function")
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hising_from_list(
                    args[0], kwargs["values"]
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "polynomial" in kwargs:
            if len(args) != 0:
                raise TypeError("Invalid argument for this function")
            if isinstance(kwargs["polynomial"], dict):
                _make_BinaryPolynomialModel_from_hising_from_dict(kwargs["polynomial"])
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")


def _make_BinaryPolynomialModel_from_hising_from_dict(polynomial: dict):
    return make_BinaryPolynomialModel(polynomial).from_hising(polynomial)


def _make_BinaryPolynomialModel_from_hising_from_list(keys: list, values: list):
    if len(keys) == 0:
        return make_BinaryPolynomialModel({}).from_hising(keys, values)

    i = 0
    label = None
    while i < len(keys):
        if len(keys[i]) > 0:
            label = keys[i][0]
            break
        i += 1

    if label is None:
        return make_BinaryPolynomialModel({(): 1.0}).from_hising(keys, values)
    else:
        if isinstance(label, list):
            label = tuple(label)
        mock_polynomial = {(label,): 1.0}
        return make_BinaryPolynomialModel(mock_polynomial).from_hising(keys, values)


def make_BinaryPolynomialModel_from_hubo(*args, **kwargs):
    if kwargs == {}:
        if len(args) == 0:
            raise TypeError("Invalid argument for this function")
        elif len(args) == 1:
            if isinstance(args[0], dict):
                return _make_BinaryPolynomialModel_from_hubo_from_dict(args[0])
            else:
                raise TypeError("Invalid argument for this function")
        elif len(args) == 2:
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(args[1], list) or isinstance(args[1], tuple)
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hubo_from_list(args[0], args[1])
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")
    else:
        if "keys" in kwargs and "values" in kwargs:
            key_condition = isinstance(kwargs["keys"], list) or isinstance(
                kwargs["keys"], tuple
            )
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hubo_from_list(
                    kwargs["keys"], kwargs["values"]
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "values" in kwargs:
            if len(args) != 1:
                raise TypeError("Invalid argument for this function")
            key_condition = isinstance(args[0], list) or isinstance(args[0], tuple)
            val_condition = isinstance(kwargs["values"], list) or isinstance(
                kwargs["values"], tuple
            )
            if key_condition and val_condition:
                return _make_BinaryPolynomialModel_from_hubo_from_list(
                    args[0], kwargs["values"]
                )
            else:
                raise TypeError("Invalid argument for this function")
        elif "polynomial" in kwargs:
            if len(args) != 0:
                raise TypeError("Invalid argument for this function")
            if isinstance(kwargs["polynomial"], dict):
                _make_BinaryPolynomialModel_from_hubo_from_dict(kwargs["polynomial"])
            else:
                raise TypeError("Invalid argument for this function")
        else:
            raise TypeError("Invalid argument for this function")


def _make_BinaryPolynomialModel_from_hubo_from_dict(polynomial: dict):
    return make_BinaryPolynomialModel(polynomial).from_hubo(polynomial)


def _make_BinaryPolynomialModel_from_hubo_from_list(keys: list, values: list):
    if len(keys) == 0:
        return make_BinaryPolynomialModel({}).from_hubo(keys, values)

    i = 0
    label = None
    while i < len(keys):
        if len(keys[i]) > 0:
            label = keys[i][0]
            break
        i += 1

    if label is None:
        return make_BinaryPolynomialModel({(): 1.0}).from_hubo(keys, values)
    else:
        if isinstance(label, list):
            label = tuple(label)
        mock_polynomial = {(label,): 1.0}
        return make_BinaryPolynomialModel(mock_polynomial).from_hubo(keys, values)


# classmethods
BinaryPolynomialModel.from_serializable = (
    lambda obj: make_BinaryPolynomialModel_from_JSON(obj).from_serializable(obj)
)
BinaryPolynomialModel.from_hising = (
    lambda *args, **kwargs: make_BinaryPolynomialModel_from_hising(*args, **kwargs)
)
BinaryPolynomialModel.from_hubo = (
    lambda *args, **kwargs: make_BinaryPolynomialModel_from_hubo(*args, **kwargs)
)
