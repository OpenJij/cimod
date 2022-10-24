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

import pprint

import cimod
import cimod.cxxcimod as cxxcimod
import dimod
import numpy as np

from cimod.vartype import to_cxxcimod
from scipy.sparse import dok_matrix, csr_matrix

from typing import Tuple, Union
from collections import defaultdict


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
        # assume that index type is int
        ind = int(0)
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


def extract_offset_and_vartype(*args, **kwargs):
    if kwargs == {}:
        if len(args) == 0:
            raise TypeError(
                f"Offset or vartype is configured incorrectly. Vartype must be set."
            )
        elif len(args) == 1:
            offset = 0.0
            [vartype] = args
        elif len(args) == 2:
            [offset, vartype] = args
        else:
            raise TypeError(
                f"Offset or vartype is configured incorrectly. Vartype must be set."
            )
    else:
        if "offset" in kwargs and "vartype" in kwargs:
            offset = kwargs["offset"]
            vartype = kwargs["vartype"]
        elif "offset" in kwargs:
            if len(args) != 1:
                raise TypeError(
                    f"Offset or vartype is configured incorrectly. Vartype must be set."
                )
            offset = kwargs["offset"]
            [vartype] = args
        elif "vartype" in kwargs:
            if len(args) >= 2:
                raise TypeError(
                    f"Offset or vartype is configured incorrectly. Vartype must be set."
                )
            elif len(args) == 0:
                offset = 0.0
            elif len(args) == 1:
                [offset] = args
            vartype = kwargs["vartype"]
        else:
            raise TypeError(
                f"Offset or vartype is configured incorrectly. Vartype must be set."
            )

    return offset, vartype


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
            # replace vartype with cxxcimod type
            vartypes = [
                dimod.SPIN,
                dimod.BINARY,
                cimod.SPIN,
                cimod.BINARY,
                "SPIN",
                "BINARY",
            ]
            args = [
                to_cxxcimod(elem)
                if not isinstance(elem, (np.ndarray, csr_matrix)) and vartypes.count(elem) != 0
                else elem
                for elem in args
            ]
            kwargs = {
                k: to_cxxcimod(v)
                if not isinstance(v, (np.ndarray, csr_matrix)) and vartypes.count(v) != 0
                else v
                for k, v in kwargs.items()
            }

            self.model_type = "cimod.BinaryQuadraticModel"
            # if linear and quadratic are given. generate matrix and initialize as matrix
            if (
                len(args) >= 2
                and isinstance(args[0], dict)
                and isinstance(args[1], dict)
            ):
                linear = args[0]
                quadratic = args[1]
                offset, vartype = extract_offset_and_vartype(*args[2:], **kwargs)

                mat, idx_to_label = self._generate_mat(linear, quadratic, False, sparse)

                if sparse is False:
                    super().__init__(mat, idx_to_label, offset, vartype, fix_format=False)
                else:
                    super().__init__(mat, idx_to_label, offset, vartype)

            else:
                super().__init__(*args, **kwargs)

        @staticmethod
        def _generate_mat(linear: dict, quadratic: dict, include_quaddiag: bool, sparse: bool) -> Tuple[Union[np.ndarray, csr_matrix], dict]:
            labels = set()

            for i in linear.keys():
                labels.add(i)

            for i, j in quadratic.keys():
                if (i == j) and (include_quaddiag == False):
                    pp_quadratic = pprint.pformat(quadratic)
                    pp_include_quaddiag = pprint.pformat(include_quaddiag)
                    pp_linear = pprint.pformat(linear)
                    print("quadratic : {0}".format(pp_quadratic))
                    print("include_quaddiag : {0}".format(pp_include_quaddiag))
                    print("linear : {0}".format(pp_linear))
                    raise RuntimeError("No self-loop allowed")

                labels.add(i)
                labels.add(j)

            idx_to_label = sorted(labels)
            label_to_idx = {elem: k for k, elem in enumerate(idx_to_label)}

            mat_size = len(idx_to_label) + 1


            # generate matrix (dense or sparse)
            if sparse is False:
                mat = np.zeros(shape=(mat_size, mat_size))
            else:
                # first defines dict and use `_update` function to make `dok_matrix`
                # then convert to `csr_matrix`
                mat = defaultdict(float)

            mat[mat_size - 1, mat_size - 1] = 1

            for i, val in linear.items():
                mat[label_to_idx[i], mat_size - 1] += val

            for (i, j), val in quadratic.items():
                idx_i = label_to_idx[i]
                idx_j = label_to_idx[j]
                if idx_i != idx_j:
                    mat[min(idx_i, idx_j), max(idx_i, idx_j)] += val
                else:
                    mat[idx_i, mat_size - 1] += val

            if sparse is False:
                return mat, idx_to_label
            else:
                dok_mat = dok_matrix((mat_size, mat_size))
                # using `_update` function skips index checking
                dok_mat._update(mat)
                csr_mat = dok_mat.tocsr()
                csr_mat.sort_indices()
                return csr_mat, idx_to_label


        @property
        def vartype(self):
            vartype = super().get_vartype()
            if vartype == cxxcimod.Vartype.SPIN:
                return dimod.SPIN
            else:
                # BINARY
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

        @property
        def sparse(self):
            return sparse

        def empty(self, vartype):
            return self.__class__(super().empty(to_cxxcimod(vartype)))

        def change_vartype(self, vartype, inplace=None):
            if inplace is None:
                super().change_vartype(to_cxxcimod(vartype))
            else:
                return self.__class__(
                    super().change_vartype(to_cxxcimod(vartype), inplace)
                )

        def __str__(self):
            return f"BinaryQuadraticModel({self.linear}, {self.quadratic}, {self.offset}, {self.vartype}, sparse={self.sparse})"

        def __repr__(self):
            return f"BinaryQuadraticModel({self.linear}, {self.quadratic}, {self.offset}, {self.vartype}, sparse={self.sparse})"

        def energy(self, sample):
            if isinstance(sample, list):
                sample = {self.variables[k]: elem for k, elem in enumerate(sample)}

            return super().energy(sample)

        def energies(self, samples_like):
            if isinstance(samples_like[0], list):
                samples_like = [
                    {self.variables[k]: elem for k, elem in enumerate(inner_array)}
                    for inner_array in samples_like
                ]

            return super().energies(samples_like)

        @classmethod
        def from_numpy_matrix(
            cls,
            mat,
            variables: list,
            offset=0.0,
            vartype="BINARY",
            fix_format=True,
            **kwargs,
        ):
            shape = np.shape(mat)
            if not (len(shape) == 2 and shape[0] == shape[1]):
                raise TypeError("numpy matrix has to be a square matrix")

            return cls(mat, variables, offset, vartype, fix_format, **kwargs)

        @classmethod
        def from_qubo(cls, Q, offset=0.0, **kwargs):
            # cxxbqm = Base.from_qubo(Q, offset)
            # return cls(cxxbqm, **kwargs)
            # separate linear and quadratic

            mat, idx_to_label = cls._generate_mat({}, Q, True, sparse)

            return cls(mat, idx_to_label, offset, "BINARY", **kwargs)

        @classmethod
        def from_ising(cls, linear, quadratic, offset=0.0, **kwargs):
            # cxxbqm = Base.from_ising(linear, quadratic, offset)
            # return cls(cxxbqm, **kwargs)
            mat, idx_to_label = cls._generate_mat(linear, quadratic, False, sparse)

            return cls(mat, idx_to_label, offset, "SPIN", **kwargs)

        @classmethod
        def from_serializable(cls, obj, **kwargs):
            cxxbqm = Base.from_serializable(obj)
            return cls(cxxbqm, **kwargs)

    return BinaryQuadraticModel


# for JSON


def make_BinaryQuadraticModel_from_JSON(obj):
    label = obj["variable_labels"][0]
    if isinstance(label, list):
        # convert to tuple
        label = tuple(label)

    mock_linear = {label: 1.0}

    if obj["version"]["bqm_schema"] == "3.0.0-dense":
        sparse = False
    elif obj["version"]["bqm_schema"] == "3.0.0":
        sparse = True
    else:
        raise TypeError("Invalid bqm_schema")

    return make_BinaryQuadraticModel(mock_linear, {}, sparse)


def BinaryQuadraticModel(linear, quadratic, *args, **kwargs):
    sparse_option = kwargs.pop("sparse", False)
    Model = make_BinaryQuadraticModel(linear, quadratic, sparse_option)

    # offset and vartype
    offset, vartype = extract_offset_and_vartype(*args, **kwargs)

    return Model(linear, quadratic, offset, vartype)


def bqm_from_numpy_matrix(
    mat, variables: list = None, offset=0.0, vartype="BINARY", **kwargs
):
    if variables is None:
        # generate array
        num_variables = mat.shape[0]
        variables = list(range(num_variables))

    sparse_option = kwargs.pop("sparse", False)

    return make_BinaryQuadraticModel(
        {variables[0]: 1.0}, {}, sparse_option
    ).from_numpy_matrix(mat, variables, offset, vartype, True, **kwargs)


BinaryQuadraticModel.from_numpy_matrix = bqm_from_numpy_matrix

def bqm_from_qubo(Q, offset=0.0, **kwargs):
    sparse_option = kwargs.pop("sparse", False)

    return make_BinaryQuadraticModel(
        {}, Q, sparse_option
    ).from_qubo(Q, offset, **kwargs)

BinaryQuadraticModel.from_qubo = bqm_from_qubo


def bqm_from_ising(linear, quadratic, offset=0.0, **kwargs):
    sparse_option = kwargs.pop("sparse", False)

    return make_BinaryQuadraticModel(
        linear, quadratic, sparse_option
    ).from_ising(linear, quadratic, offset, **kwargs)

BinaryQuadraticModel.from_ising = bqm_from_ising


BinaryQuadraticModel.from_serializable = (
    lambda obj, **kwargs: make_BinaryQuadraticModel_from_JSON(obj).from_serializable(
        obj, **kwargs
    )
)
