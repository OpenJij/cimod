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
import numpy as np

from cimod.utils.decolator import recalc
from cimod.vartype import to_cxxcimod


def make_BinaryQuadraticModel(linear, quadratic):
    """BinaryQuadraticModel factory.
       Generate BinaryQuadraticModel class with the base class specified by the arguments linear and quadratic
    Args:
        linear (dict): linear bias
        quadratic (dict): quadratic bias
    Returns:
        generated BinaryQuadraticModel class
    """
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

        if isinstance(ind, int):
            base = cxxcimod.BinaryQuadraticModel_Dict
        elif isinstance(ind, str):
            base = cxxcimod.BinaryQuadraticModel_str_Dict
        elif isinstance(ind, tuple):
            if len(ind) == 2:
                base = cxxcimod.BinaryQuadraticModel_tuple2_Dict
            elif len(ind) == 3:
                base = cxxcimod.BinaryQuadraticModel_tuple3_Dict
            elif len(ind) == 4:
                base = cxxcimod.BinaryQuadraticModel_tuple4_Dict
            else:
                raise TypeError("invalid length of tuple")
        else:
            raise TypeError("invalid types of linear and quadratic")

    # now define class
    class BinaryQuadraticModel(base):
        """Represents Binary quadratic model.
           Note that the indices are converted to the integers internally.
           The dictionaries between indices and integers are self.ind_to_num (indices -> integers) and self.num_to_ind (integers -> indices).
           Indices are listed in self._indices.
        Attributes:
            vartype (cimod.VariableType): variable type SPIN or BINARY
            linear (dict): represents linear term
            quadratic (dict): represents quadratic term
            adj (dict): represents adjacency
            indices (list): labels of each variables sorted by results variables
            ind_to_num (list): map which specifies where the index is in self._indices
            offset (float): represents constant energy term when convert to SPIN from BINARY
        """

        def __init__(self, linear, quadratic, offset, vartype):
            super().__init__(linear, quadratic, offset, to_cxxcimod(vartype))
            self._init_process()

        def _init_process(self):
            # set recalculate flag True
            # Be sure to enable this flag when variables are changed.
            self._re_calculate = True
            self._re_calculate_indices = True
            # interaction_matrix
            self._interaction_matrix = None
            # indices
            self._indices = None
            # ind to num
            self._ind_to_num = None

        @property
        def linear(self):
            return self.get_linear()

        @property
        def quadratic(self):
            return self.get_quadratic()

        @property
        def num_variables(self):
            return self.get_num_variables()

        @property
        def variables(self):
            return self.get_variables()

        @property
        def vartype(self):
            vartype = super().get_vartype()
            if vartype == cxxcimod.Vartype.SPIN:
                return dimod.SPIN
            else:
                # BINARY
                return dimod.BINARY

        @property
        def offset(self):
            return self.get_offset()

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
                self._ind_to_num = {ind: num for num, ind in enumerate(self._indices)}

                self._re_calculate_indices = False

            return self._indices, self._ind_to_num

        def interaction_matrix(self):
            """make Dense-type interaction matrix
            The Ising model: E = ΣJ_ij σiσj + Σhiσi
                Interaction matrix -> H_ij = J_ij + J_ji, H_ii = hi
            QUBO: E = Σ1/2Q_ij q_iq_j + ΣQ_ii q_i
            Returns:
                numpy.ndarray: interactioin matrix H_{ij} or Q_{ij}, energy_bias (float)
            """

            if self._re_calculate is True:

                # calculate interaction matrix
                indices, ind_to_num = self.update_indices()
                self._interaction_matrix = super().interaction_matrix(indices)
                self._re_calculate = False

            return self._interaction_matrix

        def energy(self, sample, sparse=False, convert_sample=False):
            """Determine the energy of the specified sample of a binary quadratic model.
            Args:
                sample (dict): single sample.
                sparse (bool): if true calculate energy by using adjacency matrix
                convert_sample (bool): if true the sample is automatically converted to self.vartype type.
            Returns:
                energy (float)
            """

            indices, ind_to_num = self.update_indices()

            # convert sample to dict
            if isinstance(sample, list) or isinstance(sample, np.ndarray):
                sample = {indices[i]: elem for i, elem in enumerate(sample)}

            # convert samples to SPIN or BINARY
            if convert_sample:
                for k in sample.keys():
                    if sample[k] == -1 and self.vartype == dimod.BINARY:
                        sample[k] = 0
                    if sample[k] == 0 and self.vartype == dimod.SPIN:
                        sample[k] = -1

            if sparse:
                return super().energy(sample)

            else:
                # convert to array
                if isinstance(sample, dict):
                    state = [0] * len(sample)
                    for k, v in sample.items():
                        state[ind_to_num[k]] = v
                    sample = state

                sample = np.array(sample)

                int_mat = self.interaction_matrix()

                # calculate
                if self.vartype == dimod.BINARY:
                    return (
                        np.dot(sample, np.dot(np.triu(int_mat), sample))
                        + self.get_offset()
                    )
                elif self.vartype == dimod.SPIN:
                    linear_term = np.diag(int_mat)
                    energy = (
                        np.dot(sample, np.dot(int_mat, sample)) - np.sum(linear_term)
                    ) / 2
                    energy += np.dot(linear_term, sample)
                    energy += self.get_offset()
                return energy

        def energies(self, samples_like, **kwargs):
            en_vec = []

            for elem in samples_like:
                en_vec.append(self.energy(elem, **kwargs))

            return en_vec

        @recalc
        def empty(self, *args, **kwargs):
            return super().empty(*args, **kwargs)

        @recalc
        def add_variable(self, *args, **kwargs):
            return super().add_variable(*args, **kwargs)

        @recalc
        def add_variables_from(self, *args, **kwargs):
            return super().add_variables_from(*args, **kwargs)

        @recalc
        def add_interaction(self, *args, **kwargs):
            return super().add_interaction(*args, **kwargs)

        @recalc
        def add_interactions_from(self, *args, **kwargs):
            return super().add_interactions_from(*args, **kwargs)

        @recalc
        def remove_variable(self, *args, **kwargs):
            return super().remove_variable(*args, **kwargs)

        @recalc
        def remove_variables_from(self, *args, **kwargs):
            return super().remove_variables_from(*args, **kwargs)

        @recalc
        def remove_interaction(self, *args, **kwargs):
            return super().remove_interaction(*args, **kwargs)

        @recalc
        def remove_interactions_from(self, *args, **kwargs):
            return super().remove_interactions_from(*args, **kwargs)

        @recalc
        def add_offset(self, *args, **kwargs):
            return super().add_offset(*args, **kwargs)

        @recalc
        def remove_offset(self, *args, **kwargs):
            return super().remove_offset(*args, **kwargs)

        @recalc
        def scale(self, *args, **kwargs):
            return super().scale(*args, **kwargs)

        @recalc
        def normalize(self, *args, **kwargs):
            return super().normalize(*args, **kwargs)

        @recalc
        def fix_variable(self, *args, **kwargs):
            return super().fix_variable(*args, **kwargs)

        @recalc
        def fix_variables(self, *args, **kwargs):
            return super().fix_variables(*args, **kwargs)

        @recalc
        def flip_variable(self, *args, **kwargs):
            return super().flip_variable(*args, **kwargs)

        @recalc
        def update(self, *args, **kwargs):
            return super().update(*args, **kwargs)

        @recalc
        def contract_variables(self, *args, **kwargs):
            return super().contract_variables(*args, **kwargs)

        def change_vartype(self, vartype, inplace=None):
            """
            Create a binary quadratic model with the specified vartype
            Args:
                vartype (cimod.Vartype): SPIN or BINARY
            Returns:
                A new instance of the BinaryQuadraticModel class.
            """
            cxxvartype = to_cxxcimod(vartype)
            if inplace is None:
                super().change_vartype(cxxvartype)
                return

            # in the case inplace is not None
            bqm = super().change_vartype(cxxvartype, inplace)
            self._re_calculate = True
            return BinaryQuadraticModel(
                bqm.get_linear(), bqm.get_quadratic(), bqm.get_offset(), vartype
            )

        @classmethod
        def from_qubo(cls, Q, offset=0.0, **kwargs):
            linear = {}
            quadratic = {}
            for (u, v), bias in Q.items():
                if u == v:
                    linear[u] = bias
                else:
                    quadratic[(u, v)] = bias

            return cls(linear, quadratic, offset, vartype=dimod.BINARY, **kwargs)

        @classmethod
        def from_ising(cls, linear, quadratic, offset=0.0, **kwargs):
            return cls(linear, quadratic, offset, vartype=dimod.SPIN, **kwargs)

        @classmethod
        def from_serializable(cls, obj):

            variable_labels = [
                tuple(elem) if isinstance(elem, list) else elem
                for elem in obj["variable_labels"]
            ]

            # convert to linear biases
            linear = {
                elem: obj["linear_biases"][k] for k, elem in enumerate(variable_labels)
            }

            # convert to quadratic biases
            zipped_obj = zip(
                obj["quadratic_head"], obj["quadratic_tail"], obj["quadratic_biases"]
            )
            quadratic = {
                (variable_labels[elem[0]], variable_labels[elem[1]]): elem[2]
                for elem in zipped_obj
            }

            # set offset
            offset = obj["offset"]

            # set vartype
            vartype = dimod.SPIN if obj["variable_type"] == "SPIN" else dimod.BINARY

            return cls(linear, quadratic, offset, vartype)

    return BinaryQuadraticModel


# for JSON


def make_BinaryQuadraticModel_from_JSON(obj):
    label = obj["variable_labels"][0]
    if isinstance(label, list):
        # convert to tuple
        label = tuple(label)

    mock_linear = {label: 1.0}

    return make_BinaryQuadraticModel(mock_linear, {})


def BinaryQuadraticModel(linear, quadratic, *args, **kwargs):
    Model = make_BinaryQuadraticModel(linear, quadratic)

    if len(args) == 2:
        [offset, vartype] = args
        return Model(linear, quadratic, offset, to_cxxcimod(vartype))
    elif len(args) == 1 and "vartype" in kwargs:
        [offset] = args
        vartype = kwargs["vartype"]
        return Model(linear, quadratic, offset, to_cxxcimod(vartype))
    elif len(args) == 1:
        [vartype] = args
        return Model(linear, quadratic, 0.0, to_cxxcimod(vartype))
    elif len(args) == 0 and "vartype" in kwargs:
        vartype = kwargs["vartype"]
        return Model(linear, quadratic, 0.0, to_cxxcimod(vartype))
    else:
        raise TypeError("invalid args for BinaryQuadraticModel")


# classmethods
BinaryQuadraticModel.from_qubo = (
    lambda Q, offset=0.0, **kwargs: make_BinaryQuadraticModel({}, Q).from_qubo(
        Q, offset, **kwargs
    )
)

BinaryQuadraticModel.from_ising = (
    lambda linear, quadratic, offset=0.0, **kwargs: make_BinaryQuadraticModel(
        linear, quadratic
    ).from_ising(linear, quadratic, offset, **kwargs)
)

BinaryQuadraticModel.from_serializable = (
    lambda obj: make_BinaryQuadraticModel_from_JSON(obj).from_serializable(obj)
)
