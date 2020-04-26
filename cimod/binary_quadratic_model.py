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
from cimod.vartype import to_cxxcimod
import dimod

class BinaryQuadraticModel(cxxcimod.BinaryQuadraticModel):
    """Represents Binary quadratic model. 
       Note that the indices are converted to the integers internally. 
       The dictionaries between indices and integers are self.ind_to_num (indices -> integers) and self.num_to_ind (integers -> indices).
       Indices are listed in self.indices.
    Attributes:
        var_type (cimod.VariableType): variable type SPIN or BINARY
        linear (dict): represents linear term
        quad (dict): represents quadratic term
        indices (list): labels of each variables sorted by results variables
        offset (float): represents constant energy term when convert to SPIN from BINARY
        size (int): number of variables
    """
    def __init__(self, linear, quadratic, offset=0.0,
                 var_type=dimod.SPIN, **kwargs):
        # set recalculate flag True
        # Be sure to enable this flag when variables are changed.
        self._re_calculate = True
        # interaction_matrix
        self._interaction_matrix = None

        self.indices, self.num_to_ind, self.ind_to_num = self._generate_indices_dict(linear, quadratic)

        # convert indices to integers and call super constructor
        linear      = self._conv_linear(linear, self.ind_to_num)
        quadratic   = self._conv_quadratic(quadratic, self.ind_to_num)
        super().__init__(linear, quadratic, offset, to_cxxcimod(var_type))

    def get_linear(self, original_ind=True):
        """
        get linear
        Args:
            original_ind (bool): if true returns linear with index converted to original one
        Returns:
            linear (dict)
        """
        linear = super().get_linear()

        if original_ind:
            linear = self._conv_linear(linear, self.num_to_ind)

        return linear

    def get_quadratic(self, original_ind=True):
        """
        get quadratic
        Args:
            original_ind (bool): if true returns linear with index converted to original one
        Returns:
            quadratic (dict)
        """
        quadratic = super().get_quadratic()

        if original_ind:
            quadratic = self._conv_quadratic(quadratic, self.num_to_ind)

        return quadratic

    def get_adjacency(self, original_ind=True):
        """
        get adjacency
        Args:
            original_ind (bool): if true returns linear with index converted to original one
        Returns:
            adjacency (dict)
        """
        adjacency = super().get_adjacency()

        if original_ind:
            adjacency = self._conv_adjacency(adjacency, self.num_to_ind)

        return adjacency

    @property
    def linear(self):
        return self.get_linear()

    @property
    def quadratic(self):
        return self.get_quadratic()

    @property
    def adj(self):
        return self.get_adjacency()

    @property
    def vartype(self):
        vartype = self.get_vartype()
        if vartype == cxxcimod.Vartype.SPIN:
            return dimod.SPIN
        else:
            #BINARY
            return dimod.BINARY

    @property
    def offset(self):
        return self.get_offset()

    def to_qubo(self, original_ind=True):
        """
        Convert a binary quadratic model to QUBO format.
        Args:
            original_ind (bool): if true returns linear with index converted to original one
        Returns:
            Q (dict), offset
        """
        Q, offset = super().to_qubo()

        if original_ind:
            Q = self._conv_quadratic(Q, self.num_to_ind)

        return Q, offset

    def to_ising(self, original_ind=True):
        """
        Convert a binary quadratic model to Ising format.
        Args:
            original_ind (bool): if true returns linear with index converted to original one
        Returns:
            h (dict), J (dict), offset
        """

        h, J, offset = super().to_ising()

        if original_ind:
            h = self._conv_linear(h, self.num_to_ind)
            J = self._conv_quadratic(J, self.num_to_ind)

        return h, J, offset

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
            ind_list = [self.ind_to_num[elem] for elem in self.indices]
            self._interaction_matrix = super().interaction_matrix(ind_list)
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

        if isinstance(sample, dict):
            # convert int to num
            sample = self._conv_linear(sample, self.ind_to_num)

        # convert samples to SPIN or BINARY
        if convert_sample:
            for i in range(len(sample)):
                if sample[i] == -1 and self.vartype == openjij.BINARY:
                    sample[i] = 0
                if sample[i] == 0  and self.vartype == openjij.SPIN:
                    sample[i] = -1

        if sparse:
           # convert sample to dict
           if isinstance(sample, list) or isinstance(sample, np.ndarray):
               sample = {i:elem for i,elem in enumerate(sample)}
           return super().energy(sample)

        else:
            if isinstance(sample, dict):
                state = [0] * len(sample)
                for k,v in sample.items():
                    state[k] = v
                sample = state

            sample = np.array(sample)

            int_mat = self.interaction_matrix()

            # calculate 
            if self.get_vartype() == openjij.BINARY:
                return np.dot(sample, np.dot(np.triu(int_mat), sample)) + self.get_offset()
            elif self.get_vartype() == openjij.SPIN:
                linear_term = np.diag(int_mat)
                energy = (np.dot(sample, np.dot(int_mat, sample)) -
                      np.sum(linear_term))/2
                energy += np.dot(linear_term, sample)
                energy += self.get_offset()
            return energy 

    def energies(self, samples_like, **kwargs):
        en_vec = []

        for elem in samples_like:
            en_vec.append(self.energy(elem, **kwargs))

        return en_vec

    @staticmethod
    def _generate_indices_dict(linear=None, quadratic=None):
        """
        Generate indices dictionaries.
        Args:
            linear (dict), quadratic (dict)
        Returns:
            tuple of dictionaries (indices, num_to_ind, ind_to_num)
        """
        if linear is not None:
            index_set = set(linear.keys())
        else:
            index_set = set()


        if quadratic is not None:
            for v1, v2 in quadratic.keys():
                index_set.add(v1)
                index_set.add(v2)

        indices = list(index_set)

        # generate conversion map index <-> integer
        num_to_ind = {k:val for k,val in enumerate(indices)}
        ind_to_num = {val:k for k,val in enumerate(indices)}

        return indices, num_to_ind, ind_to_num

    def _conv_linear(self, dic, conv_dict):
        """
        Convert indices of dictionary (linear)
        Args:
            dic (dict): dictionary
            conv_dict (dict): convert dict (ind_to_num or num_to_ind)
        Returns:
            dictionaries with indices converted
        """
        return {conv_dict[k]:v for k,v in dic.items()}

    def _conv_quadratic(self, dic, conv_dict):
        """
        Convert indices of dictionary (quadratic)
        Args:
            dic (dict): dictionary
            conv_dict (dict): convert dict (ind_to_num or num_to_ind)
        Returns:
            dictionaries with indices converted
        """
        return {(conv_dict[k1], conv_dict[k2]):v for (k1,k2),v in dic.items()}

    def _conv_adjacency(self, dic, conv_dict):
        """
        Convert indices of dictionary (adjacency)
        Args:
            dic (dict): dictionary
            conv_dict (dict): convert dict (ind_to_num or num_to_ind)
        Returns:
            dictionaries with indices converted
        """
        return {conv_dict[index]:{conv_dict[k]:v for k,v in adj_dic.items()} for index,adj_dic in dic.items()}
