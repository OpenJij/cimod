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

import unittest

import numpy as np
import cimod
import cimod.cxxcimod as cxxcimod
import dimod

def calculate_ising_energy(h, J, spins):
    energy = 0.0
    for (i, j), Jij in J.items():
        energy += Jij*spins[i]*spins[j]
    for i, hi in h.items():
        energy += hi * spins[i]
    return energy


def calculate_qubo_energy(Q, binary):
    energy = 0.0
    for (i, j), Qij in Q.items():
        energy += Qij*binary[i]*binary[j]
    return energy

def calculate_bpm_energy(polynomial, variables):
    energy = 0.0
    for (index, val) in polynomial.items():
        temp = 1
        for site in index:
            temp *= variables[site]
        energy += temp*val
    return energy

class VariableTypeTest(unittest.TestCase):
    def test_variable_type(self):
        spin = cimod.vartype.to_cxxcimod('SPIN')
        self.assertEqual(spin, cxxcimod.Vartype.SPIN)

        binary = cimod.vartype.to_cxxcimod('BINARY')
        self.assertEqual(binary, cxxcimod.Vartype.BINARY)


class ModelTest(unittest.TestCase):

    def setUp(self):
        self.h = {1: -2, 0: 1}
        self.J = {(2, 3): 0.5, (0, 1): -1, (1, 2): -3}
        self.spins = {0: 1, 2: 1, 3: 1, 1: -1}

        self.Q = {(0, 0): 1, (2, 0): -0.2, (1, 3): 3, (1, 2): -1}
        self.binaries = {0: 0, 1: 1, 2: 1, 3: 0}

        self.strh = {'b': -2, 'a': 1}
        self.strJ = {('c', 'd'): 0.5, ('a', 'b'): -1, ('b', 'c'): -3}
        self.strspins = {'a': 1, 'c': 1, 'd': 1, 'b': -1}

        self.strQ = {('a', 'a'): 1, ('c', 'a'): -0.2, ('b', 'd'): 3, ('b', 'c'): -1}
        self.strbinaries = {'a': 0, 'b': 1, 'c': 1, 'd': 0}

        self.tupleh = {(1,2,3): -2, (0,1,2): 1}
        self.tupleJ = {((2,3,4), (3,4,5)): 0.5, ((0,1,2), (1,2,3)): -1, ((1,2,3), (2,3,4)): -3}
        self.tuplespins = {(0,1,2): 1, (2,3,4): 1, (3,4,5): 1, (1,2,3): -1}

        self.tupleQ = {((0,1,2), (0,1,2)): 1, ((2,3,4), (0,1,2)): -0.2, ((1,2,3), (3,4,5)): 3, ((1,2,3), (2,3,4)): -1}
        self.tuplebinaries = {(0,1,2): 0, (1,2,3): 1, (2,3,4): 1, (3,4,5): 0}

    def test_bqm_constructor(self):
        # Test BinaryQuadraticModel constructor
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)
        self.assertEqual(bqm.variables, [0,1,2,3])
        self.assertEqual(bqm.num_variables, 4)

        self.assertEqual(bqm.vartype, cimod.SPIN)

        bqm_qubo = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.Q)
        self.assertEqual(bqm_qubo.vartype, cimod.BINARY)

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.strh, self.strJ, 'BINARY')
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)

        self.assertEqual(bqm.vartype, cimod.BINARY)

        bqm_qubo = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.strQ)
        self.assertEqual(bqm_qubo.vartype, cimod.BINARY)

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.tupleh, self.tupleJ, 2.0, 'SPIN')
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)

        self.assertAlmostEqual(bqm.offset, 2.0)
        self.assertEqual(bqm.vartype, cimod.SPIN)

        bqm_qubo = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)
        self.assertEqual(bqm_qubo.vartype, cimod.BINARY)

        # offset and vartype
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 3, vartype='SPIN')
        self.assertAlmostEqual(bqm.offset, 3)
        self.assertAlmostEqual(bqm.vartype, dimod.SPIN)

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        self.assertAlmostEqual(bqm.offset, 0)
        self.assertAlmostEqual(bqm.vartype, dimod.SPIN)

    def test_interaction_matrix(self):
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        ising_matrix = np.array([
            [1, -1,  0,  0],
            [-1, -2, -3, 0],
            [0, -3, 0, 0.5],
            [0, 0, 0.5, 0]
        ])
        np.testing.assert_array_equal(
            bqm.interaction_matrix(), ising_matrix
        )

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.strh, self.strJ, 'SPIN')
        ising_matrix = np.array([
            [1, -1,  0,  0],
            [-1, -2, -3, 0],
            [0, -3, 0, 0.5],
            [0, 0, 0.5, 0]
        ])
        np.testing.assert_array_equal(
            bqm.interaction_matrix(), ising_matrix
        )

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.tupleh, self.tupleJ, 'SPIN')
        ising_matrix = np.array([
            [1, -1,  0,  0],
            [-1, -2, -3, 0],
            [0, -3, 0, 0.5],
            [0, 0, 0.5, 0]
        ])
        np.testing.assert_array_equal(
            bqm.interaction_matrix(), ising_matrix
        )

    def test_bqm_calc_energy(self):
        # Test to calculate energy

        # Test Ising energy
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        ising_energy_bqm = bqm.energy(self.spins)
        true_ising_e = calculate_ising_energy(self.h, self.J, self.spins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.Q)
        qubo_energy_bqm = bqm.energy(self.binaries)
        true_qubo_e = calculate_qubo_energy(self.Q, self.binaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {0: 1, 1: 1, 2: -1, 3: 1}
        binary = {0: 1, 1: 1, 2: 0, 3: 1}
        qubo_bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.Q)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

        # Test to calculate energy

        # Test Ising energy
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.strh, self.strJ, 'SPIN')
        ising_energy_bqm = bqm.energy(self.strspins)
        true_ising_e = calculate_ising_energy(self.strh, self.strJ, self.strspins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.strQ)
        qubo_energy_bqm = bqm.energy(self.strbinaries)
        true_qubo_e = calculate_qubo_energy(self.strQ, self.strbinaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {'a': 1, 'b': 1, 'c': -1, 'd': 1}
        binary = {'a': 1, 'b': 1, 'c': 0, 'd': 1}
        qubo_bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.strQ)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

        # Test to calculate energy

        # Test Ising energy
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.tupleh, self.tupleJ, 'SPIN')
        ising_energy_bqm = bqm.energy(self.tuplespins)
        true_ising_e = calculate_ising_energy(self.tupleh, self.tupleJ, self.tuplespins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)
        qubo_energy_bqm = bqm.energy(self.tuplebinaries)
        true_qubo_e = calculate_qubo_energy(self.tupleQ, self.tuplebinaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {(0,1,2): 1, (1,2,3): 1, (2,3,4): -1, (3,4,5): 1}
        binary = {(0,1,2): 1, (1,2,3): 1, (2,3,4): 0, (3,4,5): 1}
        qubo_bqm = cimod.model.legacy.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

    def test_change_vartype(self):
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        self.assertEqual(bqm.vartype, cimod.SPIN)
        bqm2 = bqm.change_vartype('BINARY', inplace=True)
        self.assertEqual(bqm.vartype, cimod.BINARY)
        self.assertEqual(bqm.linear, bqm2.linear)
        self.assertEqual(bqm.quadratic, bqm2.quadratic)
        self.assertEqual(bqm.offset, bqm2.offset)
        self.assertEqual(bqm.vartype, bqm2.vartype)
        bqm.change_vartype('SPIN', inplace=False)
        self.assertEqual(bqm.vartype, cimod.BINARY)

    def test_serializable(self):
        bqm = cimod.model.legacy.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        serial = bqm.to_serializable()
        decode_bqm = cimod.model.legacy.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.strh, self.strJ, vartype='SPIN')
        serial = bqm.to_serializable()
        decode_bqm = cimod.model.legacy.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

        bqm = cimod.model.legacy.BinaryQuadraticModel(self.tupleh, self.tupleJ, 'SPIN')
        serial = bqm.to_serializable()
        decode_bqm = cimod.model.legacy.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    # since the latest version of dimod (>= 0.10.0) removes from_serializable function, this test is not needed.
    #def test_serializable_consistent_with_dimod(self):
    #    for (_from,_to) in [(dimod, cimod.model.legacy), (cimod.model.legacy, dimod)]:
    #        bqm = _from.BinaryQuadraticModel(self.h, self.J, vartype='SPIN')
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        self.assertEqual(bqm.linear, decode_bqm.linear)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    #        bqm = _from.BinaryQuadraticModel(self.strh, self.strJ, vartype='SPIN')
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        self.assertEqual(bqm.linear, decode_bqm.linear)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    #        bqm = _from.BinaryQuadraticModel(self.tupleh, self.tupleJ, vartype='SPIN')
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        self.assertEqual(bqm.linear, decode_bqm.linear)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)



if __name__ == '__main__':
    unittest.main()
