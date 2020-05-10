import unittest

import numpy as np
import cimod
import cxxcimod


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

    def test_bqm_constructor(self):
        # Test BinaryQuadraticModel constructor
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)

        self.assertEqual(bqm.vartype, cimod.SPIN)

        bqm_qubo = cimod.BinaryQuadraticModel.from_qubo(Q=self.Q)
        self.assertEqual(bqm_qubo.vartype, cimod.BINARY)

        bqm = cimod.BinaryQuadraticModel(self.strh, self.strJ)
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)

        self.assertEqual(bqm.vartype, cimod.SPIN)

        bqm_qubo = cimod.BinaryQuadraticModel.from_qubo(Q=self.strQ)
        self.assertEqual(bqm_qubo.vartype, cimod.BINARY)

    def test_interaction_matrix(self):
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        ising_matrix = np.array([
            [1, -1,  0,  0],
            [-1, -2, -3, 0],
            [0, -3, 0, 0.5],
            [0, 0, 0.5, 0]
        ])
        np.testing.assert_array_equal(
            bqm.interaction_matrix(), ising_matrix
        )

        bqm = cimod.BinaryQuadraticModel(self.strh, self.strJ)
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
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        ising_energy_bqm = bqm.energy(self.spins)
        true_ising_e = calculate_ising_energy(self.h, self.J, self.spins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.Q)
        qubo_energy_bqm = bqm.energy(self.binaries)
        true_qubo_e = calculate_qubo_energy(self.Q, self.binaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {0: 1, 1: 1, 2: -1, 3: 1}
        binary = {0: 1, 1: 1, 2: 0, 3: 1}
        qubo_bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.Q)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

        # Test to calculate energy

        # Test Ising energy
        bqm = cimod.BinaryQuadraticModel(self.strh, self.strJ)
        ising_energy_bqm = bqm.energy(self.strspins)
        true_ising_e = calculate_ising_energy(self.strh, self.strJ, self.strspins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.strQ)
        qubo_energy_bqm = bqm.energy(self.strbinaries)
        true_qubo_e = calculate_qubo_energy(self.strQ, self.strbinaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {'a': 1, 'b': 1, 'c': -1, 'd': 1}
        binary = {'a': 1, 'b': 1, 'c': 0, 'd': 1}
        qubo_bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.strQ)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

    def test_change_vartype(self):
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        self.assertEqual(bqm.vartype, cimod.SPIN)
        bqm.change_vartype('BINARY')
        self.assertEqual(bqm.vartype, cimod.BINARY)
        bqm.change_vartype('SPIN', inplace=False)
        self.assertEqual(bqm.vartype, cimod.BINARY)
        


if __name__ == '__main__':
    unittest.main()
