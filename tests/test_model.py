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
from scipy.sparse import csr_matrix
import cimod
import cimod.cxxcimod as cxxcimod
import dimod
import random

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

    def test_empty(self):
        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
            bqm2 = bqm.empty('SPIN')
            self.assertEqual(type(bqm2.interaction_matrix()), mat_type)
            self.assertEqual(bqm2.variables, [])
            self.assertEqual(bqm2.num_variables, 0)
            self.assertEqual(bqm2.vartype, cimod.SPIN)
            bqm2 = bqm.empty('BINARY')
            self.assertEqual(type(bqm2.interaction_matrix()), mat_type)
            self.assertEqual(bqm2.variables, [])
            self.assertEqual(bqm2.num_variables, 0)
            self.assertEqual(bqm2.vartype, cimod.BINARY)

    def test_bqm_constructor(self):
        # Test BinaryQuadraticModel constructor

        bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN')
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)
        self.assertEqual(bqm.variables, [0,1,2,3])
        self.assertEqual(bqm.num_variables, 4)

        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)
            self.assertEqual(bqm.variables, [0,1,2,3])
            self.assertEqual(bqm.num_variables, 4)

            self.assertEqual(bqm.vartype, cimod.SPIN)

            bqm_qubo = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.Q, sparse=sparse)
            self.assertEqual(bqm_qubo.vartype, cimod.BINARY)
            self.assertEqual(type(bqm_qubo.interaction_matrix()), mat_type)

            bqm = cimod.model.BinaryQuadraticModel(self.strh, self.strJ, 'BINARY', sparse=sparse)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)

            self.assertEqual(bqm.vartype, cimod.BINARY)

            bqm_qubo = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.strQ, sparse=sparse)
            self.assertEqual(bqm_qubo.vartype, cimod.BINARY)
            self.assertEqual(type(bqm_qubo.interaction_matrix()), mat_type)

            bqm = cimod.model.BinaryQuadraticModel(self.tupleh, self.tupleJ, 2.0, 'SPIN', sparse=sparse)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)

            self.assertAlmostEqual(bqm.offset, 2.0)
            self.assertEqual(bqm.vartype, cimod.SPIN)

            bqm_qubo = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.tupleQ, sparse=sparse)
            self.assertEqual(bqm_qubo.vartype, cimod.BINARY)
            self.assertEqual(type(bqm_qubo.interaction_matrix()), mat_type)

            # offset and vartype
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 3, vartype='SPIN', sparse=sparse)
            self.assertAlmostEqual(bqm.offset, 3)
            self.assertAlmostEqual(bqm.vartype, dimod.SPIN)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)

            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
            self.assertAlmostEqual(bqm.offset, 0)
            self.assertAlmostEqual(bqm.vartype, dimod.SPIN)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)

            # empty object
            bqm = cimod.model.BinaryQuadraticModel({}, {}, 'SPIN', sparse=sparse)
            self.assertAlmostEqual(bqm.offset, 0)
            self.assertAlmostEqual(bqm.vartype, dimod.SPIN)
            self.assertEqual(type(bqm.interaction_matrix()), mat_type)

    def test_bqm_calc_energy(self):
        # Test to calculate energy

        # Test Ising energy
        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
            ising_energy_bqm = bqm.energy(self.spins)
            true_ising_e = calculate_ising_energy(self.h, self.J, self.spins)
            self.assertEqual(ising_energy_bqm, true_ising_e)
            ising_energy_bqm = bqm.energies([self.spins, self.spins])
            self.assertEqual(ising_energy_bqm, [true_ising_e, true_ising_e])

            # Test QUBO energy
            bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.Q, sparse=sparse)
            qubo_energy_bqm = bqm.energy(self.binaries)
            true_qubo_e = calculate_qubo_energy(self.Q, self.binaries)
            self.assertEqual(qubo_energy_bqm, true_qubo_e)

            # QUBO == Ising
            spins = {0: 1, 1: 1, 2: -1, 3: 1}
            binary = {0: 1, 1: 1, 2: 0, 3: 1}
            qubo_bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.Q, sparse=sparse)

            qubo_energy = qubo_bqm.energy(binary)
            qubo_bqm.change_vartype('SPIN')

            self.assertEqual(qubo_energy, qubo_bqm.energy(spins))

            # QUBO == Ising
            spins = [1,1,-1,1]
            binary = [1,1,0,1]
            qubo_bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.Q, sparse=sparse)

            qubo_energy = qubo_bqm.energy(binary)
            qubo_bqm.change_vartype('SPIN')

            self.assertEqual(qubo_energy, qubo_bqm.energy(spins))

            spins = [[1,1,-1,1], [-1,-1,1,1]]
            binary = [[1,1,0,1], [0,0,1,1]]

            # Test to calculate energy
            qubo_bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.Q, sparse=sparse)

            qubo_energy = qubo_bqm.energies(binary)
            qubo_bqm.change_vartype('SPIN')

            self.assertEqual(qubo_energy, qubo_bqm.energies(spins))

            # Test Ising energy
            bqm = cimod.model.BinaryQuadraticModel(self.strh, self.strJ, 'SPIN', sparse=sparse)
            ising_energy_bqm = bqm.energy(self.strspins)
            true_ising_e = calculate_ising_energy(self.strh, self.strJ, self.strspins)
            self.assertEqual(ising_energy_bqm, true_ising_e)

            # Test QUBO energy
            bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.strQ, sparse=sparse)
            qubo_energy_bqm = bqm.energy(self.strbinaries)
            true_qubo_e = calculate_qubo_energy(self.strQ, self.strbinaries)
            self.assertEqual(qubo_energy_bqm, true_qubo_e)

            # QUBO == Ising
            spins = {'a': 1, 'b': 1, 'c': -1, 'd': 1}
            binary = {'a': 1, 'b': 1, 'c': 0, 'd': 1}
            qubo_bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.strQ, sparse=sparse)

            qubo_energy = qubo_bqm.energy(binary)
            qubo_bqm.change_vartype('SPIN')

            self.assertEqual(qubo_energy, qubo_bqm.energy(spins))

            # Test to calculate energy

            # Test Ising energy
            bqm = cimod.model.BinaryQuadraticModel(self.tupleh, self.tupleJ, 'SPIN', sparse=sparse)
            ising_energy_bqm = bqm.energy(self.tuplespins)
            true_ising_e = calculate_ising_energy(self.tupleh, self.tupleJ, self.tuplespins)
            self.assertEqual(ising_energy_bqm, true_ising_e)

            # Test QUBO energy
            bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.tupleQ, sparse=sparse)
            qubo_energy_bqm = bqm.energy(self.tuplebinaries)
            true_qubo_e = calculate_qubo_energy(self.tupleQ, self.tuplebinaries)
            self.assertEqual(qubo_energy_bqm, true_qubo_e)

            # QUBO == Ising
            spins = {(0,1,2): 1, (1,2,3): 1, (2,3,4): -1, (3,4,5): 1}
            binary = {(0,1,2): 1, (1,2,3): 1, (2,3,4): 0, (3,4,5): 1}
            qubo_bqm = cimod.model.BinaryQuadraticModel.from_qubo(Q=self.tupleQ, sparse=sparse)

            qubo_energy = qubo_bqm.energy(binary)
            qubo_bqm.change_vartype('SPIN')

            self.assertEqual(qubo_energy, qubo_bqm.energy(spins))

    def test_change_vartype(self):

        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
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
        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(self.h, self.J, 'SPIN', sparse=sparse)
            serial = bqm.to_serializable()
            decode_bqm = cimod.model.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)

            bqm = cimod.model.BinaryQuadraticModel(self.strh, self.strJ, vartype='SPIN', sparse=sparse)
            serial = bqm.to_serializable()
            decode_bqm = cimod.model.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)

            bqm = cimod.model.BinaryQuadraticModel(self.tupleh, self.tupleJ, 'SPIN', sparse=sparse)
            serial = bqm.to_serializable()
            decode_bqm = cimod.model.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)

    # since the latest version of dimod (>= 0.10.0) removes from_serializable function, this test is not needed.

    #def test_serializable_consistent_with_dimod(self):
    #    for (_from,_to) in [(dimod, cimod), (cimod, dimod)]:
    #        bqm = _from.BinaryQuadraticModel(self.h, self.J, vartype='SPIN', sparse=True)
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        l1 = {k: v for (k,v) in bqm.linear.items() if v != 0}
    #        l2 = {k: v for (k,v) in decode_bqm.linear.items() if v != 0}
    #        self.assertEqual(l1, l2)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    #        bqm = _from.BinaryQuadraticModel(self.strh, self.strJ, vartype='SPIN', sparse=True)
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        l1 = {k: v for (k,v) in bqm.linear.items() if v != 0}
    #        l2 = {k: v for (k,v) in decode_bqm.linear.items() if v != 0}
    #        self.assertEqual(l1, l2)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    #        bqm = _from.BinaryQuadraticModel(self.tupleh, self.tupleJ, vartype='SPIN', sparse=True)
    #        serial = bqm.to_serializable()
    #        decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
    #        l1 = {k: v for (k,v) in bqm.linear.items() if v != 0}
    #        l2 = {k: v for (k,v) in decode_bqm.linear.items() if v != 0}
    #        self.assertEqual(l1, l2)
    #        # order of indices in quadratic is not considered.
    #        self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
    #        self.assertEqual(bqm.offset, decode_bqm.offset)
    #        self.assertEqual(bqm.vartype, decode_bqm.vartype)


#BinaryPolynomialModel
class PolynomialModelTest(unittest.TestCase):
    def setUp(self):
        self.poly     = {(1,):1.0, (3,):3.0, (1,2):12.0, (1,3):13.0, (2,3,4):234.0, (3,5):35.0}
        self.spins    = {1:+1, 2:-1, 3:+1, 4:-1, 5:+1} 
        self.binaries = {1: 1, 2: 0, 3: 1, 4: 0, 5: 1}

        self.poly_str     = {("a",):1.0, ("c",):3.0, ("a","b"):12.0, ("a","c"):13.0, ("b","c","d"):234.0, ("c","e"):35.0}
        self.spins_str    = {"a":+1, "b":-1, "c":+1, "d":-1, "e":+1} 
        self.binaries_str = {"a": 1, "b": 0, "c": 1, "d": 0, "e": 1}

        self.poly_tuple2     = {((1,1),):1.0, ((3,3),):3.0, ((1,1),(2,2)):12.0, ((1,1),(3,3)):13.0, ((2,2),(3,3),(4,4)):234.0, ((3,3),(5,5)):35.0}
        self.spins_tuple2    = {(1,1):+1, (2,2):-1, (3,3):+1, (4,4):-1, (5,5):+1} 
        self.binaries_tuple2 = {(1,1): 1, (2,2): 0, (3,3): 1, (4,4): 0, (5,5): 1}

        self.poly_tuple3     = {((1,1,1),):1.0, ((3,3,3),):3.0, ((1,1,1),(2,2,2)):12.0, ((1,1,1),(3,3,3)):13.0, ((2,2,2),(3,3,3),(4,4,4)):234.0, ((3,3,3),(5,5,5)):35.0}
        self.spins_tuple3    = {(1,1,1):+1, (2,2,2):-1, (3,3,3):+1, (4,4,4):-1, (5,5,5):+1} 
        self.binaries_tuple3 = {(1,1,1): 1, (2,2,2): 0, (3,3,3): 1, (4,4,4): 0, (5,5,5): 1}

        self.poly_tuple4     = {((1,1,1,1),):1.0, ((3,3,3,3),):3.0, ((1,1,1,1),(2,2,2,2)):12.0, ((1,1,1,1),(3,3,3,3)):13.0, ((2,2,2,2),(3,3,3,3),(4,4,4,4)):234.0, ((3,3,3,3),(5,5,5,5)):35.0}
        self.spins_tuple4    = {(1,1,1,1):+1, (2,2,2,2):-1, (3,3,3,3):+1, (4,4,4,4):-1, (5,5,5,5):+1} 
        self.binaries_tuple4 = {(1,1,1,1): 1, (2,2,2,2): 0, (3,3,3,3): 1, (4,4,4,4): 0, (5,5,5,5): 1}

    def state_test_bpm(self, bpm, poly: dict, vartype):
        self.assertEqual(bpm.vartype, vartype) #Check Vartype
        self.assertEqual(bpm.num_interactions, len(poly)) #Check the number of the interactions
        self.assertEqual(bpm.num_variables, len(set(j for i in poly.keys() for j in i))) #Check the number of the variables
        self.assertEqual(bpm.degree, max([len(i) for i in poly.keys()])) #Check the max degree of the interactions
        self.assertEqual(bpm.get_variables(), sorted(list(set(j for i in poly.keys() for j in i)))) #Check the variables
        self.assertAlmostEqual(bpm.get_offset(), poly[()] if tuple() in poly else 0.0) #Check the offset
        for k, v in bpm.get_polynomial().items():#Check the interactions
            self.assertAlmostEqual(v, poly[k])

        num = 0
        for i in sorted(list(set(j for i in poly.keys() for j in i))):
            self.assertEqual(bpm.get_variables_to_integers(i), num)
            self.assertEqual(bpm.has_variable(i), True)
            num += 1
        
         #Check the specific interactions 
        for index in poly.keys():
            self.assertAlmostEqual(bpm.get_polynomial(index)                           , poly[index])
            self.assertAlmostEqual(bpm.get_polynomial(random.sample(index, len(index))), poly[index])
            self.assertAlmostEqual(bpm.get_polynomial(list(index))                     , poly[index])
            self.assertAlmostEqual(bpm.get_polynomial(key = index)                     , poly[index])
            if tuple(index) != ():
                self.assertAlmostEqual(bpm.get_polynomial(*index), poly[index])
            else:
                self.assertAlmostEqual(bpm.get_polynomial(index), poly[index])

    def state_test_bpm_empty(self, bpm, vartype):
        self.assertEqual(bpm.vartype, vartype)
        self.assertEqual(bpm.num_interactions, 0)
        self.assertEqual(bpm.num_variables, 0)
        self.assertEqual(bpm.degree, 0)
        self.assertEqual(bpm.get_polynomial(), {})
        self.assertEqual(bpm.get_variables_to_integers(), {})
        self.assertEqual(bpm.get_variables(), [])
        self.assertAlmostEqual(bpm.get_offset(), 0.0)

    # Test BinaryPolynomialModel constructor
    def test_construction_bpm(self):
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly       , cimod.SPIN), self.poly       , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_str   , cimod.SPIN), self.poly_str   , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN), self.poly_tuple2, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN), self.poly_tuple3, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN), self.poly_tuple4, cimod.SPIN)

        self.state_test_bpm(cimod.BinaryPolynomialModel(list(self.poly.keys())       , list(self.poly.values())       , cimod.SPIN), self.poly       , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(list(self.poly_str.keys())   , list(self.poly_str.values())   , cimod.SPIN), self.poly_str   , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(list(self.poly_tuple2.keys()), list(self.poly_tuple2.values()), cimod.SPIN), self.poly_tuple2, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(list(self.poly_tuple3.keys()), list(self.poly_tuple3.values()), cimod.SPIN), self.poly_tuple3, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(list(self.poly_tuple4.keys()), list(self.poly_tuple4.values()), cimod.SPIN), self.poly_tuple4, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({}, "SPIN")
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel([], [], "SPIN")
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_add_interaction_bpm_basic(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.add_interaction((-12345, -2, 897654321), 0.1234567)
        self.poly[tuple(sorted([-12345, -2, 897654321]))] = 0.1234567
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.add_interaction(("åß∂ƒ©", "あいうえお", "ABCD"), -123)
        self.poly_str[tuple(sorted(["åß∂ƒ©", "あいうえお", "ABCD"]))] = -123
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.add_interaction(((-11, 11), (22, 22)), -123)
        self.poly_tuple2[tuple(sorted([(-11, 11), (22, 22)]))] = -123
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.add_interaction(((-11, 11, -321), (22, 22, -321)), -123)
        self.poly_tuple3[tuple(sorted([(-11, 11, -321), (22, 22, -321)]))] = -123
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.add_interaction(((-11, 11, -321, 0), (22, 22, -321, 0)), -123)
        self.poly_tuple4[tuple(sorted([(-11, 11, -321, 0), (22, 22, -321, 0)]))] = -123
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_add_interaction_bpm_duplicate_value_1(self):
        bpm = cimod.BinaryPolynomialModel({k: v/2 for k, v in self.poly.items()}, cimod.SPIN)
        for k, v in self.poly.items():
            bpm.add_interaction(k, v/2)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: v/2 for k, v in self.poly_str.items()}, cimod.SPIN)
        for k, v in self.poly_str.items():
            bpm.add_interaction(k, v/2)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: v/2 for k, v in self.poly_tuple2.items()}, cimod.SPIN)
        for k, v in self.poly_tuple2.items():
            bpm.add_interaction(k, v/2)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: v/2 for k, v in self.poly_tuple3.items()}, cimod.SPIN)
        for k, v in self.poly_tuple3.items():
            bpm.add_interaction(k, v/2)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: v/2 for k, v in self.poly_tuple4.items()}, cimod.SPIN)
        for k, v in self.poly_tuple4.items():
            bpm.add_interaction(k, v/2)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_add_interaction_bpm_duplicate_value_2(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        for k, v in self.poly.items():
            bpm.add_interaction(k, -v)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        for k, v in self.poly_str.items():
            bpm.add_interaction(k, -v)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        for k, v in self.poly_tuple2.items():
            bpm.add_interaction(k, -v)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        for k, v in self.poly_tuple3.items():
            bpm.add_interaction(k, -v)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        for k, v in self.poly_tuple4.items():
            bpm.add_interaction(k, -v)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_add_interaction_bpm_self_loop(self):
        def check(poly, vartype):
            bpm = cimod.BinaryPolynomialModel(poly, vartype)
            for k, v in poly.items():
                key = []
                for i in k:
                    key.append(i)
                for i in k:
                    key.append(i)
                for i in k:
                    key.append(i)
                bpm.add_interaction(key, -v)
            self.state_test_bpm_empty(bpm, vartype)
        
        check(self.poly, cimod.SPIN)
        check(self.poly_str, cimod.SPIN)
        check(self.poly_tuple2, cimod.SPIN)
        check(self.poly_tuple3, cimod.SPIN)
        check(self.poly_tuple4, cimod.SPIN)

        check(self.poly, cimod.BINARY)
        check(self.poly_str, cimod.BINARY)
        check(self.poly_tuple2, cimod.BINARY)
        check(self.poly_tuple3, cimod.BINARY)
        check(self.poly_tuple4, cimod.BINARY)

    def test_get_polynomial_self_loop(self):
        bpm = cimod.BinaryPolynomialModel({(): 0.1, (1,): 0.2, (1,2): 0.3}, "SPIN")
        self.assertAlmostEqual(bpm.get_polynomial(0,0,0,0), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial(1,1,1), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial(2,2,2,1,1,1), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ("a",): 0.2, ("a","b"): 0.3}, "SPIN")
        self.assertAlmostEqual(bpm.get_polynomial("a","a","a","a"), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial("a","a","a"), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial("b","b","b","a","a","a"), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1),): 0.2, ((1,1),(2,2)): 0.3}, "SPIN")
        self.assertAlmostEqual(bpm.get_polynomial((0,0),(0,0),(0,0),(0,0)), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1),(1,1),(1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2),(2,2),(2,2),(1,1),(1,1),(1,1)), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1,1),): 0.2, ((1,1,1),(2,2,2)): 0.3}, "SPIN")
        self.assertAlmostEqual(bpm.get_polynomial((0,0,0),(0,0,0),(0,0,0),(0,0,0)), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1,1),(1,1,1),(1,1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2,2),(2,2,2),(2,2,2),(1,1,1),(1,1,1),(1,1,1)), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1,1,1),): 0.2, ((1,1,1,1),(2,2,2,2)): 0.3}, "SPIN")
        self.assertAlmostEqual(bpm.get_polynomial((0,0,0,0),(0,0,0,0),(0,0,0,0),(0,0,0,0)), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1,1,1),(1,1,1,1),(1,1,1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2,2,2),(2,2,2,2),(2,2,2,2),(1,1,1,1),(1,1,1,1),(1,1,1,1)), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, (1,): 0.2, (1,2): 0.3}, "BINARY")
        self.assertAlmostEqual(bpm.get_polynomial(()), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial(1,1), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial(2,2,2,1,1), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ("a",): 0.2, ("a","b"): 0.3}, "BINARY")
        self.assertAlmostEqual(bpm.get_polynomial(()), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial("a","a"), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial("b","b","b","a","a"), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1),): 0.2, ((1,1),(2,2)): 0.3}, "BINARY")
        self.assertAlmostEqual(bpm.get_polynomial(()), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1),(1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2),(2,2),(2,2),(1,1),(1,1)), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1,1),): 0.2, ((1,1,1),(2,2,2)): 0.3}, "BINARY")
        self.assertAlmostEqual(bpm.get_polynomial(()), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1,1),(1,1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2,2),(2,2,2),(2,2,2),(1,1,1),(1,1,1)), 0.3)

        bpm = cimod.BinaryPolynomialModel({(): 0.1, ((1,1,1,1),): 0.2, ((1,1,1,1),(2,2,2,2)): 0.3}, "BINARY")
        self.assertAlmostEqual(bpm.get_polynomial(()), 0.1)
        self.assertAlmostEqual(bpm.get_polynomial((1,1,1,1),(1,1,1,1)), 0.2)
        self.assertAlmostEqual(bpm.get_polynomial((2,2,2,2),(2,2,2,2),(2,2,2,2),(1,1,1,1),(1,1,1,1)), 0.3)

    def test_add_interactions_from_bpm_dict(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, "SPIN").empty("SPIN")
        bpm.add_interactions_from(self.poly)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, "SPIN").empty("SPIN")
        bpm.add_interactions_from(self.poly_str)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, "SPIN").empty("SPIN")
        bpm.add_interactions_from(self.poly_tuple2)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, "SPIN").empty("SPIN")
        bpm.add_interactions_from(self.poly_tuple3)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, "SPIN").empty("SPIN")
        bpm.add_interactions_from(self.poly_tuple4)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_add_interactions_from_bpm_keyvalues(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, "SPIN").empty("SPIN")
        bpm.add_interactions_from(list(self.poly.keys()), list(self.poly.values()))
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, "SPIN").empty("SPIN")
        bpm.add_interactions_from(list(self.poly_str.keys()), list(self.poly_str.values()))
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, "SPIN").empty("SPIN")
        bpm.add_interactions_from(list(self.poly_tuple2.keys()), list(self.poly_tuple2.values()))
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, "SPIN").empty("SPIN")
        bpm.add_interactions_from(list(self.poly_tuple3.keys()), list(self.poly_tuple3.values()))
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, "SPIN").empty("SPIN")
        bpm.add_interactions_from(list(self.poly_tuple4.keys()), list(self.poly_tuple4.values()))
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_add_offset_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.add_offset(3.0)
        self.assertAlmostEqual(bpm.get_offset(), 3.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 3.0)
        bpm.add_interaction((), 3.0)
        self.assertAlmostEqual(bpm.get_offset(), 6.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 6.0)
        bpm.add_offset(-6.0)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.add_offset(3.0)
        self.assertAlmostEqual(bpm.get_offset(), 3.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 3.0)
        bpm.add_interaction((), 3.0)
        self.assertAlmostEqual(bpm.get_offset(), 6.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 6.0)
        bpm.add_offset(-6.0)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.add_offset(3.0)
        self.assertAlmostEqual(bpm.get_offset(), 3.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 3.0)
        bpm.add_interaction((), 3.0)
        self.assertAlmostEqual(bpm.get_offset(), 6.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 6.0)
        bpm.add_offset(-6.0)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.add_offset(3.0)
        self.assertAlmostEqual(bpm.get_offset(), 3.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 3.0)
        bpm.add_interaction((), 3.0)
        self.assertAlmostEqual(bpm.get_offset(), 6.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 6.0)
        bpm.add_offset(-6.0)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.add_offset(3.0)
        self.assertAlmostEqual(bpm.get_offset(), 3.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 3.0)
        bpm.add_interaction((), 3.0)
        self.assertAlmostEqual(bpm.get_offset(), 6.0)
        self.assertAlmostEqual(bpm.get_polynomial(()), 6.0)
        bpm.add_offset(-6.0)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_remove_interaction_bpm_basic(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.add_interaction((11, 12, 14), -1.0)
        bpm.add_interaction((-7,)       , -2.0)
        bpm.add_interaction((2, 11)     , -3.0)
        bpm.add_interaction(()          , -4.0)
        self.assertAlmostEqual(bpm.get_polynomial(11, 14, 12), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial((-7,))     , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([11, 2])   , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])        , -4.0)
        bpm.remove_interaction((11, 12, 14))
        bpm.remove_interaction(-7)
        bpm.remove_interaction([2, 11])
        bpm.remove_interaction([])
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.add_interaction(("œ∑®", "≈ç", "≥≤µ"), -1.0)
        bpm.add_interaction(("A",)              , -2.0)
        bpm.add_interaction(("¡", "∆")          , -3.0)
        bpm.add_interaction(()                  , -4.0)
        self.assertAlmostEqual(bpm.get_polynomial("œ∑®", "≈ç", "≥≤µ"), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(("A",))            , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial(["¡", "∆"])        , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                , -4.0)
        bpm.remove_interaction(("œ∑®", "≈ç", "≥≤µ"))
        bpm.remove_interaction("A")
        bpm.remove_interaction(["¡", "∆"])
        bpm.remove_interaction([])
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.add_interaction(((11,11), (12,12), (14,14)), -1.0)
        bpm.add_interaction(((-7,-7),)                 , -2.0)
        bpm.add_interaction(((2,2), (11,11))           , -3.0)
        bpm.add_interaction(()                         , -4.0)
        self.assertAlmostEqual(bpm.get_polynomial((11,11), (12,12), (14,14)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7),))               , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2), (11,11)])         , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                       , -4.0)
        bpm.remove_interaction(((11,11), (12,12), (14,14)))
        bpm.remove_interaction((-7,-7))
        bpm.remove_interaction([(2,2), (11,11)])
        bpm.remove_interaction([])
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.add_interaction(((11,11,11), (12,12,12), (14,14,14)), -1.0)
        bpm.add_interaction(((-7,-7,-7),)                       , -2.0)
        bpm.add_interaction(((2,2,2), (11,11,11))               , -3.0)
        bpm.add_interaction(()                                  , -4.0)
        self.assertAlmostEqual(bpm.get_polynomial((11,11,11), (12,12,12), (14,14,14)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7,-7),))                     , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2,2), (11,11,11)])             , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                                , -4.0)
        bpm.remove_interaction(((11,11,11), (12,12,12), (14,14,14)))
        bpm.remove_interaction((-7,-7,-7))
        bpm.remove_interaction([(2,2,2), (11,11,11)])
        bpm.remove_interaction([])
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.add_interaction(((11,11,11,123456789), (12,12,12,123456789), (14,14,14,123456789)), -1.0)
        bpm.add_interaction(((-7,-7,-7,123456789),)                       , -2.0)
        bpm.add_interaction(((2,2,2,123456789), (11,11,11,123456789))               , -3.0)
        bpm.add_interaction(()                                  , -4.0)
        self.assertAlmostEqual(bpm.get_polynomial((14,14,14,123456789), (11,11,11,123456789), (12,12,12,123456789)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7,-7,123456789),))                     , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2,2,123456789), (11,11,11,123456789)])             , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                                , -4.0)
        bpm.remove_interaction(((11,11,11,123456789), (14,14,14,123456789), (12,12,12,123456789)))
        bpm.remove_interaction((-7,-7,-7,123456789))
        bpm.remove_interaction([(2,2,2,123456789), (11,11,11,123456789)])
        bpm.remove_interaction([])
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_remove_interaction_bpm_remove_all(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        for k in self.poly.keys():
            bpm.remove_interaction(k)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        for k in self.poly_str.keys():
            bpm.remove_interaction(*k)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        for k in self.poly_tuple2.keys():
            bpm.remove_interaction(list(k))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        for k in self.poly_tuple3.keys():
            bpm.remove_interaction(k)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        for k in self.poly_tuple4.keys():
            bpm.remove_interaction(*k)
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_remove_interaction_bpm_self_loop(self):
        def check(poly, vartype):
            bpm = cimod.BinaryPolynomialModel(poly, vartype)
            for k in poly.keys():
                key = []
                for i in k:
                    key.append(i)
                for i in k:
                    key.append(i)
                for i in k:
                    key.append(i)
                bpm.remove_interaction(key)
            self.state_test_bpm_empty(bpm, vartype)

        check(self.poly, cimod.SPIN)
        check(self.poly_str, cimod.SPIN)
        check(self.poly_tuple2, cimod.SPIN)
        check(self.poly_tuple3, cimod.SPIN)
        check(self.poly_tuple4, cimod.SPIN)
        check(self.poly, cimod.BINARY)
        check(self.poly_str, cimod.BINARY)
        check(self.poly_tuple2, cimod.BINARY)
        check(self.poly_tuple3, cimod.BINARY)
        check(self.poly_tuple4, cimod.BINARY)

    def test_remove_interactions_from_bpm_basic(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        added_J = {(11, 12, 14): -1.0, (-7,): -2.0, (2, 11): -3.0, (): -4.0}
        bpm.add_interactions_from(added_J)
        self.assertAlmostEqual(bpm.get_polynomial(11, 14, 12), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial((-7,))     , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([11, 2])   , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])        , -4.0)
        bpm.remove_interactions_from(list(added_J))
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        added_J = {("œ∑®", "≈ç", "≥≤µ"): -1.0, ("A",): -2.0, ("¡", "∆"): -3.0, (): -4.0}
        bpm.add_interactions_from(added_J)
        self.assertAlmostEqual(bpm.get_polynomial("œ∑®", "≈ç", "≥≤µ"), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(("A",))            , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial(["¡", "∆"])        , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                , -4.0)
        bpm.remove_interactions_from(list(added_J))
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        added_J = {((11,11), (12,12), (14,14)): -1.0, ((-7,-7),): -2.0, ((2,2), (11,11)): -3.0, (): -4.0}
        bpm.add_interactions_from(added_J)
        self.assertAlmostEqual(bpm.get_polynomial((11,11), (12,12), (14,14)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7),))               , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2), (11,11)])         , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                       , -4.0)
        bpm.remove_interactions_from(list(added_J))
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        added_J = {((11,11,11), (12,12,12), (14,14,14)): -1.0, ((-7,-7,-7),): -2.0, ((2,2,2), (11,11,11)): -3.0, (): -4.0}
        bpm.add_interactions_from(added_J)
        self.assertAlmostEqual(bpm.get_polynomial((11,11,11), (12,12,12), (14,14,14)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7,-7),))                     , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2,2), (11,11,11)])             , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                                , -4.0)
        bpm.remove_interactions_from(list(added_J))
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        added_J = {((11,11,11,123456789), (12,12,12,123456789), (14,14,14,123456789)): -1.0, ((-7,-7,-7,123456789),): -2.0, ((2,2,2,123456789), (11,11,11,123456789)): -3.0, (): -4.0}
        bpm.add_interactions_from(added_J)
        self.assertAlmostEqual(bpm.get_polynomial((14,14,14,123456789), (11,11,11,123456789), (12,12,12,123456789)), -1.0)
        self.assertAlmostEqual(bpm.get_polynomial(((-7,-7,-7,123456789),))                                         , -2.0)
        self.assertAlmostEqual(bpm.get_polynomial([(2,2,2,123456789), (11,11,11,123456789)])                       , -3.0)
        self.assertAlmostEqual(bpm.get_polynomial([])                                                              , -4.0)
        bpm.remove_interactions_from(list(added_J))
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_remove_interactions_from_bpm_remove_all(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.remove_interactions_from(list(self.poly.keys()))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.remove_interactions_from(tuple(self.poly_str.keys()))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.remove_interactions_from(*list(self.poly_tuple2.keys()))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.remove_interactions_from(*tuple(self.poly_tuple3.keys()))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.remove_interactions_from(list(self.poly_tuple4.keys()))
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_remove_offset_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.add_offset(100)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.remove_offset()
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.add_offset(100)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.remove_offset()
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.add_offset(100)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.remove_offset()
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.add_offset(100)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.remove_offset()
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.add_offset(100)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.remove_offset()
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_energy_bpm(self):
        #Spin
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.SPIN).energy(self.spins)       , calculate_bpm_energy(self.poly       , self.spins)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.SPIN).energy(self.spins_str)   , calculate_bpm_energy(self.poly_str   , self.spins_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN).energy(self.spins_tuple2), calculate_bpm_energy(self.poly_tuple2, self.spins_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN).energy(self.spins_tuple3), calculate_bpm_energy(self.poly_tuple3, self.spins_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN).energy(self.spins_tuple4), calculate_bpm_energy(self.poly_tuple4, self.spins_tuple4))

        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.SPIN).energy(list(self.spins.values()))       , calculate_bpm_energy(self.poly       , self.spins)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.SPIN).energy(list(self.spins_str.values()))   , calculate_bpm_energy(self.poly_str   , self.spins_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN).energy(list(self.spins_tuple2.values())), calculate_bpm_energy(self.poly_tuple2, self.spins_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN).energy(list(self.spins_tuple3.values())), calculate_bpm_energy(self.poly_tuple3, self.spins_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN).energy(list(self.spins_tuple4.values())), calculate_bpm_energy(self.poly_tuple4, self.spins_tuple4))

        #Binary
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.BINARY).energy(self.binaries)       , calculate_bpm_energy(self.poly       , self.binaries)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.BINARY).energy(self.binaries_str)   , calculate_bpm_energy(self.poly_str   , self.binaries_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY).energy(self.binaries_tuple2), calculate_bpm_energy(self.poly_tuple2, self.binaries_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY).energy(self.binaries_tuple3), calculate_bpm_energy(self.poly_tuple3, self.binaries_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY).energy(self.binaries_tuple4), calculate_bpm_energy(self.poly_tuple4, self.binaries_tuple4))

        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.BINARY).energy(list(self.binaries.values()))       , calculate_bpm_energy(self.poly       , self.binaries)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.BINARY).energy(list(self.binaries_str.values()))   , calculate_bpm_energy(self.poly_str   , self.binaries_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY).energy(list(self.binaries_tuple2.values())), calculate_bpm_energy(self.poly_tuple2, self.binaries_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY).energy(list(self.binaries_tuple3.values())), calculate_bpm_energy(self.poly_tuple3, self.binaries_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY).energy(list(self.binaries_tuple4.values())), calculate_bpm_energy(self.poly_tuple4, self.binaries_tuple4))

    def test_energies_bpm(self):
        #Spin
        spins_list = [self.spins, self.spins, self.spins, self.spins]
        anser_list = [calculate_bpm_energy(self.poly, self.spins) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly, cimod.SPIN).energies(spins_list), anser_list)

        spins_list = [self.spins_str, self.spins_str, self.spins_str, self.spins_str]
        anser_list = [calculate_bpm_energy(self.poly_str, self.spins_str) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN).energies(spins_list), anser_list)

        spins_list = [self.spins_tuple2, self.spins_tuple2, self.spins_tuple2, self.spins_tuple2]
        anser_list = [calculate_bpm_energy(self.poly_tuple2, self.spins_tuple2) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN).energies(spins_list), anser_list)

        spins_list = [self.spins_tuple3, self.spins_tuple3, self.spins_tuple3, self.spins_tuple3]
        anser_list = [calculate_bpm_energy(self.poly_tuple3, self.spins_tuple3) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN).energies(spins_list), anser_list)

        spins_list = [self.spins_tuple4, self.spins_tuple4, self.spins_tuple4, self.spins_tuple4]
        anser_list = [calculate_bpm_energy(self.poly_tuple4, self.spins_tuple4) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN).energies(spins_list), anser_list)

        #Binary
        binaries_list = [self.binaries, self.binaries, self.binaries, self.binaries]
        anser_list = [calculate_bpm_energy(self.poly, self.binaries) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly, cimod.BINARY).energies(binaries_list), anser_list)

        binaries_list = [self.binaries_str, self.binaries_str, self.binaries_str, self.binaries_str]
        anser_list = [calculate_bpm_energy(self.poly_str, self.binaries_str) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_str, cimod.BINARY).energies(binaries_list), anser_list)

        binaries_list = [self.binaries_tuple2, self.binaries_tuple2, self.binaries_tuple2, self.binaries_tuple2]
        anser_list = [calculate_bpm_energy(self.poly_tuple2, self.binaries_tuple2) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY).energies(binaries_list), anser_list)

        binaries_list = [self.binaries_tuple3, self.binaries_tuple3, self.binaries_tuple3, self.binaries_tuple3]
        anser_list = [calculate_bpm_energy(self.poly_tuple3, self.binaries_tuple3) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY).energies(binaries_list), anser_list)

        binaries_list = [self.binaries_tuple4, self.binaries_tuple4, self.binaries_tuple4, self.binaries_tuple4]
        anser_list = [calculate_bpm_energy(self.poly_tuple4, self.binaries_tuple4) for _ in range(4)]
        self.assertListEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY).energies(binaries_list), anser_list)

    def test_scale_bpm_all_scaled(self):
        d = {}
        for k, v in self.poly.items():
            d[k] = 2*v
        bpm = cimod.BinaryPolynomialModel(d, cimod.SPIN)
        bpm.scale(0.5)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        d = {}
        for k, v in self.poly_str.items():
            d[k] = 2*v
        bpm = cimod.BinaryPolynomialModel(d, cimod.SPIN)
        bpm.scale(0.5)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        d = {}
        for k, v in self.poly_tuple2.items():
            d[k] = 2*v
        bpm = cimod.BinaryPolynomialModel(d, cimod.SPIN)
        bpm.scale(0.5)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        d = {}
        for k, v in self.poly_tuple3.items():
            d[k] = 2*v
        bpm = cimod.BinaryPolynomialModel(d, cimod.SPIN)
        bpm.scale(0.5)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        d = {}
        for k, v in self.poly_tuple4.items():
            d[k] = 2*v
        bpm = cimod.BinaryPolynomialModel(d, cimod.SPIN)
        bpm.scale(0.5)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_scale_bpm_ignored_interaction(self):
        bpm = cimod.BinaryPolynomialModel({k: 2*v for k, v in self.poly.items()}, cimod.SPIN)
        bpm.scale(0.5, ((1, 2), [2, 3, 4]))
        self.assertAlmostEqual(bpm.get_polynomial(1, 2)   , 12.0*2 )
        self.assertAlmostEqual(bpm.get_polynomial(2, 3, 4), 234.0*2)
        bpm.add_interaction((1, 2)   , -12 )
        bpm.add_interaction([2, 3, 4], -234)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: 2*v for k, v in self.poly_str.items()}, cimod.SPIN)
        bpm.scale(0.5, (("a", "b"), ["b", "c", "d"]))
        self.assertAlmostEqual(bpm.get_polynomial("a", "b")     , 12.0*2 )
        self.assertAlmostEqual(bpm.get_polynomial("b", "c", "d"), 234.0*2)
        bpm.add_interaction(("a", "b")     , -12)
        bpm.add_interaction(["b", "c", "d"], -234)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: 2*v for k, v in self.poly_tuple2.items()}, cimod.SPIN)
        bpm.scale(0.5, (((1, 1), (2, 2)), [(2, 2), (3, 3), (4, 4)]))
        self.assertAlmostEqual(bpm.get_polynomial((2, 2), (1, 1))        , 12.0*2 )
        self.assertAlmostEqual(bpm.get_polynomial((4, 4), (3, 3), (2, 2)), 234.0*2)
        bpm.add_interaction([(2, 2), (1, 1)]        , -12 )
        bpm.add_interaction([(4, 4), (3, 3), (2, 2)], -234)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: 2*v for k, v in self.poly_tuple3.items()}, cimod.SPIN)
        bpm.scale(0.5, (((1, 1, 1), (2, 2, 2)), [(2, 2, 2), (3, 3, 3), (4, 4, 4)]))
        self.assertAlmostEqual(bpm.get_polynomial((2, 2, 2), (1, 1, 1))           , 12.0*2 )
        self.assertAlmostEqual(bpm.get_polynomial((4, 4, 4), (3, 3, 3), (2, 2, 2)), 234.0*2)
        bpm.add_interaction([(2, 2, 2), (1, 1, 1)]           , -12 )
        bpm.add_interaction([(4, 4, 4), (3, 3, 3), (2, 2, 2)], -234)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel({k: 2*v for k, v in self.poly_tuple4.items()}, cimod.SPIN)
        bpm.scale(0.5, (((1, 1, 1, 1), (2, 2, 2, 2)), [(2, 2, 2, 2), (3, 3, 3, 3), (4, 4, 4, 4)]))
        self.assertAlmostEqual(bpm.get_polynomial((2, 2, 2, 2), (1, 1, 1, 1))              , 12.0*2 )
        self.assertAlmostEqual(bpm.get_polynomial((4, 4, 4, 4), (3, 3, 3, 3), (2, 2, 2, 2)), 234.0*2)
        bpm.add_interaction([(2, 2, 2, 2), (1, 1, 1, 1)]              , -12 )
        bpm.add_interaction([(4, 4, 4, 4), (3, 3, 3, 3), (2, 2, 2, 2)], -234)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_scale_bpm_ignored_offset(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.add_offset(100)
        bpm.scale(0.5, ((),))
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, ignored_offset = True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, [()], True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.add_offset(100)
        bpm.scale(0.5, ((),))
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, ignored_offset = True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, [()], True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.add_offset(100)
        bpm.scale(0.5, ((),))
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, ignored_offset = True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, [()], True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.add_offset(100)
        bpm.scale(0.5, ((),))
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, ignored_offset = True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, [()], True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.add_offset(100)
        bpm.scale(0.5, ((),))
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, ignored_offset = True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        bpm.scale(0.5, [()], True)
        self.assertAlmostEqual(bpm.get_polynomial(()), 100)
        
    def test_normalize_bpm_all_normalize(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)        
        bpm.normalize((-1, +1))
        bpm.scale(234.0)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)        
        bpm.normalize((-1, +1))
        bpm.scale(234.0)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)        
        bpm.normalize((-1, +1))
        bpm.scale(234.0)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)        
        bpm.normalize((-1, +1))
        bpm.scale(234.0)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)        
        bpm.normalize((-1, +1))
        bpm.scale(234.0)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_normalize_bpm_ignored_interaction(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.normalize((-1, 1), list(self.poly))
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.normalize((-1, 1), list(self.poly_str))
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.normalize((-1, 1), list(self.poly_tuple2))
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.normalize((-1, 1), list(self.poly_tuple3))
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.normalize((-1, 1), list(self.poly_tuple4))
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

    def test_serializable_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm(bpm_from, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm(bpm_from, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm(bpm_from, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm(bpm_from, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm(bpm_from, self.poly_tuple4, cimod.SPIN)

    def test_serializable_bpm_empty(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.clear()
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm_empty(bpm_from, cimod.SPIN)
        self.assertEqual(bpm_from.index_type, bpm.index_type)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.clear()
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm_empty(bpm_from, cimod.SPIN)
        self.assertEqual(bpm_from.index_type, bpm.index_type)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.clear()
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm_empty(bpm_from, cimod.SPIN)
        self.assertEqual(bpm_from.index_type, bpm.index_type)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.clear()
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm_empty(bpm_from, cimod.SPIN)
        self.assertEqual(bpm_from.index_type, bpm.index_type)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.clear()
        bpm_from = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.state_test_bpm_empty(bpm_from, cimod.SPIN)
        self.assertEqual(bpm_from.index_type, bpm.index_type)

    def test_from_hubo_bpm_from_dict(self):
        bpm = cimod.BinaryPolynomialModel.from_hubo(self.poly)
        self.state_test_bpm(bpm, self.poly, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(self.poly_str)
        self.state_test_bpm(bpm, self.poly_str, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple2)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple3)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple4)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo({})
        self.state_test_bpm_empty(bpm, cimod.BINARY)

    def test_from_hubo_bpm_from_key_value(self):
        bpm = cimod.BinaryPolynomialModel.from_hubo(list(self.poly.keys()), list(self.poly.values()))
        self.state_test_bpm(bpm, self.poly, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(list(self.poly_str.keys()), list(self.poly_str.values()))
        self.state_test_bpm(bpm, self.poly_str, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(list(self.poly_tuple2.keys()), list(self.poly_tuple2.values()))
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(list(self.poly_tuple3.keys()), list(self.poly_tuple3.values()))
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo(list(self.poly_tuple4.keys()), list(self.poly_tuple4.values()))
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel.from_hubo([], [])
        self.state_test_bpm_empty(bpm, cimod.BINARY)

    def test_from_hising_bpm_from_dict(self):
        bpm = cimod.BinaryPolynomialModel.from_hising(self.poly)
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(self.poly_str)
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(self.poly_tuple2)
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(self.poly_tuple3)
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(self.poly_tuple4)
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising({})
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_from_hising_bpm_from_key_value(self):
        bpm = cimod.BinaryPolynomialModel.from_hising(list(self.poly.keys()), list(self.poly.values()))
        self.state_test_bpm(bpm, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(list(self.poly_str.keys()), list(self.poly_str.values()))
        self.state_test_bpm(bpm, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(list(self.poly_tuple2.keys()), list(self.poly_tuple2.values()))
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(list(self.poly_tuple3.keys()), list(self.poly_tuple3.values()))
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising(list(self.poly_tuple4.keys()), list(self.poly_tuple4.values()))
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel.from_hising([], [])
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_clear_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm.clear()
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm.clear()
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm.clear()
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm.clear()
        self.state_test_bpm_empty(bpm, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm.clear()
        self.state_test_bpm_empty(bpm, cimod.SPIN)

    def test_to_hubo_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        J_hubo = bpm.to_hubo()
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, J_hubo, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        J_hubo = bpm.to_hubo()
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, J_hubo, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        J_hubo = bpm.to_hubo()
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, J_hubo, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        J_hubo = bpm.to_hubo()
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, J_hubo, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        J_hubo = bpm.to_hubo()
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, J_hubo, cimod.BINARY)

    def test_to_hising_bpm(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.BINARY)
        J_ising = bpm.to_hising()
        bpm.change_vartype("SPIN")
        self.state_test_bpm(bpm, J_ising, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.BINARY)
        J_ising = bpm.to_hising()
        bpm.change_vartype("SPIN")
        self.state_test_bpm(bpm, J_ising, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY)
        J_ising = bpm.to_hising()
        bpm.change_vartype("SPIN")
        self.state_test_bpm(bpm, J_ising, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY)
        J_ising = bpm.to_hising()
        bpm.change_vartype("SPIN")
        self.state_test_bpm(bpm, J_ising, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY)
        J_ising = bpm.to_hising()
        bpm.change_vartype("SPIN")
        self.state_test_bpm(bpm, J_ising, cimod.SPIN)

    def test_change_vartype_bpm_spin_binary_spin(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm_binary = bpm.change_vartype("BINARY", False)
        bpm_ising  = bpm_binary.change_vartype("SPIN", False)
        self.state_test_bpm(bpm_ising, self.poly, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm_binary = bpm.change_vartype("BINARY", False)
        bpm_ising  = bpm_binary.change_vartype("SPIN", False)
        self.state_test_bpm(bpm_ising, self.poly_str, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm_binary = bpm.change_vartype("BINARY", False)
        bpm_ising  = bpm_binary.change_vartype("SPIN", False)
        self.state_test_bpm(bpm_ising, self.poly_tuple2, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm_binary = bpm.change_vartype("BINARY", False)
        bpm_ising  = bpm_binary.change_vartype("SPIN", False)
        self.state_test_bpm(bpm_ising, self.poly_tuple3, cimod.SPIN)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm_binary = bpm.change_vartype("BINARY", False)
        bpm_ising  = bpm_binary.change_vartype("SPIN", False)
        self.state_test_bpm(bpm_ising, self.poly_tuple4, cimod.SPIN)

    def test_change_vartype_bpm_binary_spin_binary(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.BINARY)
        bpm.change_vartype("SPIN")
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, self.poly, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.BINARY)
        bpm.change_vartype("SPIN")
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, self.poly_str, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY)
        bpm.change_vartype("SPIN")
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, self.poly_tuple2, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY)
        bpm.change_vartype("SPIN")
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, self.poly_tuple3, cimod.BINARY)

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY)
        bpm.change_vartype("SPIN")
        bpm.change_vartype("BINARY")
        self.state_test_bpm(bpm, self.poly_tuple4, cimod.BINARY)

class UtilTest(unittest.TestCase):
    def test_convert_sample(self):

        #BQM
        strh = {'b': -2, 'a': 1}
        strJ = {('c', 'd'): 0.5, ('a', 'b'): -1, ('b', 'c'): -3}
        strspins = {'a': 1, 'c': 1, 'd': 1, 'b': -1}
        raw_strspins = [1,-1,1,1]
        raw_strspins_binaries = [1,0,1,1]

        strQ = {('a', 'a'): 1, ('c', 'a'): -0.2, ('b', 'd'): 3, ('b', 'c'): -1}
        strbinaries = {'a': 0, 'b': 1, 'c': 1, 'd': 0}
        raw_strbinaries = [0,1,1,0]
        raw_strbinaries_spins = [-1,1,1,-1]

        for (sparse, mat_type) in [(True, csr_matrix), (False, np.ndarray)]:
            bqm = cimod.model.BinaryQuadraticModel(strh, strJ, 'SPIN', sparse=sparse)

            #spins
            true_ising_e = calculate_ising_energy(strh, strJ, strspins)
            state, energy = cimod.utils.get_state_and_energy(bqm, raw_strspins)
            self.assertDictEqual(state, strspins)
            self.assertEqual(energy, true_ising_e)

            #spins with raw binaries
            state, energy = cimod.utils.get_state_and_energy(bqm, raw_strspins_binaries, offset=10)
            self.assertDictEqual(state, strspins)
            self.assertEqual(energy, true_ising_e + 10)

            #binaries
            bqm = cimod.model.BinaryQuadraticModel.from_qubo(strQ, sparse=sparse)
            true_qubo_e = calculate_qubo_energy(strQ, strbinaries)
            state, energy = cimod.utils.get_state_and_energy(bqm, raw_strbinaries)
            self.assertDictEqual(state, strbinaries)
            self.assertEqual(energy, true_qubo_e)

            #binaries with raw spins
            state, energy = cimod.utils.get_state_and_energy(bqm, raw_strbinaries_spins, offset=12)
            self.assertDictEqual(state, strbinaries)
            self.assertEqual(energy, true_qubo_e + 12)

        #BPM
        poly_str     = {("a",):1.0, ("c",):3.0, ("a","b"):12.0, ("a","c"):13.0, ("b","c","d"):234.0, ("c","e"):35.0}
        spins_str    = {"a":+1, "b":-1, "c":+1, "d":-1, "e":+1} 
        binaries_str = {"a": 1, "b": 0, "c": 1, "d": 0, "e": 1}
        raw_spins    = [1,-1,1,-1,1]
        raw_binaries = [1,0,1,0,1]

        true_ising_e = calculate_bpm_energy(poly_str, spins_str)
        true_hubo_e = calculate_bpm_energy(poly_str, binaries_str)

        #spins
        bpm = cimod.model.BinaryPolynomialModel(poly_str, 'SPIN')
        state, energy = cimod.utils.get_state_and_energy(bpm, raw_spins)
        self.assertDictEqual(state, spins_str)
        self.assertEqual(energy, true_ising_e)

        #spins with raw binaries
        state, energy = cimod.utils.get_state_and_energy(bpm, raw_binaries, offset=20)
        self.assertDictEqual(state, spins_str)
        self.assertEqual(energy, true_ising_e + 20)

        #binaries
        bpm = cimod.model.BinaryPolynomialModel(poly_str, 'BINARY')
        state, energy = cimod.utils.get_state_and_energy(bpm, raw_binaries)
        self.assertDictEqual(state, binaries_str)
        self.assertEqual(energy, true_hubo_e)

        #spins with raw binaries
        state, energy = cimod.utils.get_state_and_energy(bpm, raw_spins, offset=30)
        self.assertDictEqual(state, binaries_str)
        self.assertEqual(energy, true_hubo_e + 30)




if __name__ == '__main__':
    unittest.main()
