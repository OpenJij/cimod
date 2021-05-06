import unittest

import numpy as np
import cimod
import cxxcimod
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

        bqm = cimod.BinaryQuadraticModel(self.tupleh, self.tupleJ)
        self.assertEqual(type(bqm.interaction_matrix()), np.ndarray)

        self.assertEqual(bqm.vartype, cimod.SPIN)

        bqm_qubo = cimod.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)
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

        bqm = cimod.BinaryQuadraticModel(self.tupleh, self.tupleJ)
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

        # Test to calculate energy

        # Test Ising energy
        bqm = cimod.BinaryQuadraticModel(self.tupleh, self.tupleJ)
        ising_energy_bqm = bqm.energy(self.tuplespins)
        true_ising_e = calculate_ising_energy(self.tupleh, self.tupleJ, self.tuplespins)
        self.assertEqual(ising_energy_bqm, true_ising_e)

        # Test QUBO energy
        bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)
        qubo_energy_bqm = bqm.energy(self.tuplebinaries)
        true_qubo_e = calculate_qubo_energy(self.tupleQ, self.tuplebinaries)
        self.assertEqual(qubo_energy_bqm, true_qubo_e)

        # QUBO == Ising
        spins = {(0,1,2): 1, (1,2,3): 1, (2,3,4): -1, (3,4,5): 1}
        binary = {(0,1,2): 1, (1,2,3): 1, (2,3,4): 0, (3,4,5): 1}
        qubo_bqm = cimod.BinaryQuadraticModel.from_qubo(Q=self.tupleQ)

        qubo_energy = qubo_bqm.energy(binary)

        self.assertEqual(qubo_energy, qubo_bqm.energy(spins, convert_sample=True))

    def test_change_vartype(self):
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        self.assertEqual(bqm.vartype, cimod.SPIN)
        bqm2 = bqm.change_vartype('BINARY')
        self.assertEqual(bqm.vartype, cimod.BINARY)
        self.assertEqual(bqm.linear, bqm2.linear)
        self.assertEqual(bqm.quadratic, bqm2.quadratic)
        self.assertEqual(bqm.offset, bqm2.offset)
        self.assertEqual(bqm.vartype, bqm2.vartype)
        bqm.change_vartype('SPIN', inplace=False)
        self.assertEqual(bqm.vartype, cimod.BINARY)

    def test_serializable(self):
        bqm = cimod.BinaryQuadraticModel(self.h, self.J)
        serial = bqm.to_serializable()
        decode_bqm = cimod.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

        bqm = cimod.BinaryQuadraticModel(self.strh, self.strJ)
        serial = bqm.to_serializable()
        decode_bqm = cimod.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

        bqm = cimod.BinaryQuadraticModel(self.tupleh, self.tupleJ)
        serial = bqm.to_serializable()
        decode_bqm = cimod.BinaryQuadraticModel.from_serializable(serial)
        self.assertEqual(bqm.linear, decode_bqm.linear)
        self.assertEqual(bqm.quadratic, decode_bqm.quadratic)
        self.assertEqual(bqm.offset, decode_bqm.offset)
        self.assertEqual(bqm.vartype, decode_bqm.vartype)

    def test_serializable_consistent_with_dimod(self):
        for (_from,_to) in [(dimod, cimod), (cimod, dimod)]:
            bqm = _from.BinaryQuadraticModel(self.h, self.J, vartype='SPIN')
            serial = bqm.to_serializable()
            decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            # order of indices in quadratic is not considered.
            self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)

            bqm = _from.BinaryQuadraticModel(self.strh, self.strJ, vartype='SPIN')
            serial = bqm.to_serializable()
            decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            # order of indices in quadratic is not considered.
            self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)

            bqm = _from.BinaryQuadraticModel(self.tupleh, self.tupleJ, vartype='SPIN')
            serial = bqm.to_serializable()
            decode_bqm = _to.BinaryQuadraticModel.from_serializable(serial)
            self.assertEqual(bqm.linear, decode_bqm.linear)
            # order of indices in quadratic is not considered.
            self.assertEqual({(min(k), max(k)):v for k,v in bqm.quadratic.items()}, {(min(k), max(k)):v for k,v in decode_bqm.quadratic.items()})
            self.assertEqual(bqm.offset, decode_bqm.offset)
            self.assertEqual(bqm.vartype, decode_bqm.vartype)


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
        self.assertEqual(bpm.vartype, vartype)
        self.assertEqual(bpm.num_interactions, len(poly))
        self.assertEqual(bpm.num_variables, len(set(j for i in poly.keys() for j in i)))
        self.assertEqual(bpm.degree, max([len(i) for i in poly.keys()]))
        self.assertDictEqual(bpm.get_polynomial(), poly)
        self.assertSetEqual(bpm.get_variables(), set(j for i in poly.keys() for j in i))
        self.assertAlmostEqual(bpm.get_offset(), poly[()] if tuple() in poly else 0.0)
        for index in poly.keys():
            self.assertAlmostEqual(bpm.get_polynomial(index), poly[index])
            self.assertAlmostEqual(bpm.get_polynomial(list(index)), poly[index])
            self.assertAlmostEqual(bpm.get_polynomial(*index), poly[index])   

    # Test BinaryPolynomialModel constructor
    def test_bpm_constructor(self):
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly       , cimod.SPIN), self.poly       , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_str   , cimod.SPIN), self.poly_str   , cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN), self.poly_tuple2, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN), self.poly_tuple3, cimod.SPIN)
        self.state_test_bpm(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN), self.poly_tuple4, cimod.SPIN)

    #Tese energy calculations
    def test_bpm_calc_energy(self):
        #Spin
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.SPIN).energy(self.spins)       , calculate_bpm_energy(self.poly       , self.spins)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.SPIN).energy(self.spins_str)   , calculate_bpm_energy(self.poly_str   , self.spins_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN).energy(self.spins_tuple2), calculate_bpm_energy(self.poly_tuple2, self.spins_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN).energy(self.spins_tuple3), calculate_bpm_energy(self.poly_tuple3, self.spins_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN).energy(self.spins_tuple4), calculate_bpm_energy(self.poly_tuple4, self.spins_tuple4))

        #Binary
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly       , cimod.BINARY).energy(self.binaries)       , calculate_bpm_energy(self.poly       , self.binaries)       )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_str   , cimod.BINARY).energy(self.binaries_str)   , calculate_bpm_energy(self.poly_str   , self.binaries_str)   )
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY).energy(self.binaries_tuple2), calculate_bpm_energy(self.poly_tuple2, self.binaries_tuple2))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY).energy(self.binaries_tuple3), calculate_bpm_energy(self.poly_tuple3, self.binaries_tuple3))
        self.assertEqual(cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY).energy(self.binaries_tuple4), calculate_bpm_energy(self.poly_tuple4, self.binaries_tuple4))

    

    '''
    # Test BinaryPolynomialModel from_hubo function
    def test_bqm_from_hubo(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.BINARY)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_hubo(self.poly)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()
        
        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.BINARY)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_hubo(self.poly_str)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.BINARY)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple2)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.BINARY)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple3)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.BINARY)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_hubo(self.poly_tuple4)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

    # Test BinaryPolynomialModel from_ising function
    def test_bqm_from_ising(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_ising(self.poly)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()
        
        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_ising(self.poly_str)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_ising(self.poly_tuple2)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_ising(self.poly_tuple3)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        bpm_from_hubo = cimod.BinaryPolynomialModel.from_ising(self.poly_tuple4)
        self.assertEqual    (bpm_from_hubo.vartype         , bpm.vartype)         #vartype
        self.assertEqual    (bpm_from_hubo.num_variables    , bpm.num_variables)    #get_length()
        self.assertSetEqual (bpm_from_hubo.get_variables()       , bpm.get_variables())       #variables
        self.assertSetEqual (bpm_from_hubo.get_variables() , bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm_from_hubo.get_polynomial()      , bpm.get_polynomial())      #polynomial
        self.assertDictEqual(bpm_from_hubo.get_polynomial(), bpm.get_polynomial())#get_polynomial()

    def test_bpm_serializable(self):
        bpm = cimod.BinaryPolynomialModel(self.poly, cimod.SPIN)
        decode_bpm = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.assertEqual    (bpm.vartype         , decode_bpm.vartype)  #vartype
        self.assertEqual    (bpm.num_variables , decode_bpm.num_variables)           #get_length()
        self.assertSetEqual (bpm.get_variables()       , decode_bpm.get_variables()) #variables
        self.assertSetEqual (bpm.get_variables() , decode_bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm.get_polynomial()      , decode_bpm.get_polynomial()) #polynomial
        self.assertDictEqual(bpm.get_polynomial(), decode_bpm.get_polynomial()) #get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_str, cimod.SPIN)
        decode_bpm = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.assertEqual    (bpm.vartype         , decode_bpm.vartype)  #vartype
        self.assertEqual    (bpm.num_variables , decode_bpm.num_variables)           #get_length()
        self.assertSetEqual (bpm.get_variables()       , decode_bpm.get_variables()) #variables
        self.assertSetEqual (bpm.get_variables() , decode_bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm.get_polynomial()      , decode_bpm.get_polynomial()) #polynomial
        self.assertDictEqual(bpm.get_polynomial(), decode_bpm.get_polynomial()) #get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple2, cimod.SPIN)
        decode_bpm = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.assertEqual    (bpm.vartype         , decode_bpm.vartype)  #vartype
        self.assertEqual    (bpm.num_variables , decode_bpm.num_variables)           #get_length()
        self.assertSetEqual (bpm.get_variables()       , decode_bpm.get_variables()) #variables
        self.assertSetEqual (bpm.get_variables() , decode_bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm.get_polynomial()      , decode_bpm.get_polynomial()) #polynomial
        self.assertDictEqual(bpm.get_polynomial(), decode_bpm.get_polynomial()) #get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple3, cimod.SPIN)
        decode_bpm = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.assertEqual    (bpm.vartype         , decode_bpm.vartype)  #vartype
        self.assertEqual    (bpm.num_variables , decode_bpm.num_variables)           #get_length()
        self.assertSetEqual (bpm.get_variables()       , decode_bpm.get_variables()) #variables
        self.assertSetEqual (bpm.get_variables() , decode_bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm.get_polynomial()      , decode_bpm.get_polynomial()) #polynomial
        self.assertDictEqual(bpm.get_polynomial(), decode_bpm.get_polynomial()) #get_polynomial()

        bpm = cimod.BinaryPolynomialModel(self.poly_tuple4, cimod.SPIN)
        decode_bpm = cimod.BinaryPolynomialModel.from_serializable(bpm.to_serializable())
        self.assertEqual    (bpm.vartype         , decode_bpm.vartype)  #vartype
        self.assertEqual    (bpm.num_variables , decode_bpm.num_variables)           #get_length()
        self.assertSetEqual (bpm.get_variables()       , decode_bpm.get_variables()) #variables
        self.assertSetEqual (bpm.get_variables() , decode_bpm.get_variables()) #GetVariables()
        self.assertDictEqual(bpm.get_polynomial()      , decode_bpm.get_polynomial()) #polynomial
        self.assertDictEqual(bpm.get_polynomial(), decode_bpm.get_polynomial()) #get_polynomial()
        '''



if __name__ == '__main__':
    unittest.main()
