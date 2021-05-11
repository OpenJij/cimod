#include "gtest/gtest.h"

#include "../src/binary_quadratic_model.hpp"
#include "../src/binary_polynomial_model.hpp"
#include "../src/binary_quadratic_model_dict.hpp"

#include "test_bqm.hpp"

#include <nlohmann/json.hpp>

#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <tuple>

using json = nlohmann::json;
using namespace cimod;

namespace
{
    TEST(DenseConstructionTest, Construction)
    {
        BQMTester<Dense>::test_DenseConstructionTest_Construction();
        BQMTester<Sparse>::test_DenseConstructionTest_Construction();
        BQMTester<Dict>::test_DenseConstructionTest_Construction();
    }

    TEST(DenseConstructionTest, ConstructionString)
    {
        BQMTester<Dense>::test_DenseConstructionTest_ConstructionString();
        BQMTester<Sparse>::test_DenseConstructionTest_ConstructionString();
        BQMTester<Dict>::test_DenseConstructionTest_ConstructionString();
    }

    TEST(DenseConstructionTest, ConstructionMatrix)
    {
        BQMTester<Dense>::test_DenseConstructionTest_ConstructionMatrix();
        BQMTester<Sparse>::test_DenseConstructionTest_ConstructionMatrix();
    }

    TEST(DenseConstructionTest, ConstructionMatrix2)
    {
        BQMTester<Dense>::test_DenseConstructionTest_ConstructionMatrix2();
        BQMTester<Sparse>::test_DenseConstructionTest_ConstructionMatrix2();
    }

    TEST(DenseBQMFunctionTest, add_variable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_add_variable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_add_variable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_add_variable();
    }

    TEST(DenseBQMFunctionTest, add_variables_from)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_add_variables_from();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_add_variables_from();
        BQMTester<Dict>::test_DenseBQMFunctionTest_add_variables_from();
    }

    TEST(DenseBQMFunctionTest, add_interaction)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_add_interaction();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_add_interaction();
        BQMTester<Dict>::test_DenseBQMFunctionTest_add_interaction();
    }

    TEST(DenseBQMFunctionTest, add_interactions_from)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_add_interactions_from();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_add_interactions_from();
        BQMTester<Dict>::test_DenseBQMFunctionTest_add_interactions_from();
    }

    TEST(DenseBQMFunctionTest, add_offset)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_add_offset();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_add_offset();
        BQMTester<Dict>::test_DenseBQMFunctionTest_add_offset();
    }

    TEST(DenseBQMFunctionTest, energy)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_energy();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_energy();
        BQMTester<Dict>::test_DenseBQMFunctionTest_energy();
    }

    TEST(DenseBQMFunctionTest, energies)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_energies();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_energies();
        BQMTester<Dict>::test_DenseBQMFunctionTest_energies();
    }

    TEST(DenseBQMFunctionTest, empty)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_empty();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_empty();
        BQMTester<Dict>::test_DenseBQMFunctionTest_empty();
    }

    TEST(DenseBQMFunctionTest, to_qubo)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_to_qubo();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_to_qubo();
        BQMTester<Dict>::test_DenseBQMFunctionTest_to_qubo();
    }

    TEST(DenseBQMFunctionTest, to_ising)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_to_ising();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_to_ising();
        BQMTester<Dict>::test_DenseBQMFunctionTest_to_ising();
    }

    TEST(DenseBQMFunctionTest, from_qubo)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_from_qubo();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_from_qubo();
        BQMTester<Dict>::test_DenseBQMFunctionTest_from_qubo();
    }

    TEST(DenseBQMFunctionTest, from_ising)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_from_ising();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_from_ising();
        BQMTester<Dict>::test_DenseBQMFunctionTest_from_ising();
    }

    TEST(DenseBQMFunctionTest, remove_variable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_remove_variable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_remove_variable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_remove_variable();
    }

    TEST(DenseBQMFunctionTest, remove_variables_from)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_remove_variables_from();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_remove_variables_from();
        BQMTester<Dict>::test_DenseBQMFunctionTest_remove_variables_from();
    }

    TEST(DenseBQMFunctionTest, remove_interaction)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_remove_interaction();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_remove_interaction();
        BQMTester<Dict>::test_DenseBQMFunctionTest_remove_interaction();
    }

    TEST(DenseBQMFunctionTest, scale)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_scale();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_scale();
        BQMTester<Dict>::test_DenseBQMFunctionTest_scale();
    }

    TEST(DenseBQMFunctionTest, normalize)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_normalize();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_normalize();
        BQMTester<Dict>::test_DenseBQMFunctionTest_normalize();
    }

    TEST(DenseBQMFunctionTest, fix_variable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_fix_variable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_fix_variable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_fix_variable();
    }

    TEST(DenseBQMFunctionTest, flip_variable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_flip_variable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_flip_variable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_flip_variable();
    }

    TEST(DenseBQMFunctionTest, flip_variable_binary)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_flip_variable_binary();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_flip_variable_binary();
        BQMTester<Dict>::test_DenseBQMFunctionTest_flip_variable_binary();
    }

    TEST(DenseBQMFunctionTest, change_vartype)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_change_vartype();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_change_vartype();
        BQMTester<Dict>::test_DenseBQMFunctionTest_change_vartype();
    }

    TEST(DenseBQMFunctionTest, to_serializable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_to_serializable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_to_serializable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_to_serializable();
    }

    TEST(DenseBQMFunctionTest, from_serializable)
    {
        BQMTester<Dense>::test_DenseBQMFunctionTest_from_serializable();
        BQMTester<Sparse>::test_DenseBQMFunctionTest_from_serializable();
        BQMTester<Dict>::test_DenseBQMFunctionTest_from_serializable();
    }

/*
//google test for binary polynomial model
TEST(ConstructionTestBPM, Construction) {
   
   Polynomial<uint32_t, double> polynomial {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
   
   Vartype vartype = Vartype::BINARY;
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);

   EXPECT_EQ(bpm.get_num_variables(), 4);
   
   //variables
   auto variables = bpm.get_variables();
   EXPECT_EQ(variables[0], 1);
   EXPECT_EQ(variables[1], 2);
   EXPECT_EQ(variables[2], 3);
   EXPECT_EQ(variables[3], 4);
      
   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(vartype, bpm.get_vartype());
   
}


TEST(ConstructionTestBPM, ConstructionString) {
   
   Polynomial<std::string, double> polynomial {
      //linear biases
      {{"a"}, 1.0}, {{"b"}, 2.0}, {{"c"}, 3.0}, {{"d"}, 4.0},
      //quadratic biases
      {{"a", "b"}, 12.0}, {{"a", "c"}, 13.0}, {{"a", "d"}, 14.0},
      {{"b", "c"}, 23.0}, {{"b", "d"}, 24.0},
      {{"c", "d"}, 34.0},
      //polynomial biases
      {{"a", "b", "c"}, 123.0}, {{"a", "b", "d"}, 124.0}, {{"a", "c", "d"}, 134.0},
      {{"b", "c", "d"}, 234.0},
      {{"a", "b", "c", "d"}, 1234.0}
   };
   
   Vartype vartype = Vartype::BINARY;
   
   BinaryPolynomialModel<std::string, double> bpm(polynomial, vartype);
   
   EXPECT_EQ(bpm.length(), 4);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count("a"), 1);
   EXPECT_EQ(bpm.get_variables().count("b"), 1);
   EXPECT_EQ(bpm.get_variables().count("c"), 1);
   EXPECT_EQ(bpm.get_variables().count("d"), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b"}), 12.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "c"}), 13.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "d"}), 14.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "c"}), 123.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "d"}), 124.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "c", "d"}), 134.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "c", "d"}), 1234.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "c"}), 23.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "d"}), 24.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "c", "d"}), 234.0);

   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("c").at({"c", "d"}), 34.0);
   
   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(vartype, bpm.get_vartype());
   
}

TEST(FunctionTestBPM, add_linear) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0},
      {{0, 1}, 0.5},
      {{0, 1, 2}, 1.5}
   };
   
   Vartype vartype = Vartype::SPIN;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   EXPECT_EQ(bpm.length(), 3);
   
   bpm.add_linear(0, 1);
   bpm.add_linear(1, 1);
   bpm.add_linear(2, 1);
   bpm.add_linear(3, 1);
   
   EXPECT_EQ(bpm.length(), 4);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0}), 1.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 2.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 1.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({3}), 1.0);

   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 0.5);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 1.5);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1}), 0.5);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 2}), 1.5);

}

TEST(FunctionTestBPM, add_interaction) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0},
      {{0, 1}, 0.5},
      {{0, 1, 2}, 1.5}
   };
   
   Vartype vartype = Vartype::SPIN;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   EXPECT_EQ(bpm.length(), 3);
   
   bpm.add_interaction({1}, 1.0);
   bpm.add_interaction({4}, 4.0);
   bpm.add_interaction({1, 2}, 12.0);
   bpm.add_interaction({0, 1, 2}, 1.5);
   bpm.add_interaction({0, 1, 2, 3}, 123.0);

   EXPECT_EQ(bpm.length(), 5);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 2.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({3}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({4}), 4.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 0.5);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1, 2}), 12.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 3.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2, 3}), 123.0);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   EXPECT_EQ(bpm.get_variables().count(4), 1);

   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1}), 0.5);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 2}), 3.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 2, 3}), 123.0);

   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2}), 12.0);

}

TEST(FunctionTestBPM, add_interaction_from) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0},
      {{0, 1}, 0.5},
      {{0, 1, 2}, 1.5}
   };
   
   Vartype vartype = Vartype::SPIN;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   Polynomial<uint32_t, double> add_polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 0.5},
      {{0, 1, 2}, 1.5}, {{0, 1, 3}, 13.0}
   };
   
   bpm.add_interactions_from(add_polynomial);
   
   EXPECT_EQ(bpm.length(), 4);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 2.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 2.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({3}), 0.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 1.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 3.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 3}), 13.0);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1}), 1.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 2}), 3.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 3}), 13.0);

}

TEST(FunctionTestBPM, energy_SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   Sample<uint32_t> sample_variables_spin_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_2{{0, +1}, {1, -1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_3{{0, -1}, {1, -1}, {2, -1}};

   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_2), -123.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_3), -81.0 );
      
}
 
TEST(FunctionTestBPM, energy_BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   Vartype vartype = Vartype::BINARY;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};

   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_3), 0.0   );
      
}

TEST(FunctionTestBPM, energies) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   std::vector<Sample<uint32_t>> sample_variables_spin {
      {{0, +1}, {1, +1}, {2, +1}},
      {{0, +1}, {1, -1}, {2, +1}},
      {{0, -1}, {1, -1}, {2, -1}}
   };
   std::vector<double> en_vec = bpm.energies(sample_variables_spin);
   
   EXPECT_DOUBLE_EQ(en_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec[1], -123.0);
   EXPECT_DOUBLE_EQ(en_vec[2], -81.0 );
   
}

TEST(FunctionTestBPM, remove_interaction) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.remove_interaction({0, 2});
   bpm.remove_interaction({0, 1, 2});

   EXPECT_TRUE(bpm.get_polynomial().find({0})       != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({1})       != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({2})       != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 1})    != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 2})    == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({1, 2})    != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 1, 2}) == bpm.get_polynomial().end());
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   
   //Adjacency
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 1}), 1);
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 2}), 0);
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 1, 2}), 0);
   EXPECT_EQ(bpm.get_adjacency().at(1).count({1, 2}), 1);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1}), 11.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2}), 12.0);

}

TEST(FunctionTestBPM, remove_variable) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.remove_variable(1);
   
   EXPECT_TRUE(bpm.get_polynomial().find({0})       != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({1})       == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({2})       != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 1})    == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 2})    != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({1, 2})    == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({0, 1, 2}) == bpm.get_polynomial().end());
   
   EXPECT_EQ(bpm.length(), 2);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 0);
   EXPECT_EQ(bpm.get_variables().count(2), 1);

   //Adjacency
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 1}), 0);
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 2}), 1);
   EXPECT_EQ(bpm.get_adjacency().count(1), 0);
   EXPECT_EQ(bpm.get_adjacency().at(0).count({0, 1, 2}), 0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 2}), 22.0);

}

TEST(FunctionTestBPM, remove_variable_String) {
   
   Polynomial<std::string, double> polynomial {
      {{"a"}, 0.0}, {{"b"}, 1.0}, {{"c"}, 2.0},
      {{"a", "b"}, 11.0}, {{"a", "c"}, 22.0}, {{"b", "c"}, 12.0},
      {{"a", "b", "c"}, 123.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<std::string, double> bpm(polynomial, vartype);
   
   bpm.remove_variable("b");
   
   EXPECT_TRUE(bpm.get_polynomial().find({"a"}) != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"c"}) != bpm.get_polynomial().end());
   
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "c"}) != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b", "c"}) == bpm.get_polynomial().end());
   
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b", "c"}) == bpm.get_polynomial().end());
   
   EXPECT_EQ(bpm.length(), 2);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count("a"), 1);
   EXPECT_EQ(bpm.get_variables().count("b"), 0);
   EXPECT_EQ(bpm.get_variables().count("c"), 1);
   
   //Adjacency
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b"}), 0);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "c"}), 1);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b", "c"}), 0);
   EXPECT_EQ(bpm.get_adjacency().count("b"), 0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "c"}), 22.0);
   
}

TEST(FunctionTestBPM, remove_variable_from) {
   
   Polynomial<std::string, double> polynomial {
      //linear biases
      {{"a"}, 1.0}, {{"b"}, 2.0}, {{"c"}, 3.0}, {{"d"}, 4.0},
      //quadratic biases
      {{"a", "b"}, 12.0}, {{"a", "c"}, 13.0}, {{"a", "d"}, 14.0},
      {{"b", "c"}, 23.0}, {{"b", "d"}, 24.0},
      {{"c", "d"}, 34.0},
      //polynomial biases
      {{"a", "b", "c"}, 123.0}, {{"a", "b", "d"}, 124.0}, {{"a", "c", "d"}, 134.0},
      {{"b", "c", "d"}, 234.0},
      {{"a", "b", "c", "d"}, 1234.0}
   };
   
   Vartype vartype = Vartype::BINARY;
   BinaryPolynomialModel<std::string, double> bpm(polynomial, vartype);
   
   std::vector<std::string> variables = {"b", "d"};
   bpm.remove_variables_from(variables);
   
   EXPECT_TRUE(bpm.get_polynomial().find({"a"}) != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"c"}) != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"d"}) == bpm.get_polynomial().end());

   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "c"}) != bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "d"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b", "c"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b", "d"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"c", "d"}) == bpm.get_polynomial().end());

   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b", "c"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b", "d"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"a", "c", "d"}) == bpm.get_polynomial().end());
   EXPECT_TRUE(bpm.get_polynomial().find({"b", "c", "d"}) == bpm.get_polynomial().end());

   EXPECT_TRUE(bpm.get_polynomial().find({"a", "b", "c", "d"}) == bpm.get_polynomial().end());
   
   EXPECT_EQ(bpm.length(), 2);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count("a"), 1);
   EXPECT_EQ(bpm.get_variables().count("b"), 0);
   EXPECT_EQ(bpm.get_variables().count("c"), 1);
   EXPECT_EQ(bpm.get_variables().count("d"), 0);
   
   //Adjacency
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b"}), 0);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "c"}), 1);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "d"}), 0);
   EXPECT_EQ(bpm.get_adjacency().count("b"), 0);
   EXPECT_EQ(bpm.get_adjacency().at("c").count({"c", "d"}), 0);

   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b", "c"}), 0);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b", "d"}), 0);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "c", "d"}), 0);
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"b", "c", "d"}), 0);
   
   EXPECT_EQ(bpm.get_adjacency().at("a").count({"a", "b", "c", "d"}), 0);

}

TEST(FunctionTestBPM, scale) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.scale(0.5);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 0.5);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 1.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 5.5);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 2}), 11.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1, 2}), 6.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 6.0);
   
}

TEST(FunctionTestBPM, normalize) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.normalize({-1, 1});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0}), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 1.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 2.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 11.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 2}), 22.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1, 2}), 12.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 12.0/22.0);
   
}

TEST(FunctionTestBPM, _generate_indices) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 12.0},
      {{0, 1, 4}, 13.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   std::vector<uint32_t> variables = bpm._generate_indices();
   EXPECT_EQ(variables.at(0), 0);
   EXPECT_EQ(variables.at(1), 1);
   EXPECT_EQ(variables.at(2), 2);
   EXPECT_EQ(variables.at(3), 4);
   
   EXPECT_TRUE(bpm.contains(0));
   EXPECT_TRUE(bpm.contains(1));
   EXPECT_TRUE(bpm.contains(2));
   EXPECT_FALSE(bpm.contains(3));
   EXPECT_TRUE(bpm.contains(4));


}

TEST(FunctionTestBPM, from_serializable) {
   
   Polynomial<std::string, double> polynomial {
      //linear biases
      {{"a"}, 1.0}, {{"b"}, 2.0}, {{"c"}, 3.0}, {{"d"}, 4.0},
      //quadratic biases
      {{"a", "b"}, 12.0}, {{"a", "c"}, 13.0}, {{"a", "d"}, 14.0},
      {{"b", "c"}, 23.0}, {{"b", "d"}, 24.0},
      {{"c", "d"}, 34.0},
      //polynomial biases
      {{"a", "b", "c"}, 123.0}, {{"a", "b", "d"}, 124.0}, {{"a", "c", "d"}, 134.0},
      {{"b", "c", "d"}, 234.0},
      {{"a", "b", "c", "d"}, 1234.0}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<std::string, double> bpm(polynomial, vartype);
   json j = bpm.to_serializable();

   BinaryPolynomialModel<std::string, double> bpm_from = BinaryPolynomialModel<std::string, double>::from_serializable(j);
   
   EXPECT_EQ(bpm.length(), bpm_from.length());
   
   //variables
   EXPECT_EQ(bpm.get_variables().count("a"), bpm_from.get_variables().count("a"));
   EXPECT_EQ(bpm.get_variables().count("b"), bpm_from.get_variables().count("b"));
   EXPECT_EQ(bpm.get_variables().count("c"), bpm_from.get_variables().count("c"));
   EXPECT_EQ(bpm.get_variables().count("d"), bpm_from.get_variables().count("d"));
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b"}), bpm_from.get_adjacency().at("a").at({"a", "b"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "c"}), bpm_from.get_adjacency().at("a").at({"a", "c"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "d"}), bpm_from.get_adjacency().at("a").at({"a", "d"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "c"}), bpm_from.get_adjacency().at("a").at({"a", "b", "c"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "d"}), bpm_from.get_adjacency().at("a").at({"a", "b", "d"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "c", "d"}), bpm_from.get_adjacency().at("a").at({"a", "c", "d"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("a").at({"a", "b", "c", "d"}), bpm_from.get_adjacency().at("a").at({"a", "b", "c", "d"}));
   
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "c"}), bpm_from.get_adjacency().at("b").at({"b", "c"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "d"}), bpm_from.get_adjacency().at("b").at({"b", "d"}));
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("b").at({"b", "c", "d"}), bpm_from.get_adjacency().at("b").at({"b", "c", "d"}));

   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at("c").at({"c", "d"}), bpm_from.get_adjacency().at("c").at({"c", "d"}));
   
   for (const auto &it: polynomial) {
      EXPECT_DOUBLE_EQ(bpm.get_polynomial().at(it.first), bpm_from.get_polynomial().at(it.first));
   }
   
   EXPECT_EQ(bpm.get_vartype(), bpm_from.get_vartype());
}

TEST(FunctionTestBPM, from_pubo) {
   Polynomial<uint32_t, double> polynomial {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_pubo(polynomial);

   EXPECT_EQ(bpm.length(), 4);

   //variables
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   EXPECT_EQ(bpm.get_variables().count(4), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2}), 12.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 3}), 13.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 4}), 14.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 3}), 123.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 4}), 124.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 3, 4}), 134.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 3, 4}), 1234.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 3}), 23.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 4}), 24.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 3, 4}), 234.0);

   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(3).at({3, 4}), 34.0);
   
   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(bpm.get_vartype(), Vartype::BINARY);
}

TEST(FunctionTestBPM, from_ising) {
   Polynomial<uint32_t, double> polynomial {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_ising(polynomial);

   EXPECT_EQ(bpm.length(), 4);

   //variables
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   EXPECT_EQ(bpm.get_variables().count(4), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2}), 12.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 3}), 13.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 4}), 14.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 3}), 123.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 4}), 124.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 3, 4}), 134.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2, 3, 4}), 1234.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 3}), 23.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 4}), 24.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(2).at({2, 3, 4}), 234.0);

   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(3).at({3, 4}), 34.0);
   
   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(bpm.get_vartype(), Vartype::SPIN);
}

TEST(FunctionTestBPM, empty) {
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};

   Vartype vartype = Vartype::BINARY;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);

   bpm.empty();
   
   EXPECT_TRUE(bpm.get_polynomial().empty());
   EXPECT_TRUE(bpm.get_variables().empty());
   EXPECT_TRUE(bpm.get_adjacency().empty());
   EXPECT_EQ(bpm.get_vartype(), Vartype::NONE);
   EXPECT_EQ(bpm.get_info(), "");
   
   //Chech if the methods in Binary Polynomial Model work properly after executing empty()
   EXPECT_EQ(bpm.length(), 0);
   bpm.remove_variable(1);
   bpm.remove_variables_from(std::vector<uint32_t>{1,2,3,4,5});
   bpm.remove_interaction(std::vector<uint32_t>{1,2});
   bpm.remove_interactions_from(std::vector<std::vector<uint32_t>>{{1,2},{1,3},{1,4}});
   bpm.scale(1.0);
   bpm.normalize();
   
   //energy
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_1), 0.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_2), 0.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_3), 0.0);
   
   //Reset polynomial model
   bpm.add_interaction({0, 1}   , 11.0, Vartype::BINARY);
   bpm.add_interaction({0, 2}   , 22.0);
   bpm.add_interaction({1, 2}   , 12.0);
   bpm.add_interaction({0, 1, 2}, 123.0);
   
   bpm.add_linear(0, 0.0);
   bpm.add_linear(1, 1.0);
   bpm.add_linear(2, 2.0);
   
   //Check energy
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_3), 0.0   );
   
   EXPECT_EQ(bpm.length(), 3);
   
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   
   //Adjacency
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1}), 11.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 2}), 22.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(0).at({0, 1, 2}), 123.0);
   EXPECT_DOUBLE_EQ(bpm.get_adjacency().at(1).at({1, 2}), 12.0);

   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(bpm.get_vartype(), Vartype::BINARY);
   
}

*/

}
