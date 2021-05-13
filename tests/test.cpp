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

//google test for binary polynomial model
bool EXPECT_CONTAIN(double val, const PolynomialValueList<double> &poly_value) {
   int count = 0;
   for (const auto &it: poly_value) {
      if (std::abs(it - val) < 0.000000000000001/*std::pow(10, -15)*/) {
         count++;
      }
   }
   if (count >= 1) {
      return true;
   }
   else {
      return false;
   }
}

void StateTestBPMUINT(const BinaryPolynomialModel<uint32_t, double> &bpm) {
   
   EXPECT_EQ(bpm.get_num_variables(), 4);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 0.0);
   
   EXPECT_EQ(bpm.get_num_interactions(), 15);
   
   EXPECT_EQ(bpm.get_degree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1}         ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2}         ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({3}         ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({4}         ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 2}      ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 3}      ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 4}      ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 3}      ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 4}      ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({3, 4}      ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 2, 3}   ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 2, 4}   ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 3, 4}   ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 3, 4}   ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 2, 3, 4}), 1234.0);
   
   //variables_to_integers
   EXPECT_EQ(bpm.get_variables_to_integers(1), 0);
   EXPECT_EQ(bpm.get_variables_to_integers(2), 1);
   EXPECT_EQ(bpm.get_variables_to_integers(3), 2);
   EXPECT_EQ(bpm.get_variables_to_integers(4), 3);

   //Polynomial Key
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1}         ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2}         ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{3}         ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{4}         ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 2}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 3}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 4}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 3}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 4}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{3, 4}      ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 2, 3}   ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 2, 4}   ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{1, 2, 3, 4}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm._get_values()));

   //sorted_variables
   auto sorted_variables = bpm.get_sorted_variables();
   EXPECT_EQ(sorted_variables[0], 1);
   EXPECT_EQ(sorted_variables[1], 2);
   EXPECT_EQ(sorted_variables[2], 3);
   EXPECT_EQ(sorted_variables[3], 4);
   
   //variables
   EXPECT_EQ(bpm.GetVariables().count(1), 1);
   EXPECT_EQ(bpm.GetVariables().count(2), 1);
   EXPECT_EQ(bpm.GetVariables().count(3), 1);
   EXPECT_EQ(bpm.GetVariables().count(4), 1);
}

void StateTestBPMINT(const BinaryPolynomialModel<int32_t, double> &bpm) {
   
   EXPECT_EQ(bpm.get_num_variables(), 4);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 0.0);
   
   EXPECT_EQ(bpm.get_num_interactions(), 15);
   
   EXPECT_EQ(bpm.get_degree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-1}            ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-2}            ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-3}            ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4}            ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-2, -1}        ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-3, -1}        ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -1}        ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-3, -2}        ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -2}        ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -3}        ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-3, -2, -1}    ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -2, -1}    ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -3, -1}    ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -3, -2}    ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({-4, -3, -2, -1}), 1234.0);
   
   //variables_to_integers
   EXPECT_EQ(bpm.get_variables_to_integers(-4), 0);
   EXPECT_EQ(bpm.get_variables_to_integers(-3), 1);
   EXPECT_EQ(bpm.get_variables_to_integers(-2), 2);
   EXPECT_EQ(bpm.get_variables_to_integers(-1), 3);
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-1}            ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-2}            ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-3}            ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4}            ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-2, -1}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-3, -1}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -1}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-3, -2}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -2}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -3}        ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-3, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -3, -1}    ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -3, -2}    ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<int32_t>{-4, -3, -2, -1}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm._get_values()));

   //sorted_variables
   auto sorted_variables = bpm.get_sorted_variables();
   EXPECT_EQ(sorted_variables[0], -4);
   EXPECT_EQ(sorted_variables[1], -3);
   EXPECT_EQ(sorted_variables[2], -2);
   EXPECT_EQ(sorted_variables[3], -1);
   
   //variables
   EXPECT_EQ(bpm.GetVariables().count(-1), 1);
   EXPECT_EQ(bpm.GetVariables().count(-2), 1);
   EXPECT_EQ(bpm.GetVariables().count(-3), 1);
   EXPECT_EQ(bpm.GetVariables().count(-4), 1);
}

void StateTestBPMString(const BinaryPolynomialModel<std::string, double> &bpm) {
   
   EXPECT_EQ(bpm.get_num_variables(), 4);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 0.0);

   EXPECT_EQ(bpm.get_num_interactions(), 15);
   
   EXPECT_EQ(bpm.get_degree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a"}               ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"b"}               ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"c"}               ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"d"}               ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "b"}          ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "c"}          ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "d"}          ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"b", "c"}          ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"b", "d"}          ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"c", "d"}          ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "b", "c"}     ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "b", "d"}     ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "c", "d"}     ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"b", "c", "d"}     ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({"a", "b", "c", "d"}), 1234.0);
   
   //variables_to_integers
   EXPECT_EQ(bpm.get_variables_to_integers("a"), 0);
   EXPECT_EQ(bpm.get_variables_to_integers("b"), 1);
   EXPECT_EQ(bpm.get_variables_to_integers("c"), 2);
   EXPECT_EQ(bpm.get_variables_to_integers("d"), 3);
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a"}               ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"b"}               ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"c"}               ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"d"}               ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "b"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"b", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"b", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"c", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "b", "c"}     ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "b", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"b", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<std::string>{"a", "b", "c", "d"}), 1);

   //Polynomial Value
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm._get_values()));

   //sorted_variables
   auto sorted_variables = bpm.get_sorted_variables();
   EXPECT_EQ(sorted_variables[0], "a");
   EXPECT_EQ(sorted_variables[1], "b");
   EXPECT_EQ(sorted_variables[2], "c");
   EXPECT_EQ(sorted_variables[3], "d");
   
   //variables
   EXPECT_EQ(bpm.GetVariables().count("a"), 1);
   EXPECT_EQ(bpm.GetVariables().count("b"), 1);
   EXPECT_EQ(bpm.GetVariables().count("c"), 1);
   EXPECT_EQ(bpm.GetVariables().count("d"), 1);
}

template<typename IndexType>
void StateTestBPMEmpty(const BinaryPolynomialModel<IndexType, double> &bpm) {
   
   EXPECT_EQ(bpm.get_num_variables(), 0);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 0.0);
   
   EXPECT_EQ(bpm.get_num_interactions(), 0);
   
   EXPECT_EQ(bpm.get_degree(), 0);
      
   //Polynomial
   EXPECT_EQ(bpm.get_polynomial().size(), 0);
   
   //variables_to_integers
   EXPECT_EQ(bpm.get_variables_to_integers().size(), 0);
   
   //Polynomial Key
   EXPECT_EQ(bpm._get_keys().size(), 0);
   
   //Polynomial Val
   EXPECT_EQ(bpm._get_values().size(), 0);

   //sorted_variables
   auto sorted_variables = bpm.get_sorted_variables();
   EXPECT_EQ(sorted_variables.size(), 0);
   
   //variables
   EXPECT_EQ(bpm.GetVariables().size(), 0);

}

Polynomial<uint32_t, double> GeneratePolynomialUINT() {
   
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
   
   return polynomial;
}

Polynomial<int32_t, double> GeneratePolynomialINT() {
   
   Polynomial<int32_t, double> polynomial {
      //linear biases
      {{-1}, 1.0}, {{-2}, 2.0}, {{-3}, 3.0}, {{-4}, 4.0},
      //quadratic biases
      {{-1, -2}, 12.0}, {{-1, -3}, 13.0}, {{-1, -4}, 14.0},
      {{-2, -3}, 23.0}, {{-2, -4}, 24.0},
      {{-3, -4}, 34.0},
      //polynomial biases
      {{-1, -2, -3}, 123.0}, {{-1, -2, -4}, 124.0}, {{-1, -3, -4}, 134.0},
      {{-2, -3, -4}, 234.0},
      {{-1, -2, -3, -4}, 1234.0}
   };
   
   return polynomial;
}

Polynomial<std::string, double> GeneratePolynomialString() {
   
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
   
   return polynomial;
}

TEST(ConstructionBPM, PolyMapUINT) {
   BinaryPolynomialModel<uint32_t, double> bpm(GeneratePolynomialUINT(), Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMUINT(bpm);
}

TEST(ConstructionBPM, PolyMapINT) {
   BinaryPolynomialModel<int32_t, double> bpm(GeneratePolynomialINT(), Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMINT(bpm);
}

TEST(ConstructionBPM, PolyMapString) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMString(bpm);
}

TEST(ConstructionBPM, PolyKeyValueUINT) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   BinaryPolynomialModel<uint32_t, double> bpm(poly_key, poly_value, Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMUINT(bpm);
}

TEST(ConstructionBPM, PolyKeyValueINT) {
   
   PolynomialKeyList<int32_t>  poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   BinaryPolynomialModel<int32_t, double> bpm(poly_key, poly_value, Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMINT(bpm);
}

TEST(ConstructionBPM, PolyKeyValueStrign) {
   
   PolynomialKeyList<std::string> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialString()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   BinaryPolynomialModel<std::string, double> bpm(poly_key, poly_value, Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMString(bpm);
}

TEST(add_interactionBPM, basic) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.add_interaction({3}, 3.0);
   bpm.add_interaction({4}, 4.0);
   bpm.add_interaction({1, 3}, 13.0);
   bpm.add_interaction({1, 4}, 14.0);
   bpm.add_interaction({2, 3}, 23.0);
   bpm.add_interaction({2, 4}, 24.0);
   bpm.add_interaction({3, 4}, 34.0);
   bpm.add_interaction({1, 2, 4}, 124.0);
   bpm.add_interaction({1, 3, 4}, 134.0);
   bpm.add_interaction({2, 3, 4}, 234.0);
   bpm.add_interaction({1, 2, 3, 4}, 1234.0);
   
   StateTestBPMUINT(bpm);

}

TEST(add_interactionBPM, duplicate_value_1) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second /= 2;
   }
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.add_interaction(it.first, it.second);
   };
   
   StateTestBPMUINT(bpm);

}

TEST(add_interactionBPM, duplicate_value_2) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.add_interaction(it.first, -it.second);
   };
   
   StateTestBPMEmpty(bpm);
}

TEST(add_interactions_fromBPM, PolyMap) {
      
   BinaryPolynomialModel<uint32_t, double> bpm({}, Vartype::SPIN);
   
   bpm.add_interactions_from(GeneratePolynomialUINT());
   
   StateTestBPMUINT(bpm);
   
}

TEST(add_interactions_fromBPM, PolyKeyValue) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
   
   BinaryPolynomialModel<uint32_t, double> bpm({}, Vartype::SPIN);
   
   bpm.add_interactions_from(poly_key, poly_value);
   
   StateTestBPMUINT(bpm);
   
}



TEST(add_offsetBPM, basic) {
   Polynomial<uint32_t, double> polynomial;
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.add_offset(3.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_offset(), 3.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}) , 3.0);
   
   bpm.add_offset(3.0);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 6.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}) , 6.0);
   
   bpm.add_offset(-6.0);
   
   StateTestBPMEmpty(bpm);
   
}

TEST(remove_interactionBPM, basic) {
   
   Polynomial<uint32_t, double> polynomial = {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //To be removed
      {{11, 12, 14}, -1.0}, {{7}, -2.0}, {{2, 11}, -9.0}, {{}, -3},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}          ), -3.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({7}         ), -2.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 11}     ), -9.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({11, 12, 14}), -1.0);
   
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{}          ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{7}         ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 11}     ), 1);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{11, 12, 14}), 1);
   
   EXPECT_TRUE(EXPECT_CONTAIN(-3.0, bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(-2.0, bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(-9.0, bpm._get_values()));
   EXPECT_TRUE(EXPECT_CONTAIN(-1.0, bpm._get_values()));

   bpm.remove_interaction({}          );
   bpm.remove_interaction({7}         );
   bpm.remove_interaction({2, 11}     );
   bpm.remove_interaction({11, 12, 14});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}          ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({7}         ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 11}     ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({11, 12, 14}), 0.0);
   
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm._get_values()));
   
   StateTestBPMUINT(bpm);
   
}

TEST(remove_interactionBPM, remove_all) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.remove_interaction(it.first);
   };
   
   StateTestBPMEmpty(bpm);
   
}

TEST(remove_interactions_fromBPM, basic) {
   
   Polynomial<uint32_t, double> polynomial {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //To be removed
      {{11, 12, 14}, -1.0}, {{7}, -2.0}, {{2, 11}, -9.0}, {{}, -3},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   PolynomialKeyList<uint32_t> removed_key_list = {
      {11, 14, 12}, {7}, {11, 2}, {}
   };
   
   bpm.remove_interactions_from(removed_key_list);
      
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}          ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({7}         ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 11}     ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({11, 12, 14}), 0.0);
   
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm._get_keys().begin(), bpm._get_keys().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm._get_values()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm._get_values()));
   
   StateTestBPMUINT(bpm);
   
}

TEST(remove_interactions_fromBPM, remove_all) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   PolynomialKeyList<uint32_t> removed_key_list;
   
   for (const auto &it: polynomial) {
      removed_key_list.push_back(it.first);
   }
   
   bpm.remove_interactions_from(removed_key_list);
   
   StateTestBPMEmpty(bpm);
   
}

TEST(remove_offsetBPM, basic) {
   Polynomial<uint32_t, double> polynomial = {
      //linear biases
      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0}, {{4}, 4.0},
      //quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      //To be removed
      {{}, 100},
      //polynomial biases
      {{1, 2, 3}, 123.0}, {{1, 2, 4}, 124.0}, {{1, 3, 4}, 134.0},
      {{2, 3, 4}, 234.0},
      {{1, 2, 3, 4}, 1234.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.remove_offset();
   
   StateTestBPMUINT(bpm);
}

TEST(energyBPM, SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   Sample<uint32_t> sample_variables_spin_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_2{{0, +1}, {1, -1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_3{{0, -1}, {1, -1}, {2, -1}};

   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_2), -123.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_spin_3), -81.0 );
   
   std::vector<int32_t> sample_vec_variables_spin_1{+1, +1, +1};
   std::vector<int32_t> sample_vec_variables_spin_2{+1, -1, +1};
   std::vector<int32_t> sample_vec_variables_spin_3{-1, -1, -1};
   
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_spin_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_spin_2), -123.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_spin_3), -81.0 );

}

TEST(energyBPM, BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::BINARY);
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};

   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_3), 0.0   );
   
   std::vector<int32_t> sample_vec_variables_binary_1{+1, +1, +1};
   std::vector<int32_t> sample_vec_variables_binary_2{+1, +0, +1};
   std::vector<int32_t> sample_vec_variables_binary_3{+0, +0, +0};
   
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.energy(sample_vec_variables_binary_3), 0.0   );
      
}

TEST(energiesBPM, SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   std::vector<Sample<uint32_t>> sample_variables_spin {
      {{0, +1}, {1, +1}, {2, +1}},
      {{0, +1}, {1, -1}, {2, +1}},
      {{0, -1}, {1, -1}, {2, -1}}
   };
   std::vector<double> en_vec = bpm.energies(sample_variables_spin);
   
   EXPECT_DOUBLE_EQ(en_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec[1], -123.0);
   EXPECT_DOUBLE_EQ(en_vec[2], -81.0 );
   
   std::vector<std::vector<int32_t>> sample_vec_variables_spin {
      {+1, +1, +1},
      {+1, -1, +1},
      {-1, -1, -1}
   };
   
   std::vector<double> en_vec_vec = bpm.energies(sample_vec_variables_spin);
   
   EXPECT_DOUBLE_EQ(en_vec_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[1], -123.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[2], -81.0 );
   
}

TEST(energiesBPM, BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   std::vector<Sample<uint32_t>> sample_variables_binary {
      {{0, +1}, {1, +1}, {2, +1}},
      {{0, +1}, {1, +0}, {2, +1}},
      {{0, +0}, {1, +0}, {2, +0}}
   };
   
   std::vector<double> en_vec = bpm.energies(sample_variables_binary);
   
   EXPECT_DOUBLE_EQ(en_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec[1], +24.0 );
   EXPECT_DOUBLE_EQ(en_vec[2], 0.0   );
   
   std::vector<std::vector<int32_t>> sample_vec_variables_binary {
      {+1, +1, +1},
      {+1, +0, +1},
      {+0, +0, +0}
   };
   
   std::vector<double> en_vec_vec = bpm.energies(sample_vec_variables_binary);
   
   EXPECT_DOUBLE_EQ(en_vec_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[1], +24.0 );
   EXPECT_DOUBLE_EQ(en_vec_vec[2], 0.0   );
   
}

TEST(scaleBPM, all_scale) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second *= 2;
   }
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.scale(0.5);
   
   StateTestBPMUINT(bpm);
   
}

TEST(scaleBPM, ignored_interaction) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second *= 2;
   }
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.scale(0.5, {{1,2}, {2, 4}, {1, 3, 4}});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 2}   ), 12.0*2 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 4}   ), 24.0*2 );
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({1, 3, 4}), 134.0*2);

   bpm.add_interaction({1, 2}   , -12.0);
   bpm.add_interaction({2, 4}   , -24.0);
   bpm.add_interaction({1, 3, 4}, -134.0);
   
   StateTestBPMUINT(bpm);

}

TEST(scaleBPM, ignored_offset) {
   
   Polynomial<uint32_t, double> polynomial {
      {{}, 100.0}
   };
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.scale(0.5, {std::vector<uint32_t>{}});
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}), 100.0);
   
   bpm.scale(0.5, {}, true);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}), 100.0);
   
   bpm.scale(0.5, {{}}, true);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}), 100.0);
   
   bpm.scale(0.5);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}), 50.0);
   
}

TEST(normalizeBPM, all_normalize) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.normalize({-1, 1});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 1.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 2.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 11.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 2}), 22.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1, 2}), 12.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 12.0/22.0);
   
}

TEST(normalizeBPM, ignored_interaction) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.normalize({-1, 1}, {{0, 1, 2}});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1}), 1.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({2}), 2.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1}), 11.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 2}), 22.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({1, 2}), 12.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial().at({0, 1, 2}), 12.0);
   
}

TEST(serializableBPM, UINT) {
   BinaryPolynomialModel<uint32_t, double> bpm(GeneratePolynomialUINT(), Vartype::SPIN);
   BinaryPolynomialModel<uint32_t, double> bpm_from = BinaryPolynomialModel<uint32_t, double>::from_serializable(bpm.to_serializable());
   StateTestBPMUINT(bpm_from);
}

TEST(serializableBPM, INT) {
   BinaryPolynomialModel<int32_t, double> bpm(GeneratePolynomialINT(), Vartype::SPIN);
   BinaryPolynomialModel<int32_t, double> bpm_from = BinaryPolynomialModel<int32_t, double>::from_serializable(bpm.to_serializable());
   StateTestBPMINT(bpm_from);
}

TEST(serializableBPM, String) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   BinaryPolynomialModel<std::string, double> bpm_from = BinaryPolynomialModel<std::string, double>::from_serializable(bpm.to_serializable());
   StateTestBPMString(bpm_from);
}

TEST(serializableBPM, StringToUINT) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   auto obj = bpm.to_serializable();
   std::vector<std::size_t> string_to_num(obj["variables"].size());
   std::iota(string_to_num.begin(), string_to_num.end(), 1);
   obj["variables"] = string_to_num;
   BinaryPolynomialModel<uint32_t, double> bpm_from = BinaryPolynomialModel<uint32_t, double>::from_serializable(obj);
   StateTestBPMUINT(bpm_from);
}

TEST(serializableBPM, StringToINT) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   auto obj = bpm.to_serializable();
   std::vector<std::size_t> string_to_num(4);
   string_to_num[0] = -1;
   string_to_num[1] = -2;
   string_to_num[2] = -3;
   string_to_num[3] = -4;
   obj["variables"] = string_to_num;
   BinaryPolynomialModel<int32_t, double> bpm_from = BinaryPolynomialModel<int32_t, double>::from_serializable(obj);
   StateTestBPMINT(bpm_from);
}

TEST(from_hubo, MapUINT) {
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hubo(GeneratePolynomialUINT());
   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   StateTestBPMUINT(bpm);
}

TEST(from_hubo, MapINT) {
   auto bpm = BinaryPolynomialModel<int32_t, double>::from_hubo(GeneratePolynomialINT());
   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   StateTestBPMINT(bpm);
}

TEST(from_hubo, MapString) {
   auto bpm = BinaryPolynomialModel<std::string, double>::from_hubo(GeneratePolynomialString());
   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   StateTestBPMString(bpm);
}

TEST(from_hubo, KeyValueUINT) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   
   StateTestBPMUINT(bpm);
}

TEST(from_hubo, KeyValueINT) {
   PolynomialKeyList<int32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<int32_t, double>::from_hubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   
   StateTestBPMINT(bpm);
}

TEST(from_hubo, KeyValueString) {
   PolynomialKeyList<std::string> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialString()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<std::string, double>::from_hubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.get_vartype());
   
   StateTestBPMString(bpm);
}

TEST(from_hising, MapUINT) {
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hising(GeneratePolynomialUINT());
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMUINT(bpm);
}

TEST(from_hising, MapINT) {
   auto bpm = BinaryPolynomialModel<int32_t, double>::from_hising(GeneratePolynomialINT());
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMINT(bpm);
}

TEST(from_hising, MapString) {
   auto bpm = BinaryPolynomialModel<std::string, double>::from_hising(GeneratePolynomialString());
   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   StateTestBPMString(bpm);
}

TEST(from_hising, KeyValueUINT) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   
   StateTestBPMUINT(bpm);
}

TEST(from_hising, KeyValueINT) {
   PolynomialKeyList<int32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<int32_t, double>::from_hising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   
   StateTestBPMINT(bpm);
}

TEST(from_hising, KeyValueString) {
   PolynomialKeyList<std::string> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialString()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<std::string, double>::from_hising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.get_vartype());
   
   StateTestBPMString(bpm);
}

TEST(clearBPM, basic) {
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };
   
   Vartype vartype = Vartype::BINARY;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);

   bpm.clear();
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};
   
   EXPECT_TRUE(bpm.get_polynomial().empty());
   EXPECT_TRUE(bpm.GetVariables().empty());
   EXPECT_EQ(bpm.get_vartype(), Vartype::BINARY);
   
   //Chech if the methods in Binary Polynomial Model work properly after executing empty()
   EXPECT_EQ(bpm.get_num_variables(), 0);
   bpm.remove_variable(1);
   bpm.remove_variables_from(std::vector<uint32_t>{1,2,3,4,5});
   bpm.remove_interaction(std::vector<uint32_t>{1,2});
   bpm.remove_interactions_from(std::vector<std::vector<uint32_t>>{{1,2},{1,3},{1,4}});
   bpm.scale(1.0);
   bpm.normalize();
   
   //energy
   EXPECT_THROW(bpm.energy(sample_variables_binary_1), std::runtime_error);
   EXPECT_THROW(bpm.energy(sample_variables_binary_2), std::runtime_error);
   EXPECT_THROW(bpm.energy(sample_variables_binary_3), std::runtime_error);
   
   //Reset polynomial model
   bpm.add_interaction({0, 1}   , 11.0, Vartype::BINARY);
   bpm.add_interaction({0, 2}   , 22.0);
   bpm.add_interaction({1, 2}   , 12.0);
   bpm.add_interaction({0, 1, 2}, 123.0);
   
   bpm.add_interaction({0}, 0.0);
   bpm.add_interaction({1}, 1.0);
   bpm.add_interaction({2}, 2.0);
   
   //Check energy
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.energy(sample_variables_binary_3), 0.0   );
   
   EXPECT_EQ(bpm.get_num_variables(), 3);
   
   EXPECT_EQ(bpm.GetVariables().count(0), 1);
   EXPECT_EQ(bpm.GetVariables().count(1), 1);
   EXPECT_EQ(bpm.GetVariables().count(2), 1);
   
   for (const auto &it: bpm.get_polynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(bpm.get_vartype(), Vartype::BINARY);
   
}

TEST(vartypeBPM, SpinBinarySpin) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hising(GeneratePolynomialUINT());
   auto bpm_binary = BinaryPolynomialModel<uint32_t, double>(bpm.to_hubo(), Vartype::BINARY);
   auto bpm_ising  = BinaryPolynomialModel<uint32_t, double>(bpm_binary.to_hising(), Vartype::SPIN);
   
   StateTestBPMUINT(bpm_ising);
   
}

TEST(vartypeBPM, BinarySPINBinary) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hubo(GeneratePolynomialUINT());
   auto bpm_ising  = BinaryPolynomialModel<uint32_t, double>(bpm.to_hising(), Vartype::SPIN);
   auto bpm_bianry = BinaryPolynomialModel<uint32_t, double>(bpm_ising.to_hubo(), Vartype::BINARY);
   
   StateTestBPMUINT(bpm_bianry);
   
}

TEST(vartypeBPM, change_vartypeSpinBinarySpin) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hising(GeneratePolynomialUINT());
   auto bpm_binary = bpm.change_vartype(Vartype::BINARY, true);
   auto bpm_ising  = bpm_binary.change_vartype(Vartype::SPIN, true);
   
   StateTestBPMUINT(bpm_ising);
   bpm.change_vartype(Vartype::SPIN);
   StateTestBPMUINT(bpm);
   
}

TEST(vartypeBPM, change_vartypeBinarySPINBinary) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::from_hubo(GeneratePolynomialUINT());
   auto bpm_ising  = bpm.change_vartype(Vartype::SPIN, true);
   auto bpm_binary = bpm_ising.change_vartype(Vartype::BINARY, true);
   
   StateTestBPMUINT(bpm_binary);
   bpm.change_vartype(Vartype::BINARY);
   StateTestBPMUINT(bpm);
   
}



}

