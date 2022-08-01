//    Copyright 2022 Jij Inc.

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <tuple>

#include <cimod/binary_quadratic_model.hpp>
#include <cimod/binary_polynomial_model.hpp>
#include <cimod/binary_quadratic_model_dict.hpp>

#include "test_bqm.hpp"

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
   
   EXPECT_EQ(bpm.GetNumVariables(), 4);

   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 0.0);
   
   EXPECT_EQ(bpm.GetNumInteractions(), 15);
   
   EXPECT_EQ(bpm.GetDegree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1}         ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2}         ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({3}         ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({4}         ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 2}      ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 3}      ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 4}      ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 3}      ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 4}      ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({3, 4}      ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 2, 3}   ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 2, 4}   ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 3, 4}   ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 3, 4}   ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 2, 3, 4}), 1234.0);
   
   //Polynomial duplicate key
   if (bpm.GetVartype() == cimod::Vartype::SPIN) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 1, 1}         ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 1, 1, 2}      ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 3, 3, 3}      ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({3, 2, 3, 2, 3, 2}), 23.0);
   }
   else if (bpm.GetVartype() == cimod::Vartype::BINARY) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 1, 1, 1}         ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 1, 1, 2, 2}      ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 3, 3, 3, 1}      ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({3, 2, 3, 2, 3, 2, 2}), 23.0);
   }
   
   //variables_to_integers
   EXPECT_EQ(bpm.GetVariablesToIntegers(1), 0);
   EXPECT_EQ(bpm.GetVariablesToIntegers(2), 1);
   EXPECT_EQ(bpm.GetVariablesToIntegers(3), 2);
   EXPECT_EQ(bpm.GetVariablesToIntegers(4), 3);

   //Polynomial Key
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1}         ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2}         ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{3}         ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{4}         ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 2}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 3}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 3}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{3, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 2, 3}   ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 2, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{1, 2, 3, 4}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.GetValueList()));

   //sorted_variables
   auto sorted_variables = bpm.GetSortedVariables();
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
   
   EXPECT_EQ(bpm.GetNumVariables(), 4);

   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 0.0);
   
   EXPECT_EQ(bpm.GetNumInteractions(), 15);
   
   EXPECT_EQ(bpm.GetDegree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1}            ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-2}            ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3}            ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4}            ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-2, -1}        ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3, -1}        ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -1}        ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3, -2}        ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -2}        ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -3}        ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3, -2, -1}    ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -2, -1}    ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -3, -1}    ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -3, -2}    ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-4, -3, -2, -1}), 1234.0);
   
   //Polynomial duplicate key
   if (bpm.GetVartype() == cimod::Vartype::SPIN) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -1, -1}            ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -1, -1, -2}        ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -3, -3, -3}        ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3, -2, -3, -2, -3, -2}), 23.0);
   }
   else if (bpm.GetVartype() == cimod::Vartype::BINARY) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -1, -1, -1}            ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -1, -1, -2, -2}        ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-1, -3, -3, -3, -1}        ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({-3, -2, -3, -2, -3, -2, -2}), 23.0);
   }
   
   //variables_to_integers
   EXPECT_EQ(bpm.GetVariablesToIntegers(-4), 0);
   EXPECT_EQ(bpm.GetVariablesToIntegers(-3), 1);
   EXPECT_EQ(bpm.GetVariablesToIntegers(-2), 2);
   EXPECT_EQ(bpm.GetVariablesToIntegers(-1), 3);
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-1}            ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-2}            ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-3}            ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4}            ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-2, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-3, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-3, -2}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -2}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -3}        ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-3, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -3, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -3, -2}    ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<int32_t>{-4, -3, -2, -1}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.GetValueList()));

   //sorted_variables
   auto sorted_variables = bpm.GetSortedVariables();
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
   
   EXPECT_EQ(bpm.GetNumVariables(), 4);

   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 0.0);

   EXPECT_EQ(bpm.GetNumInteractions(), 15);
   
   EXPECT_EQ(bpm.GetDegree(), 4);
      
   //Polynomial
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a"}               ), 1.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"b"}               ), 2.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"c"}               ), 3.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"d"}               ), 4.0   );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "b"}          ), 12.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "c"}          ), 13.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "d"}          ), 14.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"b", "c"}          ), 23.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"b", "d"}          ), 24.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"c", "d"}          ), 34.0  );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "b", "c"}     ), 123.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "b", "d"}     ), 124.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "c", "d"}     ), 134.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"b", "c", "d"}     ), 234.0 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "b", "c", "d"}), 1234.0);
   
   //Polynomial duplicate key
   if (bpm.GetVartype() == cimod::Vartype::SPIN) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "a", "a"}            ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "a", "a", "b"}        ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "c", "c", "c"}        ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"c", "b", "c", "b", "c", "b"}), 23.0);
   }
   else if (bpm.GetVartype() == cimod::Vartype::BINARY) {
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "a", "a", "a"}            ), 1.0 );
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "a", "a", "b", "b"}        ), 12.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"a", "c", "c", "c", "a"}        ), 13.0);
      EXPECT_DOUBLE_EQ(bpm.GetPolynomial({"c", "b", "c", "b", "c", "b", "b"}), 23.0);
   }
   
   //variables_to_integers
   EXPECT_EQ(bpm.GetVariablesToIntegers("a"), 0);
   EXPECT_EQ(bpm.GetVariablesToIntegers("b"), 1);
   EXPECT_EQ(bpm.GetVariablesToIntegers("c"), 2);
   EXPECT_EQ(bpm.GetVariablesToIntegers("d"), 3);
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a"}               ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"b"}               ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"c"}               ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"d"}               ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "b"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"b", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"b", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"c", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "b", "c"}     ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "b", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"b", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<std::string>{"a", "b", "c", "d"}), 1);

   //Polynomial Value
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.GetValueList()));

   //sorted_variables
   auto sorted_variables = bpm.GetSortedVariables();
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
   
   EXPECT_EQ(bpm.GetNumVariables(), 0);

   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 0.0);
   
   EXPECT_EQ(bpm.GetNumInteractions(), 0);
   
   EXPECT_EQ(bpm.GetDegree(), 0);
      
   //Polynomial
   EXPECT_EQ(bpm.GetPolynomial().size(), 0);
   
   //variables_to_integers
   EXPECT_EQ(bpm.GetVariablesToIntegers().size(), 0);
   
   //Polynomial Key
   EXPECT_EQ(bpm.GetKeyList().size(), 0);
   
   //Polynomial Val
   EXPECT_EQ(bpm.GetValueList().size(), 0);

   //sorted_variables
   auto sorted_variables = bpm.GetSortedVariables();
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
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMUINT(bpm);
}

TEST(ConstructionBPM, PolyMapINT) {
   BinaryPolynomialModel<int32_t, double> bpm(GeneratePolynomialINT(), Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMINT(bpm);
}

TEST(ConstructionBPM, PolyMapString) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
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
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
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
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
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
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMString(bpm);
}

TEST(AddInteractionBPM, basic) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.AddInteraction({3}, 3.0);
   bpm.AddInteraction({4}, 4.0);
   bpm.AddInteraction({1, 3}, 13.0);
   bpm.AddInteraction({1, 4}, 14.0);
   bpm.AddInteraction({2, 3}, 23.0);
   bpm.AddInteraction({2, 4}, 24.0);
   bpm.AddInteraction({3, 4}, 34.0);
   bpm.AddInteraction({1, 2, 4}, 124.0);
   bpm.AddInteraction({1, 3, 4}, 134.0);
   bpm.AddInteraction({2, 3, 4}, 234.0);
   bpm.AddInteraction({1, 2, 3, 4}, 1234.0);
   
   StateTestBPMUINT(bpm);

}

TEST(AddInteractionBPM, self_loop_SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.AddInteraction({3, 3, 3}, 3.0);
   bpm.AddInteraction({4, 4, 4, 4, 4}, 4.0);
   bpm.AddInteraction({1, 3, 1, 3, 3, 1}, 13.0);
   bpm.AddInteraction({1, 1, 1, 4, 4, 4}, 14.0);
   bpm.AddInteraction({3, 3, 3, 2, 2, 2}, 23.0);
   bpm.AddInteraction({2, 4}, 24.0);
   bpm.AddInteraction({3, 4}, 34.0);
   bpm.AddInteraction({1, 2, 4}, 124.0);
   bpm.AddInteraction({1, 3, 4}, 134.0);
   bpm.AddInteraction({2, 3, 4, 2, 4, 3, 3, 4, 2}, 234.0);
   bpm.AddInteraction({1, 2, 3, 4}, 1234.0);
   
   StateTestBPMUINT(bpm);

}

TEST(AddInteractionBPM, self_loop_BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::BINARY);
   
   bpm.AddInteraction({3, 3}, 3.0);
   bpm.AddInteraction({4, 4, 4}, 4.0);
   bpm.AddInteraction({1, 1, 3}, 13.0);
   bpm.AddInteraction({1, 4, 4, 4}, 14.0);
   bpm.AddInteraction({3, 3, 2}, 23.0);
   bpm.AddInteraction({2, 4}, 24.0);
   bpm.AddInteraction({3, 4}, 34.0);
   bpm.AddInteraction({1, 2, 4}, 124.0);
   bpm.AddInteraction({1, 3, 4}, 134.0);
   bpm.AddInteraction({2, 3, 4}, 234.0);
   bpm.AddInteraction({1, 2, 3, 3, 4, 1}, 1234.0);
   
   StateTestBPMUINT(bpm);

}

TEST(AddInteractionBPM, duplicate_value_1) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second /= 2;
   }
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.AddInteraction(it.first, it.second);
   };
   
   StateTestBPMUINT(bpm);

}

TEST(AddInteractionBPM, duplicate_value_2) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.AddInteraction(it.first, -it.second);
   };
   
   StateTestBPMEmpty(bpm);
}

TEST(AddInteractionsFromBPM, PolyMap) {
      
   BinaryPolynomialModel<uint32_t, double> bpm({}, Vartype::SPIN);
   
   bpm.AddInteractionsFrom(GeneratePolynomialUINT());
   
   StateTestBPMUINT(bpm);
   
}

TEST(AddInteractionsFromBPM, PolyKeyValue) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
   
   BinaryPolynomialModel<uint32_t, double> bpm({}, Vartype::SPIN);
   
   bpm.AddInteractionsFrom(poly_key, poly_value);
   
   StateTestBPMUINT(bpm);
   
}

TEST(AddOffsetBPM, basic) {
   Polynomial<uint32_t, double> polynomial;
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.AddOffset(3.0);
   
   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 3.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}) , 3.0);
   
   bpm.AddOffset(3.0);

   EXPECT_DOUBLE_EQ(bpm.GetOffset(), 6.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}) , 6.0);
   
   bpm.AddOffset(-6.0);
   
   StateTestBPMEmpty(bpm);
   
}

TEST(RemoveInteractionBPM, basic) {
   
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
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}          ), -3.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({7}         ), -2.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 11}     ), -9.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({11, 12, 14}), -1.0);
   
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{}          ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{7}         ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 11}     ), 1);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{11, 12, 14}), 1);
   
   EXPECT_TRUE(EXPECT_CONTAIN(-3.0, bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(-2.0, bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(-9.0, bpm.GetValueList()));
   EXPECT_TRUE(EXPECT_CONTAIN(-1.0, bpm.GetValueList()));

   bpm.RemoveInteraction({}          );
   bpm.RemoveInteraction({7}         );
   bpm.RemoveInteraction({2, 11}     );
   bpm.RemoveInteraction({11, 12, 14});
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}          ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({7}         ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 11}     ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({11, 12, 14}), 0.0);
   
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm.GetValueList()));
   
   StateTestBPMUINT(bpm);
   
}

TEST(RemoveInteractionBPM, remove_all) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   for (const auto &it: polynomial) {
      bpm.RemoveInteraction(it.first);
   };
   
   StateTestBPMEmpty(bpm);
   
}

TEST(RemoveInteractionBPM, self_loop_SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.RemoveInteraction({1, 1, 1});
   bpm.RemoveInteraction({2, 2, 2, 2, 2});
   bpm.RemoveInteraction({1, 1, 1, 2, 2, 2, 2, 2});
   bpm.RemoveInteraction({1, 1, 1, 2, 2, 2, 3, 3, 3});
   
   StateTestBPMEmpty(bpm);
}

TEST(RemoveInteractionBPM, self_loop_BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{1}, 1.0}, {{2}, 2.0},
      {{1, 2}, 12.0},
      {{1, 2, 3}, 123.0}
   };
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::BINARY);
   
   bpm.RemoveInteraction({1, 1});
   bpm.RemoveInteraction({2, 2, 2, 2});
   bpm.RemoveInteraction({1, 1, 1, 2, 2, 2});
   bpm.RemoveInteraction({1, 1, 2, 2, 2, 3});
   
   StateTestBPMEmpty(bpm);
}

TEST(RemoveInteractionsFromBPM, basic) {
   
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
   
   bpm.RemoveInteractionsFrom(removed_key_list);
      
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}          ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({7}         ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 11}     ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({11, 12, 14}), 0.0);
   
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm.GetKeyList().begin(), bpm.GetKeyList().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm.GetValueList()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm.GetValueList()));
   
   StateTestBPMUINT(bpm);
   
}

TEST(RemoveInteractionsFromBPM, remove_all) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   PolynomialKeyList<uint32_t> removed_key_list;
   
   for (const auto &it: polynomial) {
      removed_key_list.push_back(it.first);
   }
   
   bpm.RemoveInteractionsFrom(removed_key_list);
   
   StateTestBPMEmpty(bpm);
   
}

TEST(RemoveOffsetBPM, basic) {
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
   
   bpm.RemoveOffset();
   
   StateTestBPMUINT(bpm);
}

TEST(EnergyBPM, SPIN) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   Sample<uint32_t> sample_variables_spin_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_2{{0, +1}, {1, -1}, {2, +1}};
   Sample<uint32_t> sample_variables_spin_3{{0, -1}, {1, -1}, {2, -1}};

   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_spin_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_spin_2), -123.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_spin_3), -81.0 );
   
   std::vector<int32_t> sample_vec_variables_spin_1{+1, +1, +1};
   std::vector<int32_t> sample_vec_variables_spin_2{+1, -1, +1};
   std::vector<int32_t> sample_vec_variables_spin_3{-1, -1, -1};
   
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_spin_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_spin_2), -123.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_spin_3), -81.0 );

}

TEST(EnergyBPM, BINARY) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::BINARY);
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};

   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_3), 0.0   );
   
   std::vector<int32_t> sample_vec_variables_binary_1{+1, +1, +1};
   std::vector<int32_t> sample_vec_variables_binary_2{+1, +0, +1};
   std::vector<int32_t> sample_vec_variables_binary_3{+0, +0, +0};
   
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_vec_variables_binary_3), 0.0   );
      
}

TEST(EnergiesBPM, SPIN) {
   
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
   std::vector<double> en_vec = bpm.Energies(sample_variables_spin);
   
   EXPECT_DOUBLE_EQ(en_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec[1], -123.0);
   EXPECT_DOUBLE_EQ(en_vec[2], -81.0 );
   
   std::vector<std::vector<int32_t>> sample_vec_variables_spin {
      {+1, +1, +1},
      {+1, -1, +1},
      {-1, -1, -1}
   };
   
   std::vector<double> en_vec_vec = bpm.Energies(sample_vec_variables_spin);
   
   EXPECT_DOUBLE_EQ(en_vec_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[1], -123.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[2], -81.0 );
   
}

TEST(EnergiesBPM, BINARY) {
   
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
   
   std::vector<double> en_vec = bpm.Energies(sample_variables_binary);
   
   EXPECT_DOUBLE_EQ(en_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec[1], +24.0 );
   EXPECT_DOUBLE_EQ(en_vec[2], 0.0   );
   
   std::vector<std::vector<int32_t>> sample_vec_variables_binary {
      {+1, +1, +1},
      {+1, +0, +1},
      {+0, +0, +0}
   };
   
   std::vector<double> en_vec_vec = bpm.Energies(sample_vec_variables_binary);
   
   EXPECT_DOUBLE_EQ(en_vec_vec[0], +171.0);
   EXPECT_DOUBLE_EQ(en_vec_vec[1], +24.0 );
   EXPECT_DOUBLE_EQ(en_vec_vec[2], 0.0   );
   
}

TEST(ScaleBPM, all_scale) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second *= 2;
   }
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.Scale(0.5);
   
   StateTestBPMUINT(bpm);
   
}

TEST(ScaleBPM, ignored_interaction) {
   
   Polynomial<uint32_t, double> polynomial = GeneratePolynomialUINT();
   
   for (auto &&it: polynomial) {
      it.second *= 2;
   }
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.Scale(0.5, {{1,2}, {2, 4}, {1, 3, 4}});
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 2}   ), 12.0*2 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({2, 4}   ), 24.0*2 );
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({1, 3, 4}), 134.0*2);

   bpm.AddInteraction({1, 2}   , -12.0);
   bpm.AddInteraction({2, 4}   , -24.0);
   bpm.AddInteraction({1, 3, 4}, -134.0);
   
   StateTestBPMUINT(bpm);

}

TEST(ScaleBPM, ignored_offset) {
   
   Polynomial<uint32_t, double> polynomial {
      {{}, 100.0}
   };
      
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, Vartype::SPIN);
   
   bpm.Scale(0.5, {std::vector<uint32_t>{}});
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}), 100.0);
   
   bpm.Scale(0.5, {}, true);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}), 100.0);
   
   bpm.Scale(0.5, {{}}, true);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}), 100.0);
   
   bpm.Scale(0.5);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial({}), 50.0);
   
}

TEST(NormalizeBPM, all_normalize) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.normalize({-1, 1});
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({1}), 1.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({2}), 2.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 1}), 11.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 2}), 22.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({1, 2}), 12.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 1, 2}), 12.0/22.0);
   
}

TEST(NormalizeBPM, ignored_interaction) {
   
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, +12}
   };

   Vartype vartype = Vartype::SPIN;
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   bpm.normalize({-1, 1}, {{0, 1, 2}});
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({1}), 1.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({2}), 2.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 1}), 11.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 2}), 22.0/22.0);
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({1, 2}), 12.0/22.0);
   
   EXPECT_DOUBLE_EQ(bpm.GetPolynomial().at({0, 1, 2}), 12.0);
   
}

TEST(SerializableBPM, UINT) {
   BinaryPolynomialModel<uint32_t, double> bpm(GeneratePolynomialUINT(), Vartype::SPIN);
   BinaryPolynomialModel<uint32_t, double> bpm_from = BinaryPolynomialModel<uint32_t, double>::FromSerializable(bpm.ToSerializable());
   StateTestBPMUINT(bpm_from);
}

TEST(SerializableBPM, INT) {
   BinaryPolynomialModel<int32_t, double> bpm(GeneratePolynomialINT(), Vartype::SPIN);
   BinaryPolynomialModel<int32_t, double> bpm_from = BinaryPolynomialModel<int32_t, double>::FromSerializable(bpm.ToSerializable());
   StateTestBPMINT(bpm_from);
}

TEST(SerializableBPM, String) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   BinaryPolynomialModel<std::string, double> bpm_from = BinaryPolynomialModel<std::string, double>::FromSerializable(bpm.ToSerializable());
   StateTestBPMString(bpm_from);
}

TEST(SerializableBPM, StringToUINT) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   auto obj = bpm.ToSerializable();
   std::vector<std::size_t> string_to_num(obj["variables"].size());
   std::iota(string_to_num.begin(), string_to_num.end(), 1);
   obj["variables"] = string_to_num;
   BinaryPolynomialModel<uint32_t, double> bpm_from = BinaryPolynomialModel<uint32_t, double>::FromSerializable(obj);
   StateTestBPMUINT(bpm_from);
}

TEST(SerializableBPM, StringToINT) {
   BinaryPolynomialModel<std::string, double> bpm(GeneratePolynomialString(), Vartype::SPIN);
   auto obj = bpm.ToSerializable();
   std::vector<std::size_t> string_to_num(4);
   string_to_num[0] = -1;
   string_to_num[1] = -2;
   string_to_num[2] = -3;
   string_to_num[3] = -4;
   obj["variables"] = string_to_num;
   BinaryPolynomialModel<int32_t, double> bpm_from = BinaryPolynomialModel<int32_t, double>::FromSerializable(obj);
   StateTestBPMINT(bpm_from);
}

TEST(FromHubo, MapUINT) {
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHubo(GeneratePolynomialUINT());
   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   StateTestBPMUINT(bpm);
}

TEST(FromHubo, MapINT) {
   auto bpm = BinaryPolynomialModel<int32_t, double>::FromHubo(GeneratePolynomialINT());
   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   StateTestBPMINT(bpm);
}

TEST(FromHubo, MapString) {
   auto bpm = BinaryPolynomialModel<std::string, double>::FromHubo(GeneratePolynomialString());
   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   StateTestBPMString(bpm);
}

TEST(FromHubo, KeyValueUINT) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   
   StateTestBPMUINT(bpm);
}

TEST(FromHubo, KeyValueINT) {
   PolynomialKeyList<int32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<int32_t, double>::FromHubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   
   StateTestBPMINT(bpm);
}

TEST(FromHubo, KeyValueString) {
   PolynomialKeyList<std::string> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialString()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<std::string, double>::FromHubo(poly_key, poly_value);

   EXPECT_EQ(Vartype::BINARY, bpm.GetVartype());
   
   StateTestBPMString(bpm);
}

TEST(FromHising, MapUINT) {
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHising(GeneratePolynomialUINT());
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMUINT(bpm);
}

TEST(FromHising, MapINT) {
   auto bpm = BinaryPolynomialModel<int32_t, double>::FromHising(GeneratePolynomialINT());
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMINT(bpm);
}

TEST(FromHising, MapString) {
   auto bpm = BinaryPolynomialModel<std::string, double>::FromHising(GeneratePolynomialString());
   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   StateTestBPMString(bpm);
}

TEST(FromHising, KeyValueUINT) {
   
   PolynomialKeyList<uint32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialUINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   
   StateTestBPMUINT(bpm);
}

TEST(FromHising, KeyValueINT) {
   PolynomialKeyList<int32_t> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialINT()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<int32_t, double>::FromHising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   
   StateTestBPMINT(bpm);
}

TEST(FromHising, KeyValueString) {
   PolynomialKeyList<std::string> poly_key;
   PolynomialValueList<double> poly_value;
   
   for (const auto &it: GeneratePolynomialString()) {
      poly_key.push_back(it.first);
      poly_value.push_back(it.second);
   }
      
   auto bpm = BinaryPolynomialModel<std::string, double>::FromHising(poly_key, poly_value);

   EXPECT_EQ(Vartype::SPIN, bpm.GetVartype());
   
   StateTestBPMString(bpm);
}

TEST(ClearBPM, basic) {
   Polynomial<uint32_t, double> polynomial {
      {{0}, 0.0}, {{1}, 1.0}, {{2}, 2.0},
      {{0, 1}, 11.0}, {{0, 2}, 22.0}, {{1, 2}, 12.0},
      {{0, 1, 2}, 123.0}
   };
   
   Vartype vartype = Vartype::BINARY;

   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);

   bpm.Clear();
   
   Sample<uint32_t> sample_variables_binary_1{{0, +1}, {1, +1}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_2{{0, +1}, {1, +0}, {2, +1}};
   Sample<uint32_t> sample_variables_binary_3{{0, +0}, {1, +0}, {2, +0}};
   
   EXPECT_TRUE(bpm.GetPolynomial().empty());
   EXPECT_TRUE(bpm.GetVariables().empty());
   EXPECT_EQ(bpm.GetVartype(), Vartype::BINARY);
   
   //Chech if the methods in Binary Polynomial Model work properly after executing empty()
   EXPECT_EQ(bpm.GetNumVariables(), 0);
   bpm.RemoveVariable(1);
   bpm.RemoveVariablesFrom(std::vector<uint32_t>{1,2,3,4,5});
   bpm.RemoveInteraction(std::vector<uint32_t>{1,2});
   bpm.RemoveInteractionsFrom(std::vector<std::vector<uint32_t>>{{1,2},{1,3},{1,4}});
   bpm.Scale(1.0);
   bpm.normalize();
   
   //energy
   EXPECT_THROW(bpm.Energy(sample_variables_binary_1), std::runtime_error);
   EXPECT_THROW(bpm.Energy(sample_variables_binary_2), std::runtime_error);
   EXPECT_THROW(bpm.Energy(sample_variables_binary_3), std::runtime_error);
   
   //Reset polynomial model
   bpm.AddInteraction({0, 1}   , 11.0, Vartype::BINARY);
   bpm.AddInteraction({0, 2}   , 22.0);
   bpm.AddInteraction({1, 2}   , 12.0);
   bpm.AddInteraction({0, 1, 2}, 123.0);
   
   bpm.AddInteraction({0}, 0.0);
   bpm.AddInteraction({1}, 1.0);
   bpm.AddInteraction({2}, 2.0);
   
   //Check energy
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_1), +171.0);
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_2), +24.0 );
   EXPECT_DOUBLE_EQ(bpm.Energy(sample_variables_binary_3), 0.0   );
   
   EXPECT_EQ(bpm.GetNumVariables(), 3);
   
   EXPECT_EQ(bpm.GetVariables().count(0), 1);
   EXPECT_EQ(bpm.GetVariables().count(1), 1);
   EXPECT_EQ(bpm.GetVariables().count(2), 1);
   
   for (const auto &it: bpm.GetPolynomial()) {
      EXPECT_DOUBLE_EQ(it.second, polynomial[it.first]);
   }
   
   EXPECT_EQ(bpm.GetVartype(), Vartype::BINARY);
   
}

TEST(VartypeBPM, SpinBinarySpin) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHising(GeneratePolynomialUINT());
   auto bpm_binary = BinaryPolynomialModel<uint32_t, double>(bpm.ToHubo(), Vartype::BINARY);
   auto bpm_ising  = BinaryPolynomialModel<uint32_t, double>(bpm_binary.ToHising(), Vartype::SPIN);
   
   StateTestBPMUINT(bpm_ising);
   
}

TEST(VartypeBPM, BinarySPINBinary) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHubo(GeneratePolynomialUINT());
   auto bpm_ising  = BinaryPolynomialModel<uint32_t, double>(bpm.ToHising(), Vartype::SPIN);
   auto bpm_bianry = BinaryPolynomialModel<uint32_t, double>(bpm_ising.ToHubo(), Vartype::BINARY);
   
   StateTestBPMUINT(bpm_bianry);
   
}

TEST(VartypeBPM, ChangeVartypeSpinBinarySpin) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHising(GeneratePolynomialUINT());
   auto bpm_binary = bpm.ChangeVartype(Vartype::BINARY, true);
   auto bpm_ising  = bpm_binary.ChangeVartype(Vartype::SPIN, true);
   
   StateTestBPMUINT(bpm_ising);
   bpm.ChangeVartype(Vartype::SPIN);
   StateTestBPMUINT(bpm);
   
}

TEST(VartypeBPM, ChangeVartypeBinarySPINBinary) {
   
   auto bpm = BinaryPolynomialModel<uint32_t, double>::FromHubo(GeneratePolynomialUINT());
   auto bpm_ising  = bpm.ChangeVartype(Vartype::SPIN, true);
   auto bpm_binary = bpm_ising.ChangeVartype(Vartype::BINARY, true);
   
   StateTestBPMUINT(bpm_binary);
   bpm.ChangeVartype(Vartype::BINARY);
   StateTestBPMUINT(bpm);
   
}



}

