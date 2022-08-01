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

#pragma once

#include <gtest/gtest.h>

#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdint>
#include <string>
#include <iostream>
#include <tuple>

#include <nlohmann/json.hpp>

#include <cimod/binary_quadratic_model.hpp>

using json = nlohmann::json;
using namespace cimod;

template<typename IndexType, typename FloatType, typename DataType>
using BQM = BinaryQuadraticModel<IndexType, FloatType, DataType>;

template<typename DataType>
struct BQMTester{
    static void test_DenseConstructionTest_Construction()
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
        Quadratic<uint32_t, double> quadratic
        {
            {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
                {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
                {std::make_pair(3, 4), 34.0}
        };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BQM<uint32_t, double, DataType> bqm_k4(linear, quadratic, offset, vartype);

        Linear<uint32_t, double> bqm_linear = bqm_k4.get_linear();
        Quadratic<uint32_t, double> bqm_quadratic = bqm_k4.get_quadratic();
        double bqm_offset = bqm_k4.get_offset();
        Vartype bqm_vartype = bqm_k4.get_vartype();

        // check linear
        for(auto &it : bqm_linear)
        {
            EXPECT_EQ(it.second, linear[it.first]);
        }
        // check quadratic
        for(auto &it : bqm_quadratic)
        {
            EXPECT_EQ(it.second, quadratic[it.first]);
        }
        // check offset
        EXPECT_EQ(offset, bqm_offset);
        // check vartype
        EXPECT_EQ(vartype, bqm_vartype);

        //reorder test
        linear = Linear<uint32_t, double>{};
        quadratic = Quadratic<uint32_t, double>{{std::make_pair(3, 4), 34.0}, {std::make_pair(4, 3), 34.0}};
        bqm_k4 = BQM<uint32_t, double, DataType>(linear, quadratic, offset, vartype);
        bqm_quadratic = bqm_k4.get_quadratic();
        std::pair<uint32_t, uint32_t> indexpair = std::make_pair(3,4);
        for(auto &it : bqm_quadratic)
        {
            EXPECT_EQ(it.first , indexpair);
            EXPECT_DOUBLE_EQ(it.second, 68.0);
        }

        //another constructor
        bqm_k4 = BQM<uint32_t, double, DataType>(linear, quadratic, vartype);
        EXPECT_EQ(offset, 0.0);
    }

    static void test_DenseConstructionTest_ConstructionString()
    {
        Linear<std::string, double> linear{ {"a", 1.0}, {"b", 2.0}, {"c", 3.0}, {"d", 4.0} };
        Quadratic<std::string, double> quadratic
        {
            {std::make_pair("a", "b"), 12.0}, {std::make_pair("a", "c"), 13.0}, {std::make_pair("a", "d"), 14.0},
                {std::make_pair("b", "c"), 23.0}, {std::make_pair("b", "d"), 24.0},
                {std::make_pair("c", "d"), 34.0}
        };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BQM<std::string, double, DataType> bqm_k4(linear, quadratic, offset, vartype);

        //bqm_k4.print();

        Linear<std::string, double> bqm_linear = bqm_k4.get_linear();
        Quadratic<std::string, double> bqm_quadratic = bqm_k4.get_quadratic();
        double bqm_offset = bqm_k4.get_offset();
        Vartype bqm_vartype = bqm_k4.get_vartype();

        // check linear
        for(auto &it : bqm_linear)
        {
            EXPECT_EQ(it.second, linear[it.first]);
        }
        // check quadratic
        for(auto &it : bqm_quadratic)
        {
            EXPECT_EQ(it.second, quadratic[it.first]);
        }
        // check offset
        EXPECT_EQ(offset, bqm_offset);
        // check vartype
        EXPECT_EQ(vartype, bqm_vartype);
    }

    static void test_DenseConstructionTest_ConstructionMatrix()
    {
        using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        Matrix m = Matrix::Zero(5,5);
        m(0,1) = 1;
        m(0,2) = 2;
        m(0,3) = 3;
        m(0,4) = 4;
        m(1,0) = 5;
        m(1,2) = 6;
        m(1,3) = 7;
        m(1,4) = 8;
        m(2,0) = 9;
        m(2,1) = 10;
        m(2,3) = 11;
        m(2,4) = 12;
        m(3,0) = 13;
        m(3,1) = 14;
        m(3,2) = 15;
        m(3,4) = 16;
        m(4,0) = 17;
        m(4,1) = 18;
        m(4,2) = 19;
        m(4,3) = 20;

        auto labels = std::vector<std::string>{"a", "b", "d", "e"};

        BQM<std::string, double, DataType> bqm(m, labels, 0.0, Vartype::SPIN);

        Matrix int_mat = bqm.interaction_matrix();
        EXPECT_DOUBLE_EQ(int_mat(0,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(0,1), 6);
        EXPECT_DOUBLE_EQ(int_mat(0,2), 11);
        EXPECT_DOUBLE_EQ(int_mat(0,3), 16);
        EXPECT_DOUBLE_EQ(int_mat(0,4), 21);
        EXPECT_DOUBLE_EQ(int_mat(1,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(1,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(1,2), 16);
        EXPECT_DOUBLE_EQ(int_mat(1,3), 21);
        EXPECT_DOUBLE_EQ(int_mat(1,4), 26);
        EXPECT_DOUBLE_EQ(int_mat(2,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,3), 26);
        EXPECT_DOUBLE_EQ(int_mat(2,4), 31);
        EXPECT_DOUBLE_EQ(int_mat(3,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,3), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,4), 36);
        EXPECT_DOUBLE_EQ(int_mat(4,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,3), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,4), 1);

        //check elements

        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "b"), 6);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "d"), 11);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "e"), 16);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("b", "d"), 16);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("b", "e"), 21);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("d", "e"), 26);

        EXPECT_DOUBLE_EQ(bqm.get_linear("a"), 21);
        EXPECT_DOUBLE_EQ(bqm.get_linear("b"), 26);
        EXPECT_DOUBLE_EQ(bqm.get_linear("d"), 31);
        EXPECT_DOUBLE_EQ(bqm.get_linear("e"), 36);
    }

    static void test_DenseConstructionTest_ConstructionMatrix2()
    {
        using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
        Matrix m = Matrix::Zero(4,4);
        m(0,0) = 21;
        m(1,1) = 26;
        m(2,2) = 31;
        m(3,3) = 36;

        m(0,1) = 1;
        m(0,2) = 2;
        m(0,3) = 3;

        m(1,0) = 5;
        m(1,2) = 6;
        m(1,3) = 7;

        m(2,0) = 9;
        m(2,1) = 10;
        m(2,3) = 11;

        m(3,0) = 13;
        m(3,1) = 14;
        m(3,2) = 15;

        auto labels = std::vector<std::string>{"a", "b", "d", "e"};

        BQM<std::string, double, DataType> bqm(m, labels, 0.0, Vartype::SPIN);

        Matrix int_mat = bqm.interaction_matrix();
        EXPECT_DOUBLE_EQ(int_mat(0,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(0,1), 6);
        EXPECT_DOUBLE_EQ(int_mat(0,2), 11);
        EXPECT_DOUBLE_EQ(int_mat(0,3), 16);
        EXPECT_DOUBLE_EQ(int_mat(0,4), 21);
        EXPECT_DOUBLE_EQ(int_mat(1,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(1,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(1,2), 16);
        EXPECT_DOUBLE_EQ(int_mat(1,3), 21);
        EXPECT_DOUBLE_EQ(int_mat(1,4), 26);
        EXPECT_DOUBLE_EQ(int_mat(2,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(2,3), 26);
        EXPECT_DOUBLE_EQ(int_mat(2,4), 31);
        EXPECT_DOUBLE_EQ(int_mat(3,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,3), 0);
        EXPECT_DOUBLE_EQ(int_mat(3,4), 36);
        EXPECT_DOUBLE_EQ(int_mat(4,0), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,1), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,2), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,3), 0);
        EXPECT_DOUBLE_EQ(int_mat(4,4), 1);

        //check elements

        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "b"), 6);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "d"), 11);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("a", "e"), 16);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("b", "d"), 16);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("b", "e"), 21);
        EXPECT_DOUBLE_EQ(bqm.get_quadratic("d", "e"), 26);

        EXPECT_DOUBLE_EQ(bqm.get_linear("a"), 21);
        EXPECT_DOUBLE_EQ(bqm.get_linear("b"), 26);
        EXPECT_DOUBLE_EQ(bqm.get_linear("d"), 31);
        EXPECT_DOUBLE_EQ(bqm.get_linear("e"), 36);
    }

    static void test_DenseBQMFunctionTest_add_variable()
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 2), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        // check length
        EXPECT_EQ(bqm.get_num_variables(), 2);

        bqm.add_variable(2, 2.0);
        bqm.add_variable(1, 0.33);
        bqm.add_variable(0, 0.33);

        Linear<uint32_t, double> bqm_linear = bqm.get_linear();

        // check length
        EXPECT_EQ(bqm.get_num_variables(), 3);
        // check linear
        EXPECT_EQ(bqm_linear[0], 0.33);
        EXPECT_EQ(bqm_linear[1], 0.33);
        EXPECT_EQ(bqm_linear[2], 3.0);


        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_add_variables_from()
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        // check length
        EXPECT_EQ(bqm.get_num_variables(), 0);

        Linear<uint32_t, double> linear2 = { {0, 0.5}, {1, -1.0} };
        Linear<uint32_t, double> linear3 = { {1, -1.0}, {2, 2.0} };

        bqm.add_variables_from(linear2);
        // check variable
        EXPECT_EQ(bqm.contains(1), true);

        bqm.add_variables_from(linear3);

        // check bias
        Linear<uint32_t, double> bqm_linear = bqm.get_linear();
        EXPECT_EQ(bqm_linear[1], -2.0);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_add_interaction()
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        //bqm.print();

        bqm.add_interaction(0, 2, 2);
        bqm.add_interaction(1, 0, 0.25);
        bqm.add_interaction(1, 2, 0.25);

        EXPECT_EQ(bqm.get_num_variables(), 3);

        // check bias
        Quadratic<uint32_t, double> bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], 0.75);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_add_interactions_from()
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        Quadratic<uint32_t, double> quadratic1{ {std::make_pair(0, 1), -0.5} };
        bqm.add_interactions_from(quadratic1);

        // check bias
        Quadratic<uint32_t, double> bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], -0.5);

        Quadratic<uint32_t, double> quadratic2{ {std::make_pair(0, 1), -0.5}, {std::make_pair(0, 2), 2.0} };
        Quadratic<uint32_t, double> quadratic3{ {std::make_pair(1, 2), 2.0} };
        bqm.add_interactions_from(quadratic2);
        bqm.add_interactions_from(quadratic3);

        // check length
        EXPECT_EQ(bqm.get_num_variables(), 3);

        // check bias
        bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], -1.0);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_add_offset()
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.add_offset(1.0);

        // check offset
        EXPECT_EQ(bqm.get_offset(), 0.5);
    }

    static void test_DenseBQMFunctionTest_energy()
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        // check energy
        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        EXPECT_DOUBLE_EQ(bqm.energy(sample1), -0.5);
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        EXPECT_DOUBLE_EQ(bqm.energy(sample2), 3.5);
    }

    static void test_DenseBQMFunctionTest_energies()
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        std::vector<Sample<uint32_t>> samples{sample1, sample2};
        std::vector<double> en_vec = bqm.energies(samples);

        // check energies
        EXPECT_DOUBLE_EQ(en_vec[0], -0.5);
        EXPECT_DOUBLE_EQ(en_vec[1], 3.5);
    }

    static void test_DenseBQMFunctionTest_empty()
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };

        Linear<uint32_t, double> empty_linear = Linear<uint32_t, double>();
        Quadratic<uint32_t, double> empty_quadratic = Quadratic<uint32_t, double>();


        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        BQM<uint32_t, double, DataType> bqm_s = bqm.empty(Vartype::SPIN);

        EXPECT_EQ(bqm.get_linear(), linear);
        EXPECT_EQ(bqm.get_quadratic(), quadratic);

        EXPECT_EQ(bqm_s.get_linear(), empty_linear);
        EXPECT_EQ(bqm_s.get_quadratic(), empty_quadratic);
        EXPECT_EQ(bqm_s.get_vartype(), Vartype::SPIN);

        BQM<uint32_t, double, DataType> bqm_b = bqm.empty(Vartype::BINARY);

        EXPECT_EQ(bqm.get_linear(), linear);
        EXPECT_EQ(bqm.get_quadratic(), quadratic);

        EXPECT_EQ(bqm_b.get_linear(), empty_linear);
        EXPECT_EQ(bqm_b.get_quadratic(), empty_quadratic);
        EXPECT_EQ(bqm_b.get_vartype(), Vartype::BINARY);

    }

    static void test_DenseBQMFunctionTest_to_qubo()
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        auto t_qubo = bqm.to_qubo();
        Quadratic<uint32_t, double> Q = std::get<0>(t_qubo);
        double offset_qubo = std::get<1>(t_qubo);

        // check quadratic matrix and offset
        EXPECT_DOUBLE_EQ(Q[std::make_pair(0, 0)], 1.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(0, 1)], 2.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(1, 1)], -6.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(1, 2)], 6.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(2, 2)], -2.0);
        EXPECT_DOUBLE_EQ(offset_qubo, 2.9);
    }

    static void test_DenseBQMFunctionTest_to_ising()
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        auto t_ising = bqm.to_ising();
        Linear<uint32_t, double> h = std::get<0>(t_ising);
        Quadratic<uint32_t, double> J = std::get<1>(t_ising);
        double offset_ising = std::get<2>(t_ising);

        // check biases
        EXPECT_DOUBLE_EQ(h[0], 1.0);
        EXPECT_DOUBLE_EQ(J[std::make_pair(0, 1)], 0.5);
        EXPECT_DOUBLE_EQ(h[1], -1.0);
        EXPECT_DOUBLE_EQ(J[std::make_pair(1, 2)], 1.5);
        EXPECT_DOUBLE_EQ(h[2], 0.5);
        EXPECT_DOUBLE_EQ(offset_ising, 1.4);
    }

    static void test_DenseBQMFunctionTest_from_qubo()
    {
        Quadratic<uint32_t, double> origQ{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5}, {std::make_pair(0, 0), 1.0}, {std::make_pair(1, 1), -1.0}, {std::make_pair(2, 2), 0.5}};
        double offset = 1.4;

        auto bqm = BQM<uint32_t, double, DataType>::from_qubo(origQ, offset);

        Linear<uint32_t, double> linear = bqm.get_linear();
        Quadratic<uint32_t, double> Q = bqm.get_quadratic();
        double offset_ising = bqm.get_offset();

        // check quadratic matrix and offset
        EXPECT_DOUBLE_EQ(linear[0], 1.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(0, 1)], 0.5);
        EXPECT_DOUBLE_EQ(linear[1], -1.0);
        EXPECT_DOUBLE_EQ(Q[std::make_pair(1, 2)], 1.5);
        EXPECT_DOUBLE_EQ(linear[2], 0.5);
        EXPECT_DOUBLE_EQ(offset_ising, 1.4);
    }

    static void test_DenseBQMFunctionTest_from_ising()
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;

        auto bqm = BQM<uint32_t, double, DataType>::from_ising(linear, quadratic, offset);
        Linear<uint32_t, double> h = bqm.get_linear();
        Quadratic<uint32_t, double> J = bqm.get_quadratic();
        double offset_ising = bqm.get_offset();

        // check biases
        EXPECT_DOUBLE_EQ(h[0], 1.0);
        EXPECT_DOUBLE_EQ(J[std::make_pair(0, 1)], 0.5);
        EXPECT_DOUBLE_EQ(h[1], -1.0);
        EXPECT_DOUBLE_EQ(J[std::make_pair(1, 2)], 1.5);
        EXPECT_DOUBLE_EQ(h[2], 0.5);
        EXPECT_DOUBLE_EQ(offset_ising, 1.4);
    }

    static void test_DenseBQMFunctionTest_remove_variable()
    {
        Linear<std::string, double> linear{ {"a", 0.0}, {"b", 1.0}, {"c", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), 0.25}, {std::make_pair("a", "c"), 0.5}, {std::make_pair("b", "c"), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);
        //bqm.print();
        bqm.remove_variable("a");
        //bqm.print();

        // check variables
        EXPECT_EQ(bqm.contains("a"), false);
        EXPECT_EQ(bqm.contains("b"), true);
        EXPECT_EQ(bqm.contains("c"), true);

        EXPECT_EQ(bqm.get_linear().size(), 2);
        EXPECT_EQ(bqm.get_quadratic().size(), 1);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_remove_variables_from()
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.25}, {std::make_pair(0, 2), 0.5}, {std::make_pair(1, 2), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);
        //bqm.print();

        std::vector<uint32_t> variables = {0, 1};

        bqm.remove_variables_from(variables);
        //bqm.print();

        // check variables
        EXPECT_EQ(bqm.contains(0), false);
        EXPECT_EQ(bqm.contains(1), false);
        EXPECT_EQ(bqm.contains(2), true);

        EXPECT_EQ(bqm.get_linear().size(), 1);
        EXPECT_EQ(bqm.get_quadratic().size(), 0);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_remove_interaction()
    {
        Linear<std::string, double> linear{ {"c", 2} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}, {std::make_pair("b", "c"), 1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);
        bqm.remove_interaction("b", "c");

        EXPECT_EQ(bqm.contains("a"), true);
        EXPECT_EQ(bqm.contains("b"), true);
        EXPECT_EQ(bqm.contains("c"), true);

        // the number of elements of m_linear depends on the system
        //EXPECT_EQ(bqm.get_linear().size(), 1);
        std::cout << bqm.interaction_matrix() << std::endl;
        EXPECT_EQ(bqm.get_quadratic().size(), 1);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);

        bqm.remove_interaction("a", "b");

        EXPECT_EQ(bqm.contains("a"), false);
        EXPECT_EQ(bqm.contains("b"), false);
        EXPECT_EQ(bqm.contains("c"), true);

        // the number of elements of m_linear depends on the system
        //EXPECT_EQ(bqm.get_linear().size(), 1);
        EXPECT_EQ(bqm.get_quadratic().size(), 0);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    static void test_DenseBQMFunctionTest_scale()
    {
        Linear<std::string, double> linear{{"a", -2.0}, {"b", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}};
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.scale(0.5);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();
        double bqm_offset = bqm.get_offset();

        // check biases and offset
        EXPECT_DOUBLE_EQ(bqm_linear["a"], -1.0);
        EXPECT_DOUBLE_EQ(bqm_quadratic[std::make_pair("a", "b")], -0.5);
        EXPECT_DOUBLE_EQ(bqm_offset, 0.5);
    }

    static void test_DenseBQMFunctionTest_normalize()
    {
        Linear<std::string, double> linear{ {"a", -2.0}, {"b", 1.5} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.normalize(std::make_pair(-1.0, 1.0));
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        auto comp = [](const auto &a, const auto &b) { return std::abs(a.second) < std::abs(b.second); };
        auto lin_max = std::max_element(bqm_linear.begin(), bqm_linear.end(), comp);
        auto quad_max = std::max_element(bqm_quadratic.begin(), bqm_quadratic.end(), comp);

        // check maximum biases
        EXPECT_DOUBLE_EQ(lin_max->second, -1.0);
        EXPECT_DOUBLE_EQ(quad_max->second, -0.5);
    }

    static void test_DenseBQMFunctionTest_fix_variable()
    {
        Linear<std::string, double> linear{ {"a", -0.5}, {"b", 0.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.fix_variable("a", -1);
        auto bqm_linear = bqm.get_linear();
        double bqm_offset = bqm.get_offset();

        // check offset, linear bias and variable
        EXPECT_DOUBLE_EQ(bqm_offset, 0.5);
        EXPECT_DOUBLE_EQ(bqm_linear["b"], 1.0);
        EXPECT_EQ(bqm.contains("a"), false);
    }

    static void test_DenseBQMFunctionTest_flip_variable()
    {
        Linear<uint32_t, double> linear{{1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 0.5} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.flip_variable(1);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        // check biases
        EXPECT_DOUBLE_EQ(bqm_linear[1], -1.0);
        EXPECT_DOUBLE_EQ(bqm_linear[2], 2.0);
        EXPECT_DOUBLE_EQ(bqm_quadratic[std::make_pair(1, 2)], -0.5);
    }

    static void test_DenseBQMFunctionTest_flip_variable_binary()
    {
        Linear<uint32_t, double> linear{{1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 0.5} };
        double offset = 0.5;
        Vartype vartype = Vartype::BINARY;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        bqm.flip_variable(1);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        // check biases
        EXPECT_DOUBLE_EQ(bqm_linear[1], -1.0);
        EXPECT_DOUBLE_EQ(bqm_linear[2], 2.5);
        EXPECT_DOUBLE_EQ(bqm_quadratic[std::make_pair(1, 2)], -0.5);
        EXPECT_DOUBLE_EQ(bqm.get_offset(), 1.5);
    }

    // currently disabled
    //static void test_DenseBQMFunctionTest_contract_variables()
    //{
    //    Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
    //    Quadratic<uint32_t, double> quadratic
    //    {
    //        {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
    //        {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
    //        {std::make_pair(3, 4), 34.0}
    //    };
    //    double offset = 0.5;
    //    Vartype vartype = Vartype::SPIN;

    //    BQM<uint32_t, double> bqm(linear, quadratic, offset, vartype);
    //    bqm.contract_variables(2, 3);

    //    auto bqm_quadratic = bqm.get_quadratic();

    //    // check variable and quadratic bias
    //    EXPECT_EQ(bqm.contains(3), false);
    //    EXPECT_EQ(bqm_quadratic[std::make_pair(1, 2)], 25.0);
    //}

    static void test_DenseBQMFunctionTest_change_vartype()
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BQM<uint32_t, double, DataType> bqm(linear, quadratic, offset, vartype);

        auto bqm2 = bqm.change_vartype(Vartype::BINARY, true);
        auto lin = bqm.get_linear();
        auto quad = bqm.get_quadratic();
        auto lin2 = bqm2.get_linear();
        auto quad2 = bqm2.get_quadratic();

        // check quadratic matrix and offset
        EXPECT_EQ(quad[std::make_pair(0, 1)], 2.0);
        EXPECT_EQ(quad[std::make_pair(1, 2)], 6.0);
        EXPECT_EQ(lin[0], 1.0);
        EXPECT_EQ(lin[1], -6.0);
        EXPECT_EQ(lin[2], -2.0);
        EXPECT_EQ(bqm.get_offset(), 2.9);
        EXPECT_EQ(bqm.get_vartype(), Vartype::BINARY);

        EXPECT_EQ(quad2[std::make_pair(0, 1)], 2.0);
        EXPECT_EQ(quad2[std::make_pair(1, 2)], 6.0);
        EXPECT_EQ(lin2[0], 1.0);
        EXPECT_EQ(lin2[1], -6.0);
        EXPECT_EQ(lin2[2], -2.0);
        EXPECT_EQ(bqm2.get_offset(), 2.9);
        EXPECT_EQ(bqm2.get_vartype(), Vartype::BINARY);

        auto bqm3 = bqm.change_vartype(Vartype::SPIN, false);
        EXPECT_EQ(bqm.get_offset(), 2.9);
        EXPECT_EQ(bqm.get_vartype(), Vartype::BINARY);

    }

    static void test_DenseBQMFunctionTest_to_serializable()
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 5.0;
        Vartype vartype = Vartype::BINARY;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();

        std::cout << j << std::endl;
    }

    static void test_DenseBQMFunctionTest_from_serializable()
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 5.0;
        Vartype vartype = Vartype::BINARY;

        BQM<std::string, double, DataType> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();
        //std::cout << j << std::endl;

        BQM<std::string, double, DataType> bqm2 = BQM<std::string, double, DataType>::from_serializable(j);
        //bqm2.print();

        auto bqm_linear = bqm2.get_linear();
        auto bqm_quadratic = bqm2.get_quadratic();
        EXPECT_EQ(bqm2.get_offset(), bqm.get_offset());
        EXPECT_EQ(bqm_linear["c"], linear["c"]);
        EXPECT_EQ(bqm_quadratic[std::make_pair("b", "e")], quadratic[std::make_pair("b", "e")]);
    }
};
