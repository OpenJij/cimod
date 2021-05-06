#include "gtest/gtest.h"

#include "../src/binary_quadratic_model.hpp"
#include "../src/binary_polynomial_model.hpp"
#include "../src/binary_quadratic_model_dense.hpp"
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
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
        Quadratic<uint32_t, double> quadratic
        {
            {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
            {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
            {std::make_pair(3, 4), 34.0}
        };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm_k4(linear, quadratic, offset, vartype);

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
        bqm_k4 = BinaryQuadraticModel_Dense<uint32_t, double>(linear, quadratic, offset, vartype);
        bqm_quadratic = bqm_k4.get_quadratic();
        std::pair<uint32_t, uint32_t> indexpair = std::make_pair(3,4);
        for(auto &it : bqm_quadratic)
        {
            EXPECT_EQ(it.first , indexpair);
            EXPECT_DOUBLE_EQ(it.second, 68.0);
        }

        //another constructor
        bqm_k4 = BinaryQuadraticModel_Dense<uint32_t, double>(linear, quadratic, vartype);
        EXPECT_EQ(offset, 0.0);
    }

    TEST(DenseConstructionTest, ConstructionString)
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

        BinaryQuadraticModel_Dense<std::string, double> bqm_k4(linear, quadratic, offset, vartype);

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

    TEST(DenseConstructionTest, ConstructionMatrix)
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

        BinaryQuadraticModel_Dense<std::string, double> bqm(m, labels, 0.0, Vartype::SPIN);

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

    TEST(DenseConstructionTest, ConstructionMatrix2)
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

        BinaryQuadraticModel_Dense<std::string, double> bqm(m, labels, 0.0, Vartype::SPIN);

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

    TEST(DenseBQMFunctionTest, add_variable)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 2), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, add_variables_from)
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, add_interaction)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        //bqm.print();

        bqm.add_interaction(0, 2, 2);
        bqm.add_interaction(0, 1, 0.25);
        bqm.add_interaction(1, 2, 0.25);

        EXPECT_EQ(bqm.get_num_variables(), 3);

        // check bias
        Quadratic<uint32_t, double> bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], 0.75);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    TEST(DenseBQMFunctionTest, add_interactions_from)
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, add_offset)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        bqm.add_offset(1.0);
        
        // check offset
        EXPECT_EQ(bqm.get_offset(), 0.5);
    }

    TEST(DenseBQMFunctionTest, energy)
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        // check energy
        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        EXPECT_DOUBLE_EQ(bqm.energy(sample1), -0.5);
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        EXPECT_DOUBLE_EQ(bqm.energy(sample2), 3.5);
    }

    TEST(DenseBQMFunctionTest, energies)
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        std::vector<Sample<uint32_t>> samples{sample1, sample2};
        std::vector<double> en_vec = bqm.energies(samples);

        // check energies
        EXPECT_DOUBLE_EQ(en_vec[0], -0.5);
        EXPECT_DOUBLE_EQ(en_vec[1], 3.5);
    }

    TEST(DenseBQMFunctionTest, to_qubo)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, to_ising)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, from_qubo)
    {
        Quadratic<uint32_t, double> origQ{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5}, {std::make_pair(0, 0), 1.0}, {std::make_pair(1, 1), -1.0}, {std::make_pair(2, 2), 0.5}};
        double offset = 1.4;

        auto bqm = BinaryQuadraticModel_Dense<uint32_t, double>::from_qubo(origQ, offset);

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

    TEST(DenseBQMFunctionTest, from_ising)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;

        auto bqm = BinaryQuadraticModel_Dense<uint32_t, double>::from_ising(linear, quadratic, offset);
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

    TEST(DenseBQMFunctionTest, remove_variable)
    {
        Linear<std::string, double> linear{ {"a", 0.0}, {"b", 1.0}, {"c", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), 0.25}, {std::make_pair("a", "c"), 0.5}, {std::make_pair("b", "c"), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);
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

    TEST(DenseBQMFunctionTest, remove_variables_from)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.25}, {std::make_pair(0, 2), 0.5}, {std::make_pair(1, 2), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);
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

    TEST(DenseBQMFunctionTest, remove_interaction)
    {
        Linear<std::string, double> linear{ {"c", 2} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}, {std::make_pair("b", "c"), 1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);
        bqm.remove_interaction("b", "c");

        EXPECT_EQ(bqm.contains("a"), true);
        EXPECT_EQ(bqm.contains("b"), true);
        EXPECT_EQ(bqm.contains("c"), true);

        EXPECT_EQ(bqm.get_linear().size(), 1);
        EXPECT_EQ(bqm.get_quadratic().size(), 1);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);

        bqm.remove_interaction("a", "b");

        EXPECT_EQ(bqm.contains("a"), false);
        EXPECT_EQ(bqm.contains("b"), false);
        EXPECT_EQ(bqm.contains("c"), true);

        EXPECT_EQ(bqm.get_linear().size(), 1);
        EXPECT_EQ(bqm.get_quadratic().size(), 0);

        EXPECT_EQ(bqm.interaction_matrix().rows(), bqm.get_num_variables() + 1);
        EXPECT_EQ(bqm.interaction_matrix().cols(), bqm.get_num_variables() + 1);
    }

    TEST(DenseBQMFunctionTest, scale)
    {
        Linear<std::string, double> linear{{"a", -2.0}, {"b", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}};
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);

        bqm.scale(0.5);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();
        double bqm_offset = bqm.get_offset();

        // check biases and offset
        EXPECT_DOUBLE_EQ(bqm_linear["a"], -1.0);
        EXPECT_DOUBLE_EQ(bqm_quadratic[std::make_pair("a", "b")], -0.5);
        EXPECT_DOUBLE_EQ(bqm_offset, 0.5);
    }

    TEST(DenseBQMFunctionTest, normalize)
    {
        Linear<std::string, double> linear{ {"a", -2.0}, {"b", 1.5} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, fix_variable)
    {
        Linear<std::string, double> linear{ {"a", -0.5}, {"b", 0.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);
        
        bqm.fix_variable("a", -1);
        auto bqm_linear = bqm.get_linear();
        double bqm_offset = bqm.get_offset();

        // check offset, linear bias and variable
        EXPECT_DOUBLE_EQ(bqm_offset, 0.5);
        EXPECT_DOUBLE_EQ(bqm_linear["b"], 1.0);
        EXPECT_EQ(bqm.contains("a"), false);
    }

    TEST(DenseBQMFunctionTest, flip_variable)
    {
        Linear<uint32_t, double> linear{{1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 0.5} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);
        
        bqm.flip_variable(1);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        // check biases
        EXPECT_DOUBLE_EQ(bqm_linear[1], -1.0);
        EXPECT_DOUBLE_EQ(bqm_linear[2], 2.0);
        EXPECT_DOUBLE_EQ(bqm_quadratic[std::make_pair(1, 2)], -0.5);
    }

    TEST(DenseBQMFunctionTest, flip_variable_binary)
    {
        Linear<uint32_t, double> linear{{1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 0.5} };
        double offset = 0.5;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);
        
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
    //TEST(DenseBQMFunctionTest, contract_variables)
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

    //    BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);
    //    bqm.contract_variables(2, 3);

    //    auto bqm_quadratic = bqm.get_quadratic();

    //    // check variable and quadratic bias
    //    EXPECT_EQ(bqm.contains(3), false);
    //    EXPECT_EQ(bqm_quadratic[std::make_pair(1, 2)], 25.0);
    //}
    
    TEST(DenseBQMFunctionTest, change_vartype)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel_Dense<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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

    TEST(DenseBQMFunctionTest, to_serializable)
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 5.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel_Dense<std::string, double> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();
    }

    TEST(DenseBQMFunctionTest, from_serializable)
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 5.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();
        //std::cout << j << std::endl;

        BinaryQuadraticModel<std::string, double> bqm2 = BinaryQuadraticModel<std::string, double>::from_serializable(j);
        //bqm2.print();

        auto bqm_linear = bqm2.get_linear();
        auto bqm_quadratic = bqm2.get_quadratic();
        EXPECT_EQ(bqm2.get_offset(), bqm.get_offset());
        EXPECT_EQ(bqm_linear["c"], linear["c"]);
        EXPECT_EQ(bqm_quadratic[std::make_pair("b", "e")], quadratic[std::make_pair("b", "e")]);
    }










    TEST(ConstructionTest, Construction)
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

        BinaryQuadraticModel<uint32_t, double> bqm_k4(linear, quadratic, offset, vartype);

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
    }

    TEST(ConstructionTest, ConstructionString)
    {
        Linear<std::string, double> linear{ {"a", 1.0}, {"b", 2.0}, {"c", 3.0}, {"d", 4.0} };
        Quadratic<std::string, double> quadratic
        {
            {std::make_pair("a", "b"), 12.0}, {std::make_pair("a", "c"), 13.0}, {std::make_pair("a", "d"), 14.0},
            {std::make_pair("b", "c"), 23.0}, {std::make_pair("b", "d"), 24.0},
            {std::make_pair("d", "c"), 34.0}
        };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel<std::string, double> bqm_k4(linear, quadratic, offset, vartype, "BQM_Binary");

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

    TEST(FunctionTest, add_variable)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        // check length
        EXPECT_EQ(bqm.length(), 2);

        bqm.add_variable(2, 2.0, Vartype::SPIN);
        bqm.add_variable(1, 0.33, Vartype::SPIN);
        bqm.add_variable(0, 0.33, Vartype::BINARY);

        Linear<uint32_t, double> bqm_linear = bqm.get_linear();

        // check length
        EXPECT_EQ(bqm.length(), 3);
        // check linear
        EXPECT_EQ(bqm_linear[1], 1.33);
    }

    TEST(FunctionTest, add_variables_from)
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        // check length
        EXPECT_EQ(bqm.length(), 0);

        Linear<uint32_t, double> linear2 = { {0, 0.5}, {1, -1.0} };
        Linear<uint32_t, double> linear3 = { {1, -1.0}, {2, 2.0} };

        bqm.add_variables_from(linear2, Vartype::SPIN);
        // check variable
        EXPECT_EQ(bqm.contains(1), true);

        bqm.add_variables_from(linear3);

        // check bias
        Linear<uint32_t, double> bqm_linear = bqm.get_linear();
        EXPECT_EQ(bqm_linear[1], -2.0);
    }

    TEST(FunctionTest, add_interaction)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        //bqm.print();

        bqm.add_interaction(0, 2, 2);
        bqm.add_interaction(0, 1, 0.25);
        bqm.add_interaction(1, 2, 0.25, Vartype::BINARY);

        // check bias
        Quadratic<uint32_t, double> bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], 0.75);

        //bqm.print();
    }

    TEST(FunctionTest, add_interactions_from)
    {
        Linear<uint32_t, double> linear;
        Quadratic<uint32_t, double> quadratic;
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        Quadratic<uint32_t, double> quadratic1{ {std::make_pair(0, 1), -0.5} };
        bqm.add_interactions_from(quadratic1);

        // check bias
        Quadratic<uint32_t, double> bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], -0.5);

        Quadratic<uint32_t, double> quadratic2{ {std::make_pair(0, 1), -0.5}, {std::make_pair(0, 2), 2.0} };
        Quadratic<uint32_t, double> quadratic3{ {std::make_pair(1, 2), 2.0} };
        bqm.add_interactions_from(quadratic2);
        bqm.add_interactions_from(quadratic3, Vartype::BINARY);

        // check length
        EXPECT_EQ(bqm.length(), 3);

        // check bias
        bqm_quadratic = bqm.get_quadratic();
        EXPECT_EQ(bqm_quadratic[std::make_pair(0, 1)], -1.0);
    }

    TEST(FunctionTest, add_offset)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        bqm.add_offset(1.0);
        
        // check offset
        EXPECT_EQ(bqm.get_offset(), 0.5);
    }

    TEST(FunctionTest, energy)
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        // check energy
        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        EXPECT_EQ(bqm.energy(sample1), -0.5);
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        EXPECT_EQ(bqm.energy(sample2), 3.5);
    }

    TEST(FunctionTest, energies)
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
        Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
        std::vector<Sample<uint32_t>> samples{sample1, sample2};
        std::vector<double> en_vec = bqm.energies(samples);

        // check energies
        EXPECT_EQ(en_vec[0], -0.5);
        EXPECT_EQ(en_vec[1], 3.5);
    }

    TEST(FunctionTest, to_qubo)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        auto t_qubo = bqm.to_qubo();
        Quadratic<uint32_t, double> Q = std::get<0>(t_qubo);
        double offset_qubo = std::get<1>(t_qubo);

        // check quadratic matrix and offset
        EXPECT_EQ(Q[std::make_pair(0, 0)], 1.0);
        EXPECT_EQ(Q[std::make_pair(0, 1)], 2.0);
        EXPECT_EQ(Q[std::make_pair(1, 1)], -6.0);
        EXPECT_EQ(Q[std::make_pair(1, 2)], 6.0);
        EXPECT_EQ(Q[std::make_pair(2, 2)], -2.0);
        EXPECT_EQ(offset_qubo, 2.9);
    }

    TEST(FunctionTest, to_ising)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        auto t_ising = bqm.to_ising();
        Linear<uint32_t, double> h = std::get<0>(t_ising);
        Quadratic<uint32_t, double> J = std::get<1>(t_ising);
        double offset_ising = std::get<2>(t_ising);

        // check biases
        EXPECT_EQ(h[0], 1.0);
        EXPECT_EQ(J[std::make_pair(0, 1)], 0.5);
        EXPECT_EQ(h[1], -1.0);
        EXPECT_EQ(J[std::make_pair(1, 2)], 1.5);
        EXPECT_EQ(h[2], 0.5);
        EXPECT_EQ(offset_ising, 1.4);
    }

    TEST(FunctionTest, from_qubo)
    {
        Quadratic<uint32_t, double> origQ{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5}, {std::make_pair(0, 0), 1.0}, {std::make_pair(1, 1), -1.0}, {std::make_pair(2, 2), 0.5}};
        double offset = 1.4;

        auto bqm = BinaryQuadraticModel<uint32_t, double>::from_qubo(origQ, offset);

        Linear<uint32_t, double> linear = bqm.get_linear();
        Quadratic<uint32_t, double> Q = bqm.get_quadratic();
        double offset_ising = bqm.get_offset();

        // check quadratic matrix and offset
        EXPECT_NEAR(linear[0], 1.0, 1e-5);
        EXPECT_NEAR(Q[std::make_pair(0, 1)], 0.5, 1e-5);
        EXPECT_NEAR(linear[1], -1.0, 1e-5);
        EXPECT_NEAR(Q[std::make_pair(1, 2)], 1.5, 1e-5);
        EXPECT_NEAR(linear[2], 0.5, 1e-5);
        EXPECT_NEAR(offset_ising, 1.4, 1e-5);
    }

    TEST(FunctionTest, from_ising)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;

        auto bqm = BinaryQuadraticModel<uint32_t, double>::from_ising(linear, quadratic, offset);
        Linear<uint32_t, double> h = bqm.get_linear();
        Quadratic<uint32_t, double> J = bqm.get_quadratic();
        double offset_ising = bqm.get_offset();

        // check biases
        EXPECT_EQ(h[0], 1.0);
        EXPECT_EQ(J[std::make_pair(0, 1)], 0.5);
        EXPECT_EQ(h[1], -1.0);
        EXPECT_EQ(J[std::make_pair(1, 2)], 1.5);
        EXPECT_EQ(h[2], 0.5);
        EXPECT_EQ(offset_ising, 1.4);
    }

    TEST(FunctionTest, remove_interaction)
    {
        Linear<std::string, double> linear;
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}, {std::make_pair("b", "c"), 1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);
        //bqm.print();
        bqm.remove_interaction("b", "c");
        //bqm.print();
    }

    TEST(FunctionTest, remove_variable)
    {
        Linear<std::string, double> linear{ {"a", 0.0}, {"b", 1.0}, {"c", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), 0.25}, {std::make_pair("a", "c"), 0.5}, {std::make_pair("b", "c"), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);
        //bqm.print();
        bqm.remove_variable("a");
        //bqm.print();

        // check variables
        EXPECT_EQ(bqm.contains("a"), false);
        EXPECT_EQ(bqm.contains("b"), true);
        EXPECT_EQ(bqm.contains("c"), true);
    }

    TEST(FunctionTest, remove_variables_from)
    {
        Linear<uint32_t, double> linear{ {0, 0.0}, {1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.25}, {std::make_pair(0, 2), 0.5}, {std::make_pair(1, 2), 0.75} };
        double offset = -0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
        //bqm.print();

        std::vector<uint32_t> variables = {0, 1};

        bqm.remove_variables_from(variables);
        //bqm.print();

        // check variables
        EXPECT_EQ(bqm.contains(0), false);
        EXPECT_EQ(bqm.contains(1), false);
        EXPECT_EQ(bqm.contains(2), true);
    }

    TEST(FunctionTest, scale)
    {
        Linear<std::string, double> linear{{"a", -2.0}, {"b", 2.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0}};
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);

        bqm.scale(0.5);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();
        double bqm_offset = bqm.get_offset();

        // check biases and offset
        EXPECT_EQ(bqm_linear["a"], -1.0);
        EXPECT_EQ(bqm_quadratic[std::make_pair("a", "b")], -0.5);
        EXPECT_EQ(bqm_offset, 0.5);
    }

    TEST(FunctionTest, normalize)
    {
        Linear<std::string, double> linear{ {"a", -2.0}, {"b", 1.5} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 1.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);

        bqm.normalize(std::make_pair(-1.0, 1.0));
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        auto comp = [](const auto &a, const auto &b) { return std::abs(a.second) < std::abs(b.second); };
        auto lin_max = std::max_element(bqm_linear.begin(), bqm_linear.end(), comp);
        auto quad_max = std::max_element(bqm_quadratic.begin(), bqm_quadratic.end(), comp);

        // check maximum biases
        EXPECT_EQ(lin_max->second, -1.0);
        EXPECT_EQ(quad_max->second, -0.5);
    }

    TEST(FunctionTest, fix_variable)
    {
        Linear<std::string, double> linear{ {"a", -0.5}, {"b", 0.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "b"), -1.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);
        
        bqm.fix_variable("a", -1);
        auto bqm_linear = bqm.get_linear();
        double bqm_offset = bqm.get_offset();

        // check offset, linear bias and variable
        EXPECT_EQ(bqm_offset, 0.5);
        EXPECT_EQ(bqm_linear["b"], 1.0);
        EXPECT_EQ(bqm.contains("a"), false);
    }

    TEST(FunctionTest, flip_variable)
    {
        Linear<uint32_t, double> linear{{1, 1.0}, {2, 2.0} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 0.5} };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
        
        bqm.flip_variable(1);
        auto bqm_linear = bqm.get_linear();
        auto bqm_quadratic = bqm.get_quadratic();

        // check biases
        EXPECT_EQ(bqm_linear[1], -1.0);
        EXPECT_EQ(bqm_linear[2], 2.0);
        EXPECT_EQ(bqm_quadratic[std::make_pair(1, 2)], -0.5);
    }
    TEST(FunctionTest, contract_variables)
    {
        Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
        Quadratic<uint32_t, double> quadratic
        {
            {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
            {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
            {std::make_pair(3, 4), 34.0}
        };
        double offset = 0.5;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
        bqm.contract_variables(2, 3);

        auto bqm_quadratic = bqm.get_quadratic();

        // check variable and quadratic bias
        EXPECT_EQ(bqm.contains(3), false);
        EXPECT_EQ(bqm_quadratic[std::make_pair(1, 2)], 25.0);
    }
    TEST(FunctionTest, change_vartype)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

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
    TEST(FunctionTest, interaction_matrix)
    {
        Linear<uint32_t, double> linear{{0, 1.0}, {1, -1.0}, {2, 0.5} };
        Quadratic<uint32_t, double> quadratic{ {std::make_pair(0, 1), 0.5}, {std::make_pair(1, 0), 0.5}, {std::make_pair(1, 2), 1.5} };
        double offset = 1.4;
        Vartype vartype = Vartype::SPIN;

        BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

        std::vector<uint32_t> indices = {0,1,2};

        BinaryQuadraticModel<uint32_t, double>::Matrix mat = bqm.interaction_matrix(indices);

        // check quadratic matrix and offset
        EXPECT_EQ(mat(0,0),  1);
        EXPECT_EQ(mat(0,1),  1);
        EXPECT_EQ(mat(0,2),  0);
        EXPECT_EQ(mat(1,0),  1);
        EXPECT_EQ(mat(1,1), -1);
        EXPECT_EQ(mat(1,2),  1.5);
        EXPECT_EQ(mat(2,0),  0);
        EXPECT_EQ(mat(2,1),  1.5);
        EXPECT_EQ(mat(2,2),  0.5);

    }
    TEST(FunctionTest, to_serializable)
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();
        //std::cout << j << std::endl;
    }

    TEST(FunctionTest, from_serializable)
    {
        Linear<std::string, double> linear{ {"c", -1.0}, {"d", 1.0} };
        Quadratic<std::string, double> quadratic{ {std::make_pair("a", "d"), 2.0}, {std::make_pair("b", "e"), 5.0}, {std::make_pair("a", "c"), 3.0} };
        double offset = 0.0;
        Vartype vartype = Vartype::BINARY;

        BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);

        json j = bqm.to_serializable();
        //std::cout << j << std::endl;

        BinaryQuadraticModel<std::string, double> bqm2 = BinaryQuadraticModel<std::string, double>::from_serializable(j);
        //bqm2.print();

        auto bqm_linear = bqm2.get_linear();
        auto bqm_quadratic = bqm2.get_quadratic();
        EXPECT_EQ(bqm_linear["c"], linear["c"]);
        EXPECT_EQ(bqm_quadratic[std::make_pair("b", "e")], quadratic[std::make_pair("b", "e")]);
    }

TEST(FunctionTest, empty) {
   Linear<uint32_t, double> linear{ {1, 1.0}, {2, 1.0} };
   Quadratic<uint32_t, double> quadratic{ {std::make_pair(1, 2), 1.0} };
   double offset = 0.5;
   Vartype vartype = Vartype::SPIN;

   BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

   Sample<uint32_t> sample1{ {1, -1}, {2, -1} };
   Sample<uint32_t> sample2{ {1, 1}, {2, 1} };
   
   bqm.empty();
   //Chech if the methods in Binary Polynomial Model work properly after executing empty()
   EXPECT_EQ(bqm.length(), 0);
   EXPECT_TRUE(bqm._generate_indices().empty());
   bqm.remove_offset();
   bqm.remove_variable(1);
   bqm.remove_variables_from(std::vector<uint32_t>{1,2});
   bqm.remove_interaction(1, 2);
   bqm.remove_interactions_from(std::vector<std::pair<uint32_t, uint32_t>>{std::make_pair(1, 2), std::make_pair(1, 3)});
   bqm.scale(1.0);
   bqm.normalize();
   bqm.fix_variable(1, 1);
   bqm.change_vartype(Vartype::SPIN);
   bqm.change_vartype(Vartype::BINARY);
   bqm.change_vartype(Vartype::NONE);
   bqm.to_qubo();
   bqm.to_ising();

   EXPECT_EQ(bqm.get_vartype(), Vartype::NONE);
   EXPECT_EQ(bqm.get_info(), "");
   EXPECT_DOUBLE_EQ(bqm.get_offset(), 0.0);
   EXPECT_TRUE(bqm.get_linear().empty());
   EXPECT_TRUE(bqm.get_quadratic().empty());
   EXPECT_TRUE(bqm.get_adjacency().empty());
   
   //energy
   EXPECT_DOUBLE_EQ(bqm.energy(sample1), 0.0);
   EXPECT_DOUBLE_EQ(bqm.energy(sample2), 0.0);
   
   //Reset quadratic model
   bqm.add_interaction(1, 2, 1.0, Vartype::SPIN);
   bqm.add_variable(1, 1.0);
   bqm.add_variable(2, 1.0);
   bqm.add_offset(0.5);
   
   //energy
   EXPECT_DOUBLE_EQ(bqm.energy(sample1), -0.5);
   EXPECT_DOUBLE_EQ(bqm.energy(sample2), 3.5);
   
   // check linear
   for(auto &it : bqm.get_linear()) {
      EXPECT_EQ(it.second, linear[it.first]);
   }
   
   // check quadratic
   for(auto &it : bqm.get_quadratic()) {
      EXPECT_DOUBLE_EQ(it.second, quadratic[it.first]);
   }
   // check offset
   EXPECT_DOUBLE_EQ(offset, bqm.get_offset());
   // check vartype
   EXPECT_EQ(vartype, bqm.get_vartype());
   
   //check adj
   EXPECT_DOUBLE_EQ(bqm.get_adjacency().at(1).at(2), 1.0);
   
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
