#include "gtest/gtest.h"

#include "../src/binary_quadratic_model.hpp"
#include "../src/binary_polynomial_model.hpp"
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
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1}         ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2}         ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{3}         ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{4}         ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 2}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 3}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 3}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{3, 4}      ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 2, 3}   ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 2, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 3, 4}   ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{1, 2, 3, 4}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.get_polynomial_value()));

   //sorted_variables
   auto sorted_variables = bpm._generate_sorted_variables();
   EXPECT_EQ(sorted_variables[0], 1);
   EXPECT_EQ(sorted_variables[1], 2);
   EXPECT_EQ(sorted_variables[2], 3);
   EXPECT_EQ(sorted_variables[3], 4);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   EXPECT_EQ(bpm.get_variables().count(3), 1);
   EXPECT_EQ(bpm.get_variables().count(4), 1);
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
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-1}            ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-2}            ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-3}            ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4}            ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-2, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-3, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -1}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-3, -2}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -2}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -3}        ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-3, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -2, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -3, -1}    ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -3, -2}    ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<int32_t>{-4, -3, -2, -1}), 1);

   //Polynomial Val
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.get_polynomial_value()));

   //sorted_variables
   auto sorted_variables = bpm._generate_sorted_variables();
   EXPECT_EQ(sorted_variables[0], -4);
   EXPECT_EQ(sorted_variables[1], -3);
   EXPECT_EQ(sorted_variables[2], -2);
   EXPECT_EQ(sorted_variables[3], -1);
   
   //variables
   EXPECT_EQ(bpm.get_variables().count(-1), 1);
   EXPECT_EQ(bpm.get_variables().count(-2), 1);
   EXPECT_EQ(bpm.get_variables().count(-3), 1);
   EXPECT_EQ(bpm.get_variables().count(-4), 1);
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
   
   //Polynomial Key
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a"}               ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"b"}               ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"c"}               ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"d"}               ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "b"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"b", "c"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"b", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"c", "d"}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "b", "c"}     ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "b", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"b", "c", "d"}     ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<std::string>{"a", "b", "c", "d"}), 1);

   //Polynomial Value
   EXPECT_TRUE(EXPECT_CONTAIN(1.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(2.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(3.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(4.0   , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(12.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(13.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(14.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(23.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(24.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(34.0  , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(123.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(124.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(134.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(234.0 , bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(1234.0, bpm.get_polynomial_value()));

   //sorted_variables
   auto sorted_variables = bpm._generate_sorted_variables();
   EXPECT_EQ(sorted_variables[0], "a");
   EXPECT_EQ(sorted_variables[1], "b");
   EXPECT_EQ(sorted_variables[2], "c");
   EXPECT_EQ(sorted_variables[3], "d");
   
   //variables
   EXPECT_EQ(bpm.get_variables().count("a"), 1);
   EXPECT_EQ(bpm.get_variables().count("b"), 1);
   EXPECT_EQ(bpm.get_variables().count("c"), 1);
   EXPECT_EQ(bpm.get_variables().count("d"), 1);
}

template<typename IndexType>
void StateTestBPMEmpty(const BinaryPolynomialModel<IndexType, double> &bpm) {
   
   EXPECT_EQ(bpm.get_num_variables(), 0);

   EXPECT_DOUBLE_EQ(bpm.get_offset(), 0.0);
   
   EXPECT_EQ(bpm.get_num_interactions(), 0);
   
   EXPECT_EQ(bpm.get_degree(), 0);
      
   //Polynomial
   EXPECT_EQ(bpm.get_polynomial().size(), 0);
   
   //Polynomial Key
   EXPECT_EQ(bpm.get_polynomial_key().size(), 0);
   
   //Polynomial Val
   EXPECT_EQ(bpm.get_polynomial_value().size(), 0);

   //sorted_variables
   auto sorted_variables = bpm._generate_sorted_variables();
   EXPECT_EQ(sorted_variables.size(), 0);
   
   //variables
   EXPECT_EQ(bpm.get_variables().size(), 0);

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
   
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{}          ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{7}         ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 11}     ), 1);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{11, 12, 14}), 1);
   
   EXPECT_TRUE(EXPECT_CONTAIN(-3.0, bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(-2.0, bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(-9.0, bpm.get_polynomial_value()));
   EXPECT_TRUE(EXPECT_CONTAIN(-1.0, bpm.get_polynomial_value()));

   bpm.remove_interaction({}          );
   bpm.remove_interaction({7}         );
   bpm.remove_interaction({2, 11}     );
   bpm.remove_interaction({11, 12, 14});
   
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({}          ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({7}         ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({2, 11}     ), 0.0);
   EXPECT_DOUBLE_EQ(bpm.get_polynomial({11, 12, 14}), 0.0);
   
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm.get_polynomial_value()));
   
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
   
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{}          ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{7}         ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{2, 11}     ), 0);
   EXPECT_EQ(std::count(bpm.get_polynomial_key().begin(), bpm.get_polynomial_key().end(), std::vector<uint32_t>{11, 12, 14}), 0);
   
   EXPECT_FALSE(EXPECT_CONTAIN(-3.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-2.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-9.0, bpm.get_polynomial_value()));
   EXPECT_FALSE(EXPECT_CONTAIN(-1.0, bpm.get_polynomial_value()));
   
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
   EXPECT_TRUE(bpm.get_variables().empty());
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
   
   EXPECT_EQ(bpm.get_variables().count(0), 1);
   EXPECT_EQ(bpm.get_variables().count(1), 1);
   EXPECT_EQ(bpm.get_variables().count(2), 1);
   
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

