//    Copyright 2020 Jij Inc.

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

/**
 * @mainpage cimod
 *
 * @section s_overview Overview
 * cimod is a C++ library for a binary quadratic model.
 * This library provides a binary quadratic model class which contains an Ising model or a quadratic unconstrained binary optimization (QUBO) model.
 * It also provides utilities for constructing a model and transforming to some other interfaces.
 * This library is created based on dimod (https://github.com/dwavesystems/dimod).
 *
 * @section s_bqm Binary quadratic model
 * A binary quadratic model class can contain an Ising model or a QUBO model.
 *
 * @subsection ss_ising Ising model
 * An energy of an Ising model \f$E_{\mathrm{Ising}}\f$ is represented by
 * \f[
 * E_{\mathrm{Ising}} = \sum_{i} h_i s_i + \sum_{i \neq j} J_{ij} s_i s_j + \delta_{\mathrm{Ising}},
 * \f]
 * where \f$s_i \in \{+1, -1\}\f$ denotes a spin at the site \f$i\f$, \f$h_i\f$ denotes an external magnetic field parameter, \f$J_{ij}\f$ denotes an interaction parameter and \f$\delta_{\mathrm{Ising}}\f$ denotes an offset.
 * Note that this library assumes that the interaction is not symmetric, i.e., \f$J_{ij} \neq J_{ji}\f$.
 *
 * @subsection ss_qubo QUBO model
 * An evaluation function of a QUBO model \f$E_{\mathrm{QUBO}}\f$ is represented by
 * \f[
 * E_{\mathrm{QUBO}} = \sum_{i, j} Q_{ij} x_i x_j + \delta_{\mathrm{QUBO}},
 * \f]
 * where \f$x_i \in \{0, 1\}\f$ denotes a decision variable, \f$Q_{ij}\f$ denotes a quadratic bias and \f$\delta_{\mathrm{QUBO}}\f$ denotes an offset.
 * Note that this library assumes that the quadratic bias is not symmetric, i.e., \f$Q_{ij} \neq Q_{ji}\f$ if \f$i \neq j\f$.
 *
 * @section s_example Example
 * @code
 * #include "src/binary_quadratic_model.hpp"
 *
 * using namespace cimod;
 * int main()
 * {
 * // Set linear biases and quadratic biases
 * Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
 * Quadratic<uint32_t, double> quadratic
 * {
 *      {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
 *      {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
 *      {std::make_pair(3, 4), 34.0}
 *  };
 *
 * // Set offset
 * double offset = 0.0;
 *
 * // Set variable type
 * Vartype vartype = Vartype::BINARY;
 * // Create a BinaryQuadraticModel instance
 * BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
 *
 * // Print informations of bqm
 * bqm.print();
 *
 * return 0;
 * }
 * @endcode
 */

/**
 * @file binary_quadratic_model.hpp
 * @author Fumiya Watanabe
 * @brief BinaryQuadraticModel class
 * @version 1.0.0
 * @date 2020-03-24
 *
 * @copyright Copyright (c) Jij Inc. 2020
 *
 */

#ifndef binary_polynomial_model_hpp
#define binary_polynomial_model_hpp

#include "vartypes.hpp"
#include "hash.hpp"
#include "utilities.hpp"
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cimod {

//! @brief Type alias for variable list
//! @tparam IndexType
template <typename IndexType>
using Variable = std::unordered_set<IndexType>;

//! @brief Type alias for polynomial bias
//! @tparam IndexType
//! @tparam FloatType
template <typename IndexType, typename FloatType>
using Polynomial = std::unordered_map<std::vector<IndexType>, FloatType, vector_hash>;

//! @brief Type alias for adjacency list
//! @tparam IndexType
//! @tparam FloatType
template <typename IndexType, typename FloatType>
using Adjacency_Poly = std::unordered_map<IndexType, Polynomial<IndexType, FloatType>>;

//! @brief Type alias for sample
//! @tparam IndexType
template <typename IndexType>
using Sample = std::unordered_map<IndexType, int32_t>;

//! @brief Class for binary polynomial model
//! @tparam IndexType
//! @tparam FloatType
template <typename IndexType, typename FloatType>
class BinaryPolynomialModel {
   
public:
      
   //! @brief BinaryPolynomialModel constructor
   //! @param polynomial
   //! @param vartype
   //! @param info
   BinaryPolynomialModel(const Polynomial<IndexType, FloatType> &polynomial,
                         const Vartype vartype,
                         const std::string info = ""): m_vartype(vartype), m_info(info) {
      add_interactions_from(polynomial);
   };
   
   //! @brief Copy constructor of BinaryPolynomialModel
   //! @param BinaryPolynomialModel
   BinaryPolynomialModel(const BinaryPolynomialModel&) = default;
   
   //! @brief Move constructor of BinaryPolynomialModel
   //! @param BinaryPolynomialModel
   BinaryPolynomialModel(BinaryPolynomialModel&&) = default;
   
   //! @brief Generate variable list associated with the input interactions
   //! @return Sorted variable list as std::vector<IndexType>
   std::vector<IndexType> generate_variables() const {
      std::vector<IndexType> ret;
      for (auto &it_variables: m_variables) {
         ret.push_back(it_variables);
      }
       std::sort(ret.begin(), ret.end());
       return ret;
   }
   
   //! @brief Return the number of variables.
   //! @return The number of variables.
   size_t length() const {
      return m_variables.size();
   }
   
   //! @brief Check if the variable v is contained in the variable list
   //! @return True if the variable list contain v, otherwise false
   bool contains(const IndexType &v) {
      if (m_variables.count(v) != 0) {
         return true;
      }
      else {
         return false;
      }
   }
   
   //! @brief Get the variable object
   //! @return Variable list
   const Variable<IndexType> &get_variables() const  {
      return this->m_variables;
   }
   
   //! @brief Get the Polynomial object
   //! @return Polynomial bias
   const Polynomial<IndexType, FloatType> &get_polynomial() const {
      return this->m_polynomial;
   }
   
   //! @brief Get the Adjacency object
   //! @return Adjacency list
   const Adjacency_Poly<IndexType, FloatType> &get_adjacency() const {
      return this->m_adj;
   }
   
   //! @brief Get the vartype object, which represents the type of the model
   //! @return Vartype
   const Vartype &get_vartype() const {
      return this->m_vartype;
   }
   
   //! @brief Get the info object
   //! @return Information
   const std::string &get_info() const {
      return this->m_info;
   }
   
   //! @brief Print information of binary polynomial model
   void print() const {
      
      std::vector<IndexType> variables = generate_variables();
      
      std::cout << "[BinaryPolynomialModel]" << std::endl;
      
      std::cout << "Variables = " << std::endl;
      for (auto &it_variables: variables) {
         std::cout << it_variables << ", ";
      }
      std::cout << std::endl;
      
      // Print linear, which is stored in m_polynomial with the first index std::vector of the size = 1
      std::cout << "polynomial(linear) = " << std::endl;
      for (auto &it_variables: variables) {
         std::cout << it_variables << ": " << m_polynomial.at(std::vector<IndexType>{it_variables}) << std::endl;
      }
      
      // Print polynomial, which is stored in m_polynomial with the first index std::vector of the size > 1
      std::cout << "polynomial(interaction) = " << std::endl;
      for(auto &it_polynomial : m_polynomial) {
         if (it_polynomial.first.size() > 1) {
            std::cout << "(";
            for (auto &it_index : it_polynomial.first) {
               std::cout << it_index << ", ";
            }
            std::cout << "): " << it_polynomial.second << ", ";
         }
      }
      std::cout << std::endl;
      
      // Print adjacency
      std::cout << "adjacency = " << std::endl;
      for (auto &it_variables: variables) {
         std::cout << it_variables << ": {";
         if (m_adj.count(it_variables) > 0) {
            for (auto &it_adj: m_adj.at(it_variables)) {
               std::cout << "(";
               for (auto &it_interaction: it_adj.first) {
                  std::cout << it_interaction << ", ";
               }
               std::cout << "): " << it_adj.second << ", ";
            }
         }
         std::cout << "}" << std::endl;
      }
      
      // Print vartype
      std::cout << "vartype = ";
      if(m_vartype == Vartype::SPIN) {
         std::cout << "Spin" << std::endl;
      }
      else if(m_vartype == Vartype::BINARY) {
         std::cout << "Binary" << std::endl;
      }
      else {
         std::cout << "Unknown vartype" << std::endl;
      }
      
      // Print info
      std::cout << "info = ";
      std::cout << "\"" << m_info << "\"" << std::endl;
      
   }
   
   //! @brief Create an empty binary polynomial model
   void empty() {
      m_variables  = {};
      m_polynomial = {};
      m_vartype    = Vartype::NONE;
      m_info       = "";
   }
   
   //! @brief Add variable v and correponding linear bias to binary polynomial model
   //! @param v
   //! @param bias
   //! @param vartype
   void add_linear(const IndexType &v, const FloatType &bias, const Vartype vartype = Vartype::NONE) {
      
      if (vartype != Vartype::NONE) {
         if ((m_vartype == Vartype::SPIN) && (vartype == Vartype::BINARY)) {
            std::cerr << "Cannot convert the vartype=SPIN to vartype=BINARY" << std::endl;
         }
         else if((m_vartype == Vartype::BINARY) && (vartype == Vartype::SPIN)) {
            std::cerr << "Cannot convert the vartype=BINARY to vartype=SPIN" << std::endl;
         }
         else {
            std::cerr << "Unknown vartype" << std::endl;
         }
      }
      
      // Insert or assign the bias
      FloatType value = 0;
      std::vector<IndexType> index = {v};
      if (m_polynomial.count(index) != 0) {
         value = m_polynomial[index];
      }
      insert_or_assign(m_polynomial, index, value + bias);
   }
      
   //! @Add interaction corresponding bias to binary polynomial model
   //! @param u
   //! @param bias
   //! @param vartype
   void add_interaction(const std::vector<IndexType> &u, const FloatType &bias, const Vartype vartype = Vartype::NONE) {
      
      //Check the input interaction is valid
      for (auto &it : u) {
         if (std::count(u.begin(), u.end(), it) != 1) {
            std::cerr << "No self-loops allowed, therefore (";
            for (auto &it_print : u) {
               std::cerr << it_print << ", ";
            }
            std::cerr << ") is not an allowed interaction" << std::endl;
            return;
         }
      }
      
      if((vartype != Vartype::NONE) && (vartype != m_vartype)) {
         if((m_vartype == Vartype::SPIN) && (vartype == Vartype::BINARY)) {
            std::cerr << "Cannot convert vartype=SPIN to vartype=BINARY" << std::endl;
            return;
         }
         else if((m_vartype == Vartype::BINARY) && (vartype == Vartype::SPIN)) {
            std::cerr << "Cannot convert vartype=BINARY to vartype=SPIN" << std::endl;
            return;
         }
         else {
            std::cerr << "Unknown vartype" << std::endl;
            return;
         }
      }
      
      FloatType value = 0;
      if (m_polynomial.count(u) != 0) {
         value = m_polynomial[u];
      }
      insert_or_assign(m_polynomial, u, value + bias);
      
      if (u.size() > 1) {
         update_adjacency(u);
      }
      
      for (auto &it: u) {
         if (m_variables.count(it) == 0) {
            m_variables.emplace(it);
            add_linear(it, 0);
         }
      }
      
   };
   
   //! @Add interactions corresponding biases to binary polynomial model
   //! @param polynomial
   //! @param vartype
   void add_interactions_from(const Polynomial<IndexType, FloatType> &polynomial, const Vartype vartype = Vartype::NONE) {
      for(auto &it : polynomial) {
         add_interaction(it.first, it.second, vartype);
      }
   }
   
   //! @Remove variable v and corresponding interactions from binary polynomial model
   //! @param v
   void remove_variable(const IndexType &v) {
      std::vector<std::vector<IndexType>> interactions;
      for (auto &it_polynomial: m_polynomial) {
         if (std::count(it_polynomial.first.begin(), it_polynomial.first.end(), v) != 0) {
            interactions.push_back(it_polynomial.first);
         }
      }
      remove_interactions_from(interactions);
      m_adj.erase(v);
      m_variables.erase(v);
   }
   
   //! @Remove variables v and corresponding interactions from BinaryPolynomialModel
   //! @param variables
   void remove_variables_from(const std::vector<IndexType> &variables) {
      for(auto &it : variables) {
         remove_variable(it);
      }
   };
   
   //! @brief Remove interaction specified by variables from BinaryPolynomialModel
   //! @param variables
   void remove_interaction(const std::vector<IndexType> &variables) {
      if (m_polynomial.count(variables) != 0) {
         m_polynomial.erase(variables);
         remove_adjacency(variables);
      }
   };
   
   //! @brief Remove all the specified interactions from BinaryPolynomialModel
   //! @param variable_array
   void remove_interactions_from(const std::vector<std::vector<IndexType>> &variable_array) {
      for(auto &it : variable_array) {
         remove_interaction(it);
      }
   };

   //! @brief Multiply by the specified scalar all the biases and offset of a binary quadratic model.
   //! @param scalar
   //! @param ignored_interactions
   //! @param ignored_offset
   void scale(const FloatType &scalar,
              const std::vector<std::vector<IndexType>> &ignored_interactions = {},
              const bool ignored_offset = false
              ) {
      
      for (auto &it_polynomial: m_polynomial) {
         if (std::find(ignored_interactions.begin(), ignored_interactions.end(), it_polynomial.first) != ignored_interactions.end()
             || ignored_interactions.empty()) {
            it_polynomial.second *= scalar;
         }
      }
   }
   
   //! @brief Normalizes the biases of the binary polynomial model such that they fall in the provided range(s), and adjusts the offset appropriately.
   //! @param bias_range
   //! @param ignored_variables
   //! @param ignored_offset
   void normalize(const std::pair<FloatType, FloatType>     &bias_range = {1.0, 1.0},
                  const std::vector<std::vector<IndexType>> &ignored_variables = {},
                  const bool ignored_offset = false
                  ) {
      
      auto comp = [](const auto &a, const auto &b) { return a.second < b.second; };
      auto min = std::min_element(m_polynomial.begin(), m_polynomial.end(), comp)->second/bias_range.first;
      auto max = std::max_element(m_polynomial.begin(), m_polynomial.end(), comp)->second/bias_range.second;
      
      FloatType inv_scale = (min < max) ? max : min;
      
      if (inv_scale != 0.0) {
         scale(1.0/inv_scale, ignored_variables, ignored_offset);
         for (auto &it_polynomial: m_polynomial) {
            if (it_polynomial.first.size() > 1) {update_adjacency(it_polynomial.first); };
         }
      }
   }
   
   void update(const BinaryPolynomialModel &bpm, const bool ignore_info = true) {
      add_interactions_from(bpm.get_polynomial(), bpm.get_vartype());
      if(!ignore_info) {
         m_info = bpm.get_info();
      };
   };
   

   FloatType energy(const Sample<IndexType> &sample) {
      FloatType en = 0.0;
      for (auto &it_polynomial: m_polynomial) {
         IndexType multiple_variable = 1;
         for (auto &it_vec: it_polynomial.first) {
            if (check_vartype(sample.at(it_vec), m_vartype)) {
               multiple_variable *= sample.at(it_vec);
            }
            else {
               multiple_variable = 0;
               break;
            }
         }
         en += static_cast<FloatType>(multiple_variable)*it_polynomial.second;
      }
      return en;
   }
   
   std::vector<FloatType> energies(const std::vector<Sample<IndexType>> &samples_like) {
      std::vector<FloatType> en_vec;
      for(auto &it: samples_like) {
          en_vec.push_back(energy(it));
      }
      return en_vec;
   }
      
protected:
   
   //! @brief Variable list corresponding to the interactions
   Variable<IndexType> m_variables;
   
   //! @brief Polynomial biases as a dictionary.
   Polynomial<IndexType, FloatType> m_polynomial;
   
   //! @brief The model's type. SPIN or BINARY
   Vartype m_vartype = Vartype::NONE;
   
   //! @brief A place to store miscellaneous data about the binary Polynomial model as a whole.
   std::string m_info;
   
   //! @brief The model's interactions as nested dictionaries.
   Adjacency_Poly<IndexType, FloatType> m_adj;
   
   //! @brief Add the adjacency to the adjacency list
   void update_adjacency(const std::vector<IndexType> &u) {
      if(m_polynomial.count(u) != 0) {
         for (auto &it: u) {
            if (m_adj[it].count(u) == 0) {
               insert_or_assign(m_adj[it], u, m_polynomial[u]);
               break;
            }
         }
         //int min = *std::min_element(u.begin(), u.end());
         //insert_or_assign(m_adj[min], u, m_polynomial[u]);
      }
   }
   
   //! @brief Remove the adjacency from the adjacency list
   void remove_adjacency(const std::vector<IndexType> &interaction) {
      int min = *std::min_element(interaction.begin(), interaction.end());
      m_adj[min].erase(interaction);
   }
   
};

}


#endif /* binary_polynomial_model_hpp */
