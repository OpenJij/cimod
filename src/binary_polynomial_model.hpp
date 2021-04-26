//    Copyright 2021 Jij Inc.

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
 * @section s_bpm Binary polynomial model
 * A binary polynomial model, which can be regarded as an extended model of the binary quadratic model, can handle Ising and PUBO models.
 * @subsection ss_bpm_Ising Ising model
 * An energy of an "extended" Ising model \f$E_{\mathrm{Ising}}\f$ is represented by
 * \f[
 * E_{\mathrm{Ising}} = \sum_{i} h_i s_i + \sum_{i \neq j} J_{ij} s_i s_j +  \sum_{i \neq j \neq k} J_{ijk} s_i s_j s_k + \ldots
 * \f]
 * Here \f$s_i \in \{+1, -1\}\f$ denotes the spin at the site \f$i\f$, \f$ h_i \f$ denotes the external magnetic field at the site \f$ i \f$, and \f$J_{ijk\ldots}\f$ represents the interaction between the sites.
 * Note that \f$ i \neq j \neq k \f$ means \f$ i \neq j \f$, \f$ j \neq k \f$, and \f$ i \neq k \f$.
 * This library assumes that the interaction is not symmetric. For example, \f$J_{ij} \neq J_{ji}\f$ for \f$  i\neq j\f$, \f$J_{ijk} \neq J_{jik}\f$ for \f$ i \neq j \neq k \f$, and so on.
 *
 * @subsection ss_bpm_pubo PUBO model
 * An energy of an "extended" QUBO model \f$ E_{\mathrm{PUBO}}\f$, here we call polynomial unconstrained binary optimization (PUBO), is represented by
 * \f[
 * E_{\mathrm{PUBO}} = \sum_{i \neq j} Q_{ij} x_i x_j +  \sum_{i \neq j \neq k} Q_{ijk} x_i x_j x_k + \ldots
 * \f]
 * Here \f$ x_i \in \{0, 1\} \f$ denotes the spin at the site \f$ i \f$ and \f$Q_{ijk\ldots}\f$ represents the interaction between the sites.
 * Note that \f$ i \neq j \neq k \f$ means \f$ i \neq j \f$, \f$ j \neq k \f$, and \f$ i \neq k \f$.
 * This library assumes that the interaction is not symmetric. For example, \f$Q_{ij} \neq Q_{ji}\f$ for \f$  i\neq j\f$, \f$Q_{ijk} \neq Q_{jik}\f$ for \f$ i \neq j \neq k \f$, and so on.
 *
 * @section s_example Example
 * @code
 * #include "../src/binary_quadratic_model.hpp"
 * #include "../src/binary_polynomial_model.hpp"
 *
 * using namespace cimod;
 *
 * int main() {
 *
 *   // Set linear biases and quadratic biases
 *   Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
 *   Quadratic<uint32_t, double> quadratic {
 *        {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
 *        {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
 *        {std::make_pair(3, 4), 34.0}
 *    };
 *
 *   // Set variable type
 *   Vartype vartype = Vartype::BINARY;
 *
 *   // Set offset
 *   double offset = 0.0;
 *
 *   // Create a BinaryQuadraticModel instance
 *   BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
 *
 *   // Print informations of bqm
 *   bqm.print();
 *
 *   //Set polynomial biases
 *   Polynomial<uint32_t, double> polynomial {
 *      //Linear biases
 *      {{1}, 1.0}, {{2}, 2.0}, {{3}, 3.0},
 *      //Quadratic biases
 *      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{2, 3}, 23.0},
 *      //Polynomial bias
 *      {{1, 2, 3}, 123.0}
 *   };
 *
 *   // Create a BinaryPolynominalModel instance
 *   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
 *
 *   // Print informations of bpm
 *   bpm.print();
 *
 *   return 0;
 * }
 *
 * @endcode
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
#include <bitset>

namespace cimod {


//! @brief Class for binary polynomial model.
//! @tparam IndexType
//! @tparam FloatType
template <typename IndexType, typename FloatType>
class BinaryPolynomialModel {

   using InteractionKey = std::vector<IndexType>;

   using InteractionValue = FloatType;

   using InteractionMap = std::unordered_map<InteractionKey, InteractionValue, vector_hash>;

   using InteractionKeyList = std::vector<InteractionKey>;

   using InteractionValueList = std::vector<InteractionValue>;
   
public:
      
   //! @brief BinaryPolynomialModel constructor.
   //! @param polynomial
   //! @param vartype
   //! @param info
   BinaryPolynomialModel(const InteractionMap &interaction_map,
                         const Vartype vartype,
                         const std::string info = ""): vartype_(vartype), info_(info) {
      add_interactions_from(interaction_map);
   }

   void add_interaction(InteractionKey &key, const InteractionValue &value, const Vartype vartype = Vartype::NONE) {
      if (std::abs(value) < 0.0) {
         return;
      }
      std::sort(key.begin(), key.end());
      CheckKeySelfLoop(key);

      if (vartype_ == vartype || vartype == Vartype::NONE) {
         SetKeyAndValue(key, value);
      }
      else {
         const std::size_t original_key_size     = key.size();
         const std::size_t changed_key_list_size = IntegerPower(2, original_key_size);

         if (vartype_ == Vartype::SPIN && vartype == Vartype::BINARY) {
            for (std::size_t i = 0; i < original_key_size; ++i) {
               const auto changed_key = GenerateChangedKey(key, i);
               int sign = ((original_key_size - changed_key.size())%2 == 0) ? 1.0 : -1.0;
               SetKeyAndValue(changed_key, value*IntegerPower(2, changed_key.size())*sign);
            }
         }
         else if (vartype_ == Vartype::BINARY && vartype == Vartype::SPIN) {
            FloatType changed_value = value*(1.0/IntegerPower(2, original_key_size));
            for (std::size_t i = 0; i < original_key_size; ++i) {
               SetKeyAndValue(GenerateChangedKey(key, i), changed_value);
            }
         }
         else {
            throw std::runtime_error("Unknown vartype error");
         }
      }
   }

   void add_interactions_from(InteractionMap &interaction_map, const Vartype vartype = Vartype::NONE) {
      for (auto &&it: interaction_map) {
         add_interaction(it.first, it.second, vartype);
      }
   }




   


  
      
protected:
   
   std::unordered_set<IndexType> variable_set_;
   
   InteractionKeyList interaction_key_list_;

   InteractionValueList interaction_value_list_;

   std::unordered_map<InteractionKey, std::size_t, vector_hash> interaction_key_inv_;
   
   Vartype vartype_ = Vartype::NONE;
   
   std::string info_;

   void CheckKeySelfLoop(InteractionKey &key) {
      //key is assumed to be sorted
      for (std::size_t i = 0; i < key.size() - 1; ++i) {
         if (key[i] == key[i + 1]) {
            throw std::runtime_error("No self-loops allowed");
         }
      }
   }

   void SetKeyAndValue(const InteractionKey &key, const InteractionValue &value) {
      if (interaction_key_inv_.count(key) == 0) {
         interaction_key_inv_[key] = interaction_value_list_.size();
         interaction_key_list_.push_back(key);
         interaction_value_list_.push_back(value);
         variable_set_.insert(key.begin(), key.end());
      }
      else {
         interaction_value_list_[interaction_key_inv_[key]] += value;
      }
   }

   std::size_t IntegerPower(std::size_t base, std::size_t exponent) {
      std::size_t val = 1;
      for (std::size_t i = 0; i < exponent; ++i) {
         val *= base;
      }
      return val;
   }

   InteractionKey GenerateChangedKey(const InteractionKey &original_key, const std::size_t num_of_key) {
      const std::size_t original_key_size = original_key.size();
      std::bitset<original_key_size> bs(num_of_key);
      InteractionKey changed_key;
      for (std::size_t i = 0; i < original_key_size; ++i) {
         if (bs[i]) {
            changed_key.push_back(original_key[i]);
         }
      }
      return changed_key;
   }

   
};

}


#endif /* binary_polynomial_model_hpp */
