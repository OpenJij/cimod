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

template <typename IndexType, typename FloatType>
using Polynomial = std::unordered_map<std::vector<IndexType>, FloatType, vector_hash>;

template <typename IndexType>
using PolynomialKeyList = std::vector<std::vector<IndexType>>;

template <typename FloatType>
using PolynomialValueList = std::vector<FloatType>;

//! @brief Class for binary polynomial model.
//! @tparam IndexType
//! @tparam FloatType
template <typename IndexType, typename FloatType>
class BinaryPolynomialModel {
   
   
public:
   
   //! @brief BinaryPolynomialModel constructor.
   //! @param poly_map
   //! @param vartype
   BinaryPolynomialModel(const Polynomial<IndexType, FloatType> &poly_map,
                         const Vartype vartype): vartype_(vartype) {
      add_interactions_from(poly_map);
   }
   
   BinaryPolynomialModel(PolynomialKeyList<IndexType> &key_list,
                         const PolynomialValueList<FloatType> &value_list,
                         const Vartype vartype): vartype_(vartype) {
      add_interactions_from(key_list, value_list);
   }
   
   BinaryPolynomialModel(const PolynomialKeyList<IndexType> &key_list,
                         const PolynomialValueList<FloatType> &value_list,
                         const Vartype vartype): vartype_(vartype) {
      add_interactions_from(key_list, value_list);
   }
   
   Polynomial<IndexType, FloatType> get_polynomial() const {
      Polynomial<IndexType, FloatType> poly_map;
      for (std::size_t i = 0; i < poly_key_list_.size(); ++i) {
         poly_map[poly_key_list_[i]] = poly_value_list_[i];
      }
      return poly_map;
   }
   
   FloatType get_polynomial(std::vector<IndexType> &key) const {
      std::sort(key.begin(), key.end());
      CheckKeySelfLoop(key);
      if (poly_key_inv_.count(key) != 0) {
         return poly_value_list_[poly_key_inv_[key]];
      }
      else {
         return 0;
      }
   }
   
   FloatType get_polynomial(const std::vector<IndexType> &key) const {
      std::vector<IndexType> copied_key = key;
      return get_polynomial(copied_key);
   }
   
   const PolynomialKeyList<IndexType> &get_polynomial_key() const {
      return poly_key_list_;
   }
   
   const PolynomialValueList<FloatType> &get_polynomial_value() const {
      return poly_value_list_;
   }

   std::vector<IndexType> get_variables() const {
      std::vector<IndexType> variable_list;
      variable_list.reserve(get_num_variables());
      for (const auto &index: each_variable_num_) {
         variable_list.push_back(index.first);
      }
      std::sort(variable_list.begin(), variable_list.end());
      return variable_list;
   }
   
   FloatType get_offset() const {
      return get_polynomial(std::vector<IndexType>{});
   }
   
   Vartype get_vartype() const {
      return vartype_;
   }
   
   std::size_t get_num_interactions() const {
      return poly_key_list_.size();
   }
   
   std::size_t get_num_variables() const {
      return each_variable_num_.size();
   }
   
   BinaryPolynomialModel empty(const Vartype vartype) const {
      return BinaryPolynomialModel({}, vartype);
   }
   
   void clear() {
      each_variable_num_.clear();
      PolynomialKeyList<IndexType>().swap(poly_key_list_);
      PolynomialValueList<FloatType>().swap(poly_value_list_);
      poly_key_inv_.clear();
      info_ = "";
      degree_ = 0;
   }
   
   void remove_interaction(std::vector<IndexType> &key) {
      std::sort(key.begin(), key.end());
      if (poly_key_inv_.count(key) == 0) {
         return;
      }

      for (const auto &index: key) {
         if (each_variable_num_[index] >= 2) {
            each_variable_num_[index]--;
         }
         else if (each_variable_num_[index] == 1) {
            each_variable_num_.erase(index);
         }
      }

      std::size_t inv = poly_key_inv_[key];

      std::swap(poly_key_inv_[key], poly_key_inv_[poly_key_list_.back()]);
      poly_key_inv_.erase(key);

      std::swap(poly_key_list_[inv], poly_key_list_.back());
      poly_key_list_.pop_back();

      std::swap(poly_value_list_[inv], poly_value_list_.back());
      poly_value_list_.pop_back();

   }

   void remove_interaction(const std::vector<IndexType> &key) {
      std::vector<IndexType> copied_key = key;
      remove_interaction(copied_key);
   }

   void remove_interactions_from(PolynomialKeyList<IndexType> &key_list) {
      for (auto &&key: key_list) {
         remove_interaction(key);
      }
   }

   void remove_interactions_from(const PolynomialKeyList<IndexType> &key_list) {
      for (const auto &key: key_list) {
         remove_interaction(key);
      }
   }

   void remove_offset() {
      remove_interaction(FloatType{});
   }

   void remove_variable(const IndexType &index) {
      for (auto &&key: poly_key_list_) {
         if (std::binary_search(key.begin(), key.end(), index)) {
            remove_interaction(key);
         }
      }
   }

   void remove_variables_from(const std::vector<IndexType> &key) {
      for (const auto &index: key) {
         remove_variable(index);
      }
   }
   
   void add_interaction(std::vector<IndexType> &key, const FloatType &value, const Vartype vartype = Vartype::NONE) {
      if (std::abs(value) < 0.0) {
         return;
      }
      std::sort(key.begin(), key.end());
      CheckKeySelfLoop(key);
      UpdateDegree(key.size());
      if (vartype_ == vartype || vartype == Vartype::NONE) {
         SetKeyAndValue(key, value);
      }
      else {
         const std::size_t original_key_size     = key.size();
         const std::size_t changed_key_list_size = IntegerPower(2, original_key_size);
         
         if (vartype_ == Vartype::SPIN && vartype == Vartype::BINARY) {
            for (std::size_t i = 0; i < changed_key_list_size; ++i) {
               const auto changed_key = GenerateChangedKey(key, i);
               int sign = ((original_key_size - changed_key.size())%2 == 0) ? 1.0 : -1.0;
               SetKeyAndValue(changed_key, value*IntegerPower(2, changed_key.size())*sign);
            }
         }
         else if (vartype_ == Vartype::BINARY && vartype == Vartype::SPIN) {
            FloatType changed_value = value*(1.0/changed_key_list_size);
            for (std::size_t i = 0; i < changed_key_list_size; ++i) {
               SetKeyAndValue(GenerateChangedKey(key, i), changed_value);
            }
         }
         else {
            throw std::runtime_error("Unknown vartype error");
         }
      }
   }
   
   void add_interaction(const std::vector<IndexType> &key, const FloatType &value, const Vartype vartype = Vartype::NONE) {
      std::vector<IndexType> copied_key = key;
      add_interaction(copied_key, value, vartype);
   }
   
   void add_interactions_from(const Polynomial<IndexType, FloatType> &interaction_map, const Vartype vartype = Vartype::NONE) {
      for (const auto &it: interaction_map) {
         std::vector<IndexType> copied_key = it.first;
         add_interaction(copied_key, it.second, vartype);
      }
   }
   
   void add_interactions_from(PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list, const Vartype vartype = Vartype::NONE) {
      if (key_list.size() != value_list.size()) {
         throw std::runtime_error("The sizes of key_list and value_list must match each other");
      }
      for (std::size_t i = 0; i < key_list.size(); ++i) {
         add_interaction(key_list[i], value_list[i], vartype);
      }
   }
   
   void add_interactions_from(const PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list, const Vartype vartype = Vartype::NONE) {
      PolynomialKeyList<IndexType> copied_key_list = key_list;
      add_interactions_from(copied_key_list, value_list, vartype);
   }
   
   void add_offset(FloatType offset) {
      add_interaction(std::vector<IndexType>{}, offset);
   }
   
   FloatType energy(const std::vector<int32_t> &sample, bool omp_flag = true) const {
      if (sample.size() != get_num_variables()) {
         throw std::runtime_error("The size of sample must be equal to num_variables");
      }
      std::size_t num_interactions = get_num_interactions();
      FloatType val = 0.0;
      
      if (omp_flag) {
#pragma omp parallel for reduction (+: val)
         for (std::size_t i = 0; i < num_interactions; ++i) {
            int32_t spin_multiple = 1;
            for (const auto &index: poly_key_list_[i]) {
               spin_multiple *= sample[index];
               if (std::abs(spin_multiple) <= 0.0) {
                  break;
               }
            }
            val += spin_multiple*poly_value_list_[i];
         }
      }
      else {
         for (std::size_t i = 0; i < num_interactions; ++i) {
            int32_t spin_multiple = 1;
            for (const auto &index: poly_key_list_[i]) {
               spin_multiple *= sample[index];
               if (std::abs(spin_multiple) <= 0.0) {
                  break;
               }
            }
            val += spin_multiple*poly_value_list_[i];
         }
      }
      return val;
   }
   
   PolynomialValueList<FloatType> energies(const std::vector<std::vector<int32_t>> &samples) {
      PolynomialValueList<FloatType> val_list(samples.size());
#pragma omp parallel for
      for (std::size_t i = 0; i < samples.size(); ++i) {
         val_list[i] = energy(samples[i], false);
      }
      return val_list;
   }
   
   void scale(const FloatType scalar,
              const PolynomialKeyList<IndexType> &ignored_interactions = {},
              const bool ignore_offset = false) {
      
      std::unordered_set<std::size_t> ignored_key_index_set;
      
      for (const auto &key: ignored_interactions) {
         if (poly_key_inv_.count(key) != 0) {
            ignored_key_index_set.emplace(poly_key_inv_[key]);
         }
      }
      
      std::size_t num_interactions = get_num_interactions();
      
      for (std::size_t i = 0; i < num_interactions; ++i) {
         if (ignored_key_index_set.count(i) != 0) {
            poly_value_list_[i] *= scalar;
         }
      }
      
      if (!ignore_offset && poly_key_inv_.count(std::vector<IndexType>{}) != 0) {
         poly_value_list_[poly_key_inv_[std::vector<IndexType>{}]] *= scalar;
      }
   }
   
   void normalize(const std::pair<FloatType, FloatType> &range,
                  const PolynomialKeyList<IndexType> &ignored_interactions = {},
                  const bool ignore_offset = false) {
      
      if (get_num_interactions() == 0) {
         return;
      }
      
      FloatType max_poly_value = poly_value_list_[0];
      FloatType min_poly_value = poly_value_list_[0];
      
      for (const auto &poly_value: poly_value_list_) {
         if (max_poly_value < poly_value) {
            max_poly_value = poly_value;
         }
         if (min_poly_value > poly_value) {
            min_poly_value = poly_value;
         }
      }
      FloatType scalar_for_min = range.first /min_poly_value;
      FloatType scalar_for_max = range.secons/max_poly_value;
      
      scale(std::min(scalar_for_min, scalar_for_max), ignored_interactions, ignore_offset);
      
   }
   
   BinaryPolynomialModel change_vartype(const Vartype &vartype, const bool inplace = true) {
      
      if (vartype == Vartype::SPIN) {
         if (inplace) {
            *this = ToSpin();
            return *this;
         }
         else {
            return ToSpin();
         }
      }
      else if (vartype == Vartype::BINARY) {
         if (inplace) {
            *this = ToBinary();
            return *this;
         }
         else {
            return ToBinary();
         }
      }
      else {
         throw std::runtime_error("Unknown vartype error");
      }  
   }

   Polynomial<IndexType, FloatType> to_hubo() {
      Polynomial<IndexType, FloatType> poly_map;
      std::size_t num_interactions = get_num_interactions();
      for (std::size_t i = 0; i < num_interactions; ++i) {
         const std::vector<IndexType>   original_key          = poly_key_list_[i];
         const FloatType original_value        = poly_value_list_[i];
         const std::size_t   original_key_size     = original_key.size();
         const std::size_t   changed_key_list_size = IntegerPower(2, original_key_size);
         
         for (std::size_t j = 0; j < changed_key_list_size; ++j) {
            const auto changed_key = GenerateChangedKey(original_key, j);
            int sign = ((original_key_size - changed_key.size())%2 == 0) ? 1.0 : -1.0;
            poly_map[changed_key] = original_value*IntegerPower(2, changed_key.size())*sign;
         }
      }
      return poly_map;
   }

   Polynomial<IndexType, FloatType> to_hising() {
      Polynomial<IndexType, FloatType> poly_map;
      std::size_t num_interactions = get_num_interactions();
      for (std::size_t i = 0; i < num_interactions; ++i) {
         const std::vector<IndexType>   original_key          = poly_key_list_[i];
         const FloatType original_value        = poly_value_list_[i];
         const std::size_t   original_key_size     = original_key.size();
         const std::size_t   changed_key_list_size = IntegerPower(2, original_key_size);
         const FloatType     changed_value         = original_value*(1.0/changed_key_list_size);
         
         for (std::size_t j = 0; j < changed_key_list_size; ++j) {
            poly_map[GenerateChangedKey(original_key , j)] = changed_value;
         }
      }
      return poly_map;
   }

   nlohmann::json to_serializable() const {
      nlohmann::json output;
      if (vartype_ == Vartype::BINARY) {
         output["vartype"] = "BINARY";
      }
      else if (vartype_ == Vartype::SPIN) {
         output["vartype"] = "SPIN";
      }
      else {
         throw std::runtime_error("Variable type must be SPIN or BINARY.");
      }

      output["poly_key_list"]   = poly_key_list_;
      output["poly_value_list"] = poly_value_list_;
      output["type"]            = "BinaryPolynomialModel";

      return output;
   }

   template <typename IndexType_serial = IndexType, typename FloatType_serial = FloatType>
   static BinaryPolynomialModel<IndexType_serial, FloatType_serial> from_serializable(nlohmann::json &input) {
      if(input["type"] != "BinaryPolynomialModel") {
         throw std::runtime_error("Type must be \"BinaryQuadraticModel\".\n");
      }

      Vartype vartype;
      if (input["vartype"] == "SPIN") {
         vartype = Vartype::SPIN;
      }
      else if (input["vartype"] == "BINARY") {
         vartype = Vartype::BINARY;
      }
      else {
         throw std::runtime_error("Variable type must be SPIN or BINARY.");
      }
      return BinaryPolynomialModel<IndexType_serial, FloatType_serial>(input["poly_key_list"], input["poly_value_list"], vartype);
   }

   static BinaryPolynomialModel from_hubo(const Polynomial<IndexType, FloatType> &poly_map) {
      return BinaryPolynomialModel<IndexType, FloatType>(poly_map, Vartype::BINARY);
   }
   
   static BinaryPolynomialModel from_hubo(const PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list) {
      return BinaryPolynomialModel<IndexType, FloatType>(key_list, value_list, Vartype::BINARY);
   }

   static BinaryPolynomialModel from_hubo(PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list) {
      return BinaryPolynomialModel<IndexType, FloatType>(key_list, value_list, Vartype::BINARY);
   }

   static BinaryPolynomialModel from_hising(const Polynomial<IndexType, FloatType> &poly_map) {
      return BinaryPolynomialModel<IndexType, FloatType>(poly_map, Vartype::SPIN);
   }

   static BinaryPolynomialModel from_hising(const PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list) {
      return BinaryPolynomialModel<IndexType, FloatType>(key_list, value_list, Vartype::SPIN);
   }

   static BinaryPolynomialModel from_hising(PolynomialKeyList<IndexType> &key_list, const PolynomialValueList<FloatType> &value_list) {
      return BinaryPolynomialModel<IndexType, FloatType>(key_list, value_list, Vartype::SPIN);
   }


protected:
   
   std::unordered_map<IndexType, std::size_t> each_variable_num_;
   
   PolynomialKeyList<IndexType> poly_key_list_;
   
   PolynomialValueList<FloatType> poly_value_list_;
   
   std::unordered_map<std::vector<IndexType>, std::size_t, vector_hash> poly_key_inv_;
   
   Vartype vartype_ = Vartype::NONE;
   
   std::string info_ = "";
   
   std::size_t degree_ = 0;
   
   void CheckKeySelfLoop(std::vector<IndexType> &key) const {
      //key is assumed to be sorted
      for (std::size_t i = 0; i < key.size() - 1; ++i) {
         if (key[i] == key[i + 1]) {
            throw std::runtime_error("No self-loops allowed");
         }
      }
   }
   
   void SetKeyAndValue(const std::vector<IndexType> &key, const FloatType &value) {
      if (poly_key_inv_.count(key) == 0) {
         poly_key_inv_[key] = poly_value_list_.size();
         poly_key_list_.push_back(key);
         poly_value_list_.push_back(value);
      }
      else {
         poly_value_list_[poly_key_inv_[key]] += value;
      }
      for (const auto &index: key) {
         each_variable_num_[index]++;
      }
   }
   
   std::size_t IntegerPower(std::size_t base, std::size_t exponent) const {
      std::size_t val = 1;
      for (std::size_t i = 0; i < exponent; ++i) {
         val *= base;
      }
      return val;
   }
   
   std::vector<IndexType> GenerateChangedKey(const std::vector<IndexType> &original_key, const std::size_t num_of_key) const {
      if (original_key.size() >= UINT16_MAX) {
         throw std::runtime_error("Too large degree of the interaction");
      }
      const std::size_t original_key_size = original_key.size();
      std::bitset<UINT16_MAX> bs(num_of_key);
      std::vector<IndexType> changed_key;
      for (std::size_t i = 0; i < original_key_size; ++i) {
         if (bs[i]) {
            changed_key.push_back(original_key[i]);
         }
      }
      return changed_key;
   }
   
   void UpdateDegree(std::size_t degree) {
      if (degree_ < degree) {
         degree_ = degree;
      }
   }

   BinaryPolynomialModel ToSpin() const {
      if (vartype_ == Vartype::SPIN) {
         return *this;
      }
      
      PolynomialKeyList<IndexType>   new_key_list;
      PolynomialValueList<FloatType> new_value_list;
      std::size_t num_interactions = get_num_interactions();
      
      for (std::size_t i = 0; i < num_interactions; ++i) {
         const std::vector<IndexType>   original_key          = poly_key_list_[i];
         const FloatType original_value        = poly_value_list_[i];
         const std::size_t   original_key_size     = original_key.size();
         const std::size_t   changed_key_list_size = IntegerPower(2, original_key_size);
         
         FloatType changed_value = original_value*(1.0/changed_key_list_size);
         
         for (std::size_t j = 0; j < changed_key_list_size; ++j) {
            new_key_list.push_back(GenerateChangedKey(original_key , j));
            new_value_list.push_back(changed_value);
         }
      }
      return BinaryPolynomialModel(new_key_list, new_value_list, Vartype::SPIN);
   }
   
   BinaryPolynomialModel ToBinary() const {
      
      if (vartype_ == Vartype::BINARY) {
         return *this;
      }
      
      PolynomialKeyList<IndexType>   new_key_list;
      PolynomialValueList<FloatType> new_value_list;
      std::size_t num_interactions = get_num_interactions();
      
      for (std::size_t i = 0; i < num_interactions; ++i) {
         const std::vector<IndexType>   original_key          = poly_key_list_[i];
         const FloatType original_value        = poly_value_list_[i];
         const std::size_t   original_key_size     = original_key.size();
         const std::size_t   changed_key_list_size = IntegerPower(2, original_key_size);
         
         for (std::size_t j = 0; j < changed_key_list_size; ++j) {
            const auto changed_key = GenerateChangedKey(original_key, j);
            int sign = ((original_key_size - changed_key.size())%2 == 0) ? 1.0 : -1.0;
            new_key_list.push_back(changed_key);
            new_value_list.push_back(original_value*IntegerPower(2, changed_key.size())*sign);
         }
      }
      return BinaryPolynomialModel(new_key_list, new_value_list, Vartype::BINARY);
   }
   
   
   
};

}


#endif /* binary_polynomial_model_hpp */
