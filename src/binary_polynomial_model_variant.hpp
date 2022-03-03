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

#ifndef binary_polynomial_model_variant_hpp
#define binary_polynomial_model_variant_hpp

#include "vartypes.hpp"
#include "hash.hpp"
#include "utilities.hpp"

#include <unordered_map>
#include <variant>

namespace cimod {

template <typename FloatType>
class BinaryPolynomialModelVariant {
   
public:
   using IndexType = std::variant<std::int64_t, std::string, std::vector<std::variant<std::int64_t, std::string>>>;
   
   BinaryPolynomialModelVariant(const Vartype vartype): vartype_(vartype) {}
   
   
   void AddInteraction(std::vector<IndexType> &key, const FloatType value) {
      FormatPolynomialKey(&key, vartype_);
      poly_map_[key] += value;
      variables_.insert(key.begin(), key.end());
   }
      
   void AddInteraction(std::vector<IndexType> &&key, const FloatType value) {
      FormatPolynomialKey(&key, vartype_);
      poly_map_[key] += value;
      variables_.insert(key.begin(), key.end());
   }
   
   void AddInteractionsFrom(const std::unordered_map<std::vector<IndexType>, FloatType, VariantVectorHash> &map) {
      for (const auto &it: map) {
         //auto &v = const_cast<std::vector<IndexType>&>(it.first);
         auto v = it.first;
         FormatPolynomialKey(&v, vartype_);
         poly_map_[v] += it.second;
         variables_.insert(v.begin(), v.end());
      }
   }
   
   void AddInteractionsFrom(std::vector<std::vector<IndexType>> &keys, const std::vector<FloatType> &values) {
      if (keys.size() != values.size()) {
         throw std::runtime_error("The size of keys and values must be equal");
      }
      for (std::size_t i = 0; i < keys.size(); ++i) {
         FormatPolynomialKey(&keys[i], vartype_);
         poly_map_[keys[i]] += values[i];
         variables_.insert(keys[i].begin(), keys[i].end());
      }
   }
   
   void AddInteractionsFrom(std::vector<std::pair<std::vector<IndexType>, FloatType>> &pairs) {
      for (std::size_t i = 0; i < pairs.size(); ++i) {
         FormatPolynomialKey(&pairs[i].first, vartype_);
         poly_map_[pairs[i].first] += pairs[i].second;
         variables_.insert(pairs[i].first.begin(), pairs[i].first.end());
      }
   }
   
  
private:
   std::unordered_set<IndexType, VariantHash> variables_;
   std::unordered_map<std::vector<IndexType>, FloatType, VariantVectorHash> poly_map_;
   
   std::unordered_set<std::int64_t> variables_i_;
   std::unordered_map<std::vector<std::int64_t>, FloatType, vector_hash> poly_map_i_;
   
   Vartype vartype_ = Vartype::NONE;
   
   
   
};

}


#endif /* binary_polynomial_model_hpp */
