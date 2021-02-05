//
//  binary_polynomial_model.h
//  cimod_Xcode
//
//  Created by Kohei Suzuki on 2021/02/05.
//

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

/**
 * @brief Type alias for linear bias
 *
 * @tparam IndexType
 */
template <typename IndexType, typename FloatType>
using Linear = std::unordered_map<IndexType, FloatType>;

/**
 * @brief Type alias for quadratic bias
 *
 * @tparam IndexType
 */
template <typename IndexType, typename FloatType>
using Polynominal = std::unordered_map<std::vector<IndexType>, FloatType, vector_hash>;

/**
 * @brief Type alias for adjacency list
 *
 * @tparam IndexType
 */
template <typename IndexType, typename FloatType>
using Adjacency_Poly = std::unordered_map<IndexType, std::unordered_map<std::vector<IndexType>, FloatType, vector_hash>>;

/**
 * @brief Type alias for sample
 *
 * @tparam IndexType
 */
template <typename IndexType>
using Sample = std::unordered_map<IndexType, int32_t>;

/**
 * @brief Class for binary quadratic model.
 */
template <typename IndexType, typename FloatType>
class BinaryPolynomialModel {
   
protected:
   /**
    * @brief Linear biases as a dictionary.
    *
    */
   Linear<IndexType, FloatType> m_linear;
   
   /**
    * @brief Polynominal biases as a dictionary.
    *
    */
   Polynominal<IndexType, FloatType> m_polynominal;
   
   /**
    * @brief The energy offset associated with the model.
    *
    */
   FloatType m_offset;
   
   /**
    * @brief The model's type.
    *
    */
   Vartype m_vartype = Vartype::NONE;
   
   /**
    * @brief A place to store miscellaneous data about the binary Polynominal model as a whole.
    *
    */
   std::string m_info;
   
   /**
    * @brief The model's interactions as nested dictionaries.
    *
    */
   Adjacency_Poly<IndexType, FloatType> m_adj;
   
   /**
    * @brief Add the adjacency to the adjacency list
    *
    */
   void update_adjacency(const std::vector<IndexType> &u) {
      if(m_polynominal.count(u)!=0) {
         int min = *std::min_element(u.begin(), u.end());
         insert_or_assign(m_adj[min], u, m_polynominal[u]);
      }
   }
   
   /* FIX ME
    void remove_adjacency(const IndexType &u, const IndexType &v) {
    auto k = m_adj[u].erase(v);
    };
    */
   
public:
   /**
    * @brief BinaryQuadraticModel constructor.
    *
    * @param linear
    * @param polynominal
    * @param offset
    * @param vartype
    * @param info
    */
   BinaryPolynomialModel
   (
    const Linear<IndexType, FloatType> &linear,
    const Polynominal<IndexType, FloatType> &polynominal,
    const FloatType &offset,
    const Vartype vartype,
    const std::string info = ""
    ):
   m_offset(offset),
   m_vartype(vartype),
   m_info(info)
   {
      add_variables_from(linear);
      add_interactions_from(polynominal);
   };
   
   BinaryPolynomialModel(const BinaryPolynomialModel&) = default;
   BinaryPolynomialModel(BinaryPolynomialModel&&) = default;
   
   /**
    * @brief Return the number of variables.
    *
    * @return The number of variables.
    */
   size_t length() const {
      return m_linear.size();
   };
   
   /**
    * @brief Return true if the variable contains v.
    *
    * @return Return true if the variable contains v.
    * @param v
    */
   bool contains(const IndexType &v) {
      if(m_linear.count(v)!=0) {
         return true;
      }
      else {
         return false;
      }
   };
   
   const Linear<IndexType, FloatType>& get_linear() const {
      return this->m_linear;
   };
   
   const Polynominal<IndexType, FloatType>& get_polynominal() const {
      return this->m_polynominal;
   };
   
   const Adjacency_Poly<IndexType, FloatType>& get_adjacency() const {
      return this->m_adj;
   };
   
   const FloatType& get_offset() const {
      return this->m_offset;
   };
   
   const Vartype& get_vartype() const {
      return this->m_vartype;
   };
   
   const std::string& get_info() const {
      return this->m_info;
   };
   
   void print() {
      std::cout << "[BinaryPolynomialModel]" << std::endl;
      
      // Print linear
      std::cout << "linear = " << std::endl;
      for(auto &it : m_linear) {
         std::cout << "" << it.first << ": " << it.second << std::endl;
      }
      
      
      // Print polynominal
      std::cout << "polynominal = " << std::endl;
      for(auto &it : m_polynominal)
      {
         std::cout << "(";
         for (auto &it_index : it.first) {
            std::cout << it_index << ", ";
         }
         std::cout << "): " << it.second << ", ";;
      }
      std::cout << std::endl;
      
      
      // Print adjacency
      std::cout << "adjacency = " << std::endl;
      for(auto &it_src : m_linear) {
         std::cout << it_src.first << ": {";
         for (auto &it_dst : m_adj[it_src.first]) {
            std::cout << "(";
            for (auto &it_vec : it_dst.first) {
               std::cout << it_vec << ", ";
            }
            std::cout << "): " << it_dst.second << ", ";
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
      
      // Print info
      std::cout << "info = ";
      std::cout << "\"" << m_info << "\"" << std::endl;
      
   };
   
   void empty() {
      m_linear = {};
      m_polynominal = {};
      m_offset = 0.0;
      m_vartype = Vartype::NONE;
      m_info = "";
   }
   
   void add_variable
   (
    const IndexType &v,
    const FloatType &bias,
    const Vartype vartype = Vartype::NONE
    ) {
      FloatType b = bias;
      
      // handle the case that a different vartype is provided
      if((vartype!=Vartype::NONE)&&(vartype!=m_vartype)) {
         if((m_vartype == Vartype::SPIN)&&(vartype == Vartype::BINARY))
         {
            b /= 2;
            m_offset += b;
         }
         else if((m_vartype == Vartype::BINARY)&&(vartype == Vartype::SPIN))
         {
            m_offset -= b;
            b *= 2;
         }
         else
         {
            std::cerr << "Unknown vartype" << std::endl;
         }
      }
      
      // Insert or assign the bias
      FloatType value = 0;
      if(m_linear.count(v) != 0) {
         value = m_linear[v];
      }
      insert_or_assign(m_linear, v, value + b);
   };
   
   void add_variables_from
   (
    const Linear<IndexType, FloatType> &linear,
    const Vartype vartype = Vartype::NONE
    ) {
      for(auto &it : linear)
      {
         add_variable(it.first, it.second, vartype);
      }
   };
   
   void add_interaction(const std::vector<IndexType> &u, const FloatType &bias, const Vartype vartype = Vartype::NONE) {
      
      //Check the input interaction is valid
      if (u.size() == 1) {
         std::cerr << "On-site potential" << " (" << u[0] << "): " << bias << " is not allowed" << std::endl;
         return;
      }
      
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
      
      FloatType b = bias;
      
      if((vartype!=Vartype::NONE) && (vartype!=m_vartype)) {
         if((m_vartype == Vartype::SPIN) && (vartype == Vartype::BINARY)) {
            //convert from binary to spin
            b /= 4;
            add_offset(b);
            for (auto &it : u) { add_variable(it, b); }
         }
         else if((m_vartype == Vartype::BINARY)&&(vartype == Vartype::SPIN)) {
            //convert from spin to binary
            add_offset(b);
            for (auto &it : u) { add_variable(it, -2 * b); }
            b *= 4;
         }
         else {
            std::cerr << "Unknown vartype" << std::endl;
         }
      }
      else {
         for (auto &it : u) {
            if (m_linear.count(it) == 0) {
               add_variable(it, 0);
            }
         }
      }
      
      FloatType value = 0;
      if (m_polynominal.count(u) != 0) {
         value = m_polynominal[u];
      }
      insert_or_assign(m_polynominal, u, value + b);
      update_adjacency(u);
      
   };
   
   void add_interactions_from
   (
    const Polynominal<IndexType, FloatType> &polynominal,
    const Vartype vartype = Vartype::NONE
    ) {
      for(auto &it : polynominal) {
         add_interaction(it.first, it.second, vartype);
      }
   }
   
   void remove_variable(const IndexType &v) {
      std::vector<std::pair<IndexType, IndexType>> interactions;
      for(auto &it : m_polynominal) {
         if(it.first.first == v || it.first.second == v) {
            interactions.push_back(it.first);
         }
      }
      remove_interactions_from(interactions);
      m_linear.erase(v);
      m_adj.erase(v);
   }
   
   void remove_variables_from(const std::vector<IndexType> &variables) {
      for(auto &it : variables) {
         remove_variable(it);
      }
   };
   
   /*
    void remove_interaction
    (
    const IndexType &u,
    const IndexType &v
    )
    {
    auto p = std::make_pair(u, v);
    if(m_polynominal.count(p)!=0)
    {
    auto k = m_polynominal.erase(p);
    remove_adjacency(u, v);
    }
    };
    
    void remove_interactions_from
    (
    const std::vector<std::pair<IndexType, IndexType>> &interactions
    )
    {
    for(auto &it : interactions)
    {
    remove_interaction(it.first, it.second);
    }
    };
    */
   
   void add_offset(const FloatType &offset) {
      m_offset += offset;
   };
   
   void remove_offset() {
      add_offset(-m_offset);
   };
   
   /*
    void scale
    (
    const FloatType &scalar,
    const std::vector<IndexType> &ignored_variables = {},
    const std::vector<std::pair<IndexType, IndexType>> &ignored_interactions = {},
    const bool ignored_offset = false
    )
    {
    // scaling linear
    for(auto &it : m_linear)
    {
    if(std::find(ignored_variables.begin(), ignored_variables.end(), it.first) != ignored_variables.end() || ignored_variables.empty())
    {
    it.second *= scalar;
    }
    }
    
    // scaling polynominal
    for(auto &it : m_polynominal)
    {
    if(std::find(ignored_interactions.begin(), ignored_interactions.end(), it.first) != ignored_interactions.end() || ignored_variables.empty())
    {
    it.second *= scalar;
    }
    }
    
    // scaling offset
    if(!ignored_offset)
    {
    m_offset *= scalar;
    }
    };
    */
   
   /*
    void normalize
    (
    const std::pair<FloatType, FloatType> &bias_range = {1.0, 1.0},
    const bool use_quadratic_range = false,
    const std::pair<FloatType, FloatType> &quadratic_range = {1.0, 1.0},
    const std::vector<IndexType> &ignored_variables = {},
    const std::vector<std::pair<IndexType, IndexType>> &ignored_interactions = {},
    const bool ignored_offset = false
    )
    {
    // parse range
    std::pair<FloatType, FloatType> l_range = bias_range;
    std::pair<FloatType, FloatType> q_range;
    if(!use_quadratic_range)
    {
    q_range = bias_range;
    }
    else
    {
    q_range = quadratic_range;
    }
    
    // calculate scaling value
    auto comp = [](const auto &a, const auto &b) { return a.second < b.second; };
    auto it_lin_min = std::min_element(m_linear.begin(), m_linear.end(), comp);
    auto it_lin_max = std::max_element(m_linear.begin(), m_linear.end(), comp);
    auto it_quad_min = std::min_element(m_polynominal.begin(), m_polynominal.end(), comp);
    auto it_quad_max = std::max_element(m_polynominal.begin(), m_polynominal.end(), comp);
    
    std::vector<FloatType> v_scale =
    {
    it_lin_min->second / l_range.first,
    it_lin_max->second / l_range.second,
    it_quad_min->second / q_range.first,
    it_quad_max->second / q_range.second
    };
    FloatType inv_scale = *std::max_element(v_scale.begin(), v_scale.end());
    
    // scaling
    if(inv_scale != 0.0)
    {
    scale(1.0 / inv_scale, ignored_variables, ignored_interactions, ignored_offset);
    }
    };
    */
   
   void fix_variable(const IndexType &v, const int32_t &value) {
      std::vector<std::pair<IndexType, IndexType>> interactions;
      for(auto &it : m_polynominal) {
         if(it.first.first == v) {
            add_variable(it.first.second, value*it.second);
            interactions.push_back(it.first);
         }
         else if(it.first.second == v) {
            add_variable(it.first.first, value*it.second);
            interactions.push_back(it.first);
         }
      }
      remove_interactions_from(interactions);
      add_offset(m_linear[v]*value);
      remove_variable(v);
   };
   
   void fix_variables(const std::vector<std::pair<IndexType, int32_t>> &fixed) {
      for(auto &it : fixed) {
         fix_variable(it.first, it.second);
      }
   };
   
   
   void flip_variable(const IndexType &v) {
      // check variable
      if(m_linear.count(v)==0) {
         std::cout << v << " is not a variable in the binary polynominal model." << std::endl;
         return;
      }
      
      if(m_vartype == Vartype::SPIN) {
         m_linear[v] *= -1.0;
         for(auto &it : m_polynominal) {
            if(it.first.first == v || it.first.second == v) {
               it.second *= -1.0;
               update_adjacency(it.first.first, it.first.second);
            }
         }
      }
      else if(m_vartype == Vartype::BINARY) {
         add_offset(m_linear[v]);
         m_linear[v] *= -1.0;
         
         for(auto &it : m_polynominal) {
            if(it.first.first == v) {
               m_linear[it.first.second] += it.second;
               it.second *= -1.0;
               update_adjacency(it.first.first, it.first.second);
            }
            else if(it.first.second == v) {
               m_linear[it.first.first] += it.second;
               it.second *= -1.0;
               update_adjacency(it.first.first, it.first.second);
            }
         }
      }
   };
   
   
   void update(const BinaryPolynomialModel &bpm, const bool ignore_info = true) {
      add_variables_from(bpm.get_linear(), bpm.get_vartype());
      add_interactions_from(bpm.get_quadratic(), bpm.get_vartype());
      add_offset(bpm.get_offset());
      if(!ignore_info) {
         m_info = bpm.get_info();
      };
   };
};

}


#endif /* binary_polynomial_model_hpp */
