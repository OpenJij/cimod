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
 * @file binary_quadratic_model_dense.hpp
 * @author Kohji Nishimura
 * @brief Dense BinaryQuadraticModel class
 * @version 1.0.0
 * @date 2020-03-24
 * 
 * @copyright Copyright (c) Jij Inc. 2020
 * 
 */

#ifndef BINARY_QUADRATIC_MODEL_DENSE_HPP__
#define BINARY_QUADRATIC_MODEL_DENSE_HPP__

#include "disable_eigen_warning.hpp"

#include "vartypes.hpp"
#include "hash.hpp"
#include "utilities.hpp"
#include "json.hpp"

#include <algorithm>
#include <limits>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <functional>
#include <Eigen/Dense>

#include "binary_quadratic_model.hpp"

namespace cimod
{

/**
 * @brief Class for dense binary quadratic model.
 */

template <typename IndexType, typename FloatType>
class BinaryQuadraticModel_Dense
{
public:
/**
 * @brief Eigen Matrix
 */
using Matrix = Eigen::Matrix<FloatType, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Vector = Eigen::Matrix<FloatType, Eigen::Dynamic, 1>;

protected:
    
    /**
     * @brief quadratic dense-type matrix
     * The stored matrix has the following triangular form:
     *
     * \f[
     * \begin{pmatrix}
     * J_{0,0} & J_{0,1} & \cdots & J_{0,N-1} & h_{0}\\
     * 0 & J_{1,1} & \cdots & J_{1,N-1} & h_{1}\\
     * \vdots & \vdots & \vdots & \vdots & \vdots \\
     * 0 & 0 & \cdots & J_{N-1,N-1} & h_{N-1}\\
     * 0 & 0 & \cdots & 0 & 1 \\
     * \end{pmatrix}
     * \f]
     */
    Matrix _quadmat;

    /**
     * @brief vector for converting index to label
     * the list is asssumed to be sorted
     */
    std::vector<IndexType> _idx_to_label;

    /**
     * @brief dict for converting label to index
     */
    std::unordered_map<IndexType, size_t> _label_to_idx;

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
     * @brief set _label_to_idx from _idx_to_label
     */
    inline void _set_label_to_idx(){
        //reset
        _label_to_idx.clear();
        //initialize
        for(size_t i=0; i<_idx_to_label.size(); i++){
            _label_to_idx[_idx_to_label[i]] = i;
        }
    }

    /**
     * @brief get reference of _quadmat(i,j)
     *
     * @param label_i
     * @param label_j
     *
     * @return reference of _quadmat(i,j)
     */
    inline FloatType& _mat(IndexType label_i, IndexType label_j){
        size_t i = _label_to_idx.at(label_i);
        size_t j = _label_to_idx.at(label_j);
        if(i != j)
            return _quadmat(std::min(i, j), std::max(i, j));
        else
            throw std::runtime_error("No self-loop (mat(i,i)) allowed");
    }

    /**
     * @brief get reference of _quadmat(i,i)
     *
     * @param label_i
     *
     * @return reference of _quadmat(i,i)
     */
    inline FloatType& _mat(IndexType label_i){
        size_t i = _label_to_idx.at(label_i);
        return _quadmat(i, _quadmat.rows()-1);
    }

    /**
     * @brief get reference of _quadmat(i,j)
     *
     * @param label_i
     * @param label_j
     *
     * @return reference of _quadmat(i,j)
     */
    const inline FloatType& _mat(IndexType label_i, IndexType label_j) const{
        size_t i = _label_to_idx.at(label_i);
        size_t j = _label_to_idx.at(label_j);

        if(i != j)
            return _quadmat(std::min(i, j), std::max(i, j));
        else
            throw std::runtime_error("No self-loop (mat(i,i)) allowed");
    }

    /**
     * @brief get reference of _quadmat(i,i)
     *
     * @param label_i
     *
     * @return reference of _quadmat(i,i)
     */
    const inline FloatType& _mat(IndexType label_i) const{
        size_t i = _label_to_idx.at(label_i);
        return _quadmat(i, _quadmat.rows()-1);
    }

    /**
     * @brief calculate maximum element in linear term
     *
     * @return 
     */
    const inline FloatType _max_linear() const{
        size_t N = _quadmat.rows();
        return _quadmat.block(0,N-1,N-1,1).maxCoeff();
    }

    /**
     * @brief calculate maximum element in quadratic term
     *
     * @return 
     */
    const inline FloatType _max_quadratic() const{
        size_t N = _quadmat.rows();
        return _quadmat.block(0,0,N-1,N-1).maxCoeff();
    }

    /**
     * @brief calculate minimum element in linear term
     *
     * @return 
     */
    const inline FloatType _min_linear() const{
        size_t N = _quadmat.rows();
        return _quadmat.block(0,N-1,N-1,1).minCoeff();
    }

    /**
     * @brief calculate minimum element in quadratic term
     *
     * @return 
     */
    const inline FloatType _min_quadratic() const{
        size_t N = _quadmat.rows();
        return _quadmat.block(0,0,N-1,N-1).minCoeff();
    }

    /**
     * @brief insert row and column that corresponds to added label into _quadmat
     *
     * @param label_i
     */
    inline void _insert_label_into_mat(IndexType label_i){
        size_t i = _label_to_idx.at(label_i);
        //define temp mat
        size_t N = _quadmat.rows()+1;
        Matrix tempmat = Matrix::Zero(N, N);
        //copy elements to new matrix
        tempmat.block(0,0,i,i)              = _quadmat.block(0,0,i,i);
        tempmat.block(0,i+1,i,N-i-1)        = _quadmat.block(0,i,i,N-i-1);
        tempmat.block(i+1,i+1,N-i-1,N-i-1)  = _quadmat.block(i,i,N-i-1,N-i-1);

        _quadmat = tempmat;
    }

    /**
     * @brief delete row and column that corresponds to existing label from _quadmat
     *
     * @param label_i
     */
    inline void _delete_label_from_mat(IndexType label_i){
        size_t i = _label_to_idx.at(label_i);
        //define temp mat
        size_t N = _quadmat.rows();
        Matrix tempmat = Matrix::Zero(N-1, N-1);
        //copy elements to new matrix
        tempmat.block(0,0,i,i)              = _quadmat.block(0,0,i,i);
        tempmat.block(0,i,i,N-i-1)          = _quadmat.block(0,i+1,i,N-i-1);
        tempmat.block(i,i,N-i-1,N-i-1)      = _quadmat.block(i+1,i+1,N-i-1,N-i-1);

        _quadmat = tempmat;
    }

    /**
     * @brief add new label
     * if label_i already exists, this process is skipped.
     *
     * @param label_i
     */
    inline void _add_new_label(IndexType label_i){
        if(_label_to_idx.find(label_i) == _label_to_idx.end()){
            //add label_i
            _idx_to_label.push_back(label_i);
            std::sort(_idx_to_label.begin(), _idx_to_label.end());
            _set_label_to_idx();

            _insert_label_into_mat(label_i);
        }
    }

    /**
     * @brief delete label
     * if label_i does not exist, this process is skipped.
     *
     * @param label_i
     * @param force_delete if true, delete label whenever there exists corresponding nonzero elements in the matrix.
     * otherwise, delete label only if there are no corresponding nonzero elements in the matrix.
     */
    inline void _delete_label(IndexType label_i, bool force_delete = true){
        auto position = std::find(_idx_to_label.begin(), _idx_to_label.end(), label_i);
        if(position != _idx_to_label.end()){
            if(force_delete == false){
                //check if there are corresponding nonzero elements
                size_t i = std::distance(_idx_to_label.begin(), position);
                if(_quadmat.col(i).squaredNorm() > std::numeric_limits<FloatType>::epsilon() ||
                   _quadmat.row(i).squaredNorm() > std::numeric_limits<FloatType>::epsilon() ){
                    // exists nonzero elements
                    return;
                }
            }
            //delete from matrix first
            _delete_label_from_mat(label_i);
            //add label_i
            _idx_to_label.erase(position);
            // already sorted
            //std::sort(_idx_to_label.begin(), _idx_to_label.end());
            _set_label_to_idx();
        }
    }

    /**
     * @brief initialize matrix with linear and quadratic dicts
     *
     * @param linear
     * @param quadratic
     */
    inline void _initialize_quadmat(const Linear<IndexType, FloatType> &linear, const Quadratic<IndexType, FloatType> &quadratic){
        //gather labels
        std::unordered_set<IndexType> labels;
        
        for(const auto& kv : linear){
            labels.insert(kv.first);
        }

        for(const auto& kv : quadratic){
            labels.insert(kv.first.first);
            labels.insert(kv.first.second);
        }

        // init label <-> index conversion variables
        _idx_to_label = std::vector<IndexType>(labels.begin(), labels.end());
        std::sort(_idx_to_label.begin(), _idx_to_label.end());
        _set_label_to_idx();

        //initialize _quadmat
        size_t mat_size = _idx_to_label.size() + 1;
        _quadmat = Matrix(mat_size, mat_size);
        _quadmat.fill(0);
        _quadmat(mat_size-1, mat_size-1) = 1;

        //copy linear and quadratic to _quadmat
        for(const auto& kv : linear){
            IndexType key = kv.first;
            FloatType val = kv.second;
            _mat(key) += val;
        }

        for(const auto& kv : quadratic){
            std::pair<IndexType, IndexType> key = kv.first;
            FloatType val = kv.second;
            _mat(key.first, key.second) += val;
        }
    }

    /**
     * @brief initialize matrix with matrix and labels
     * the form of matrix is assumed to be the following two forms:
     *
     * \f[
     * \begin{pmatrix}
     * J_{0,0} & J_{0,1} & \cdots & J_{0,N-1} & h_{0}\\
     * J_{1,0} & J_{1,1} & \cdots & J_{1,N-1} & h_{1}\\
     * \vdots & \vdots & \vdots & \vdots & \vdots \\
     * J_{N-1,0} & J_{N-1,1} & \cdots & J_{N-1,N-1} & h_{N-1}\\
     * h_{0} & h_{1} & \cdots & h_{N-1} & 1 \\
     * \end{pmatrix}
     * \f]
     *
     * or
     *
     * \f[
     * \begin{pmatrix}
     * h_{0} & J_{0,1} & \cdots & J_{0,N-1}\\
     * J_{1,0} & h_{1} & \cdots & J_{1,N-1}\\
     * \vdots & \vdots & \vdots & \vdots\\
     * J_{N-1,0} & J_{N-1,1} & \cdots & h_{N-1}\\
     * \end{pmatrix}
     * \f]
     *
     *
     * @param mat
     * @param labels
     */
    inline void _initialize_quadmat(const Matrix& mat, const std::vector<IndexType>& labels_vec){

        //initlaize label <-> index dict
        std::unordered_set<IndexType> labels(labels_vec.begin(), labels_vec.end());
        _idx_to_label = std::vector<IndexType>(labels.begin(), labels.end());
        std::sort(_idx_to_label.begin(), _idx_to_label.end());
        _set_label_to_idx();

        if(mat.rows() != mat.cols()){
            throw std::runtime_error("matrix must be a square matrix");
        }

        if((size_t)mat.rows() == _idx_to_label.size() + 1){
            size_t mat_size = _idx_to_label.size() + 1;
            _quadmat = Matrix::Zero(mat_size, mat_size);
            //insert elements
            _quadmat += mat.template triangularView<Eigen::StrictlyUpper>();
            _quadmat += mat.template triangularView<Eigen::StrictlyLower>().transpose();
            _quadmat(mat_size-1, mat_size-1) = 1;
        }
        else if((size_t)mat.rows() == _idx_to_label.size()){
            //convert matrix
            size_t mat_size = _idx_to_label.size() + 1;
            _quadmat = Matrix::Zero(mat_size, mat_size);
            _quadmat.block(0,0,mat_size-1,mat_size-1) += mat.template triangularView<Eigen::StrictlyUpper>();
            _quadmat.block(0,0,mat_size-1,mat_size-1) += mat.template triangularView<Eigen::StrictlyLower>().transpose();
            //local fields
            Vector loc = mat.diagonal();
            _quadmat.block(0,mat_size-1,mat_size-1,1) += loc;
            _quadmat(mat_size-1, mat_size-1) = 1;
        }
        else{
            throw std::runtime_error("the number of variables and dimension do not match.");
        }
    }


    inline Linear<IndexType, FloatType> _generate_linear() const{
        Linear<IndexType, FloatType> ret_linear;
        for(size_t i=0; i<_idx_to_label.size(); i++){
            FloatType val = _quadmat(i, _idx_to_label.size());
            if(val != 0)
                ret_linear[_idx_to_label[i]] = val;
        }

        return ret_linear;
    }

    inline Quadratic<IndexType, FloatType> _generate_quadratic() const{
        Quadratic<IndexType, FloatType> ret_quadratic;
        for(size_t i=0; i<_idx_to_label.size(); i++){
            for(size_t j=i+1; j<_idx_to_label.size(); j++){
                FloatType val = _quadmat(i, j);
                if(val != 0)
                    ret_quadratic[std::make_pair(_idx_to_label[i], _idx_to_label[j])] = val;
            }
        }

        return ret_quadratic;
    }



    /**
     * @brief change internal variable from Ising to QUBO ones
     * The following conversion is applied:
     *
     * \f[
     * \mathrm{offset} += \sum_{i<j} J_{ij} - \sum_{i}h_{i}
     * \f]
     * \f[
     * Q_ii += -2\left(\sum_{j}J_{ji}+\sum_{j}J_{ij}\right) + 2h_{i}
     * \f]
     * \f[
     * Q_{ij} = 4J_{ij}
     * \f]
     */
    inline void _spin_to_binary(){
        size_t num_variables = _idx_to_label.size();
        m_vartype = Vartype::BINARY;
        //calc col(row)wise-sum ((num_variables, 1))
        Vector colwise_sum = _quadmat.block(0,0,num_variables,num_variables).colwise().sum();
        Vector rowwise_sum = _quadmat.block(0,0,num_variables,num_variables).rowwise().sum();
        
        Vector local_field = _quadmat.block(0,num_variables,num_variables,1);

        //offset
        m_offset += colwise_sum.sum() - local_field.sum();

        //local field
        _quadmat.block(0,num_variables,num_variables,1)
            = 2 * local_field - 2 * (colwise_sum + rowwise_sum);

        //quadratic
        _quadmat.block(0,0,num_variables,num_variables) *= 4;
    }

    /**
     * @brief change internal variable from QUBO to Ising ones
     * The following conversion is applied:
     *
     * \f[
     * \mathrm{offset} += \frac{1}{4}\sum_{i<j} Q_{ij} + \frac{1}{2}\sum_{i}Q_{ii}
     * \f]
     * \f[
     * h_i += \frac{1}{4}\left(\sum_{j}Q_{ji}+\sum_{j}Q_{ij}\right) + \frac{1}{2}Q_{ii}
     * \f]
     * \f[
     * J_{ij} = \frac{1}{4}Q_{ij}
     * \f]
     */
    inline void _binary_to_spin(){
        size_t num_variables = _idx_to_label.size();
        m_vartype = Vartype::SPIN;
        //calc col(row)wise-sum ((num_variables, 1))
        Vector colwise_sum = _quadmat.block(0,0,num_variables,num_variables).colwise().sum();
        Vector rowwise_sum = _quadmat.block(0,0,num_variables,num_variables).rowwise().sum();
        
        Vector local_field = _quadmat.block(0,num_variables,num_variables,1);

        //offset
        m_offset += 0.25 * colwise_sum.sum() + 0.5 * local_field.sum();

        //local field
        _quadmat.block(0,num_variables,num_variables,1)
            = 0.5 * local_field + 0.25 * (colwise_sum + rowwise_sum);

        //quadratic
        _quadmat.block(0,0,num_variables,num_variables) *= 0.25;
    }


public:
    /**
     * @brief BinaryQuadraticModel_Dense constructor.
     * 
     * @param linear
     * @param quadratic
     * @param offset
     * @param vartype
     */
    BinaryQuadraticModel_Dense
    (
        const Linear<IndexType, FloatType> &linear,
        const Quadratic<IndexType, FloatType> &quadratic,
        const FloatType &offset,
        const Vartype vartype
    ):
        m_offset(offset),
        m_vartype(vartype)
    {
        _initialize_quadmat(linear, quadratic);
    }

    /**
     * @brief BinaryQuadraticModel_Dense constructor.
     *
     * @param linear
     * @param quadratic
     * @param vartype
     *
     */
    BinaryQuadraticModel_Dense
    (
        const Linear<IndexType, FloatType> &linear,
        const Quadratic<IndexType, FloatType> &quadratic,
        const Vartype vartype
    ): BinaryQuadraticModel_Dense(linear, quadratic, 0.0, vartype){}


    /**
     * @brief BinaryQuadraticModel_Dense constructor (with matrix);
     *
     * @param mat
     * @param labels_vec
     * @param offset
     * @param vartype
     *
     */
    BinaryQuadraticModel_Dense
    (
        const Matrix& mat,
        const std::vector<IndexType>& labels_vec,
        const FloatType &offset,
        const Vartype vartype
    ):
        m_offset(offset),
        m_vartype(vartype)
    {
        _initialize_quadmat(mat, labels_vec);
    }

    /**
     * @brief BinaryQuadraticModel_Dense constructor (with matrix);
     *
     * @param mat
     * @param labels_vec
     * @param vartype
     *
     */
    BinaryQuadraticModel_Dense
    (
        const Matrix& mat,
        const std::vector<IndexType>& labels_vec,
        const Vartype vartype
    ): BinaryQuadraticModel_Dense(mat, labels_vec, 0.0, vartype){}


    BinaryQuadraticModel_Dense(const BinaryQuadraticModel_Dense&) = default;

    /**
     * @brief get the number of variables
     *
     * @return The number of variables.
     */
    size_t get_num_variables() const{
        return _idx_to_label.size();
    }

    /**
     * @brief Return the number of variables.
     * @deprecated use get_num_variables instead.
     * 
     * @return The number of variables.
     */
    size_t length() const
    {
        return get_num_variables();
    }

    /**
     * @brief Return true if the variable contains v.
     * 
     * @return Return true if the variable contains v.
     * @param v
     */
    bool contains(const IndexType &v) const
    {
        if(_label_to_idx.find(v) != _label_to_idx.end())
        {
            return true;
        }
        else
        {
            return false;
        }
        
    }
    
    /**
     * @brief Get the element of linear object
     * 
     * @return A linear bias.
     */
    FloatType get_linear(IndexType label_i) const{
        return _mat(label_i);
    }

    /**
     * @brief Get linear object
     *
     * @return A linear object
     */
    Linear<IndexType, FloatType> get_linear() const{
        return _generate_linear();
    }

    /**
     * @brief Get the element of quadratic object
     * 
     * @return A quadratic bias.
     */
    FloatType get_quadratic(IndexType label_i, IndexType label_j) const
    {
        return _mat(label_i, label_j);
    }

    /**
     * @brief Get uadratic object
     * 
     * @return A quadratic object.
     */
    Quadratic<IndexType, FloatType> get_quadratic() const
    {
        return _generate_quadratic();
    }

    /**
     * @brief Get the offset
     * 
     * @return An offset.
     */
    FloatType get_offset() const
    {
        return this->m_offset;
    }

    /**
     * @brief Get the vartype object
     * 
     * @return Type of the model.
     */
    Vartype get_vartype() const
    {
        return this->m_vartype;
    }

    /**
     * @brief Get variables
     *
     * @return variables
     */
    const std::vector<IndexType>& get_variables() const{
        return this->_idx_to_label;
    }


    /**
     * @brief Create an empty BinaryQuadraticModel
     * 
     * @return empty object
     */
    BinaryQuadraticModel_Dense<IndexType, FloatType> empty(Vartype vartype)
    {
        return BinaryQuadraticModel_Dense<IndexType, FloatType>(
                Linear<FloatType, IndexType>(),
                Quadratic<FloatType, IndexType>(),
                0.0,
                vartype
                );
    }

    /* Update methods */

    /**
     * @brief Add variable v and/or its bias to a binary quadratic model.
     * 
     * @param v
     * @param bias
     */
    void add_variable
    (
        const IndexType &v,
        const FloatType &bias
    )
    {
        //add new label if not exist
        _add_new_label(v);
        _mat(v) += bias;
    }

    /**
     * @brief Add variables and/or linear biases to a binary quadratic model.
     * 
     * @param linear
     */
    void add_variables_from
    (
        const Linear<IndexType, FloatType> &linear
    )
    {
        for(auto &it : linear)
        {
            add_variable(it.first, it.second);
        }
    }

    /**
     * @brief Add an interaction and/or quadratic bias to a binary quadratic model.
     * 
     * @param u
     * @param v
     * @param bias
     */
    void add_interaction
    (
        const IndexType &u,
        const IndexType &v,
        const FloatType &bias
    )
    {
        //add labels u and v
        _add_new_label(u);
        _add_new_label(v);
        _mat(u, v) += bias;
    }

    /**
     * @brief Add interactions and/or quadratic biases to a binary quadratic model.
     * 
     * @param quadratic
     */
    void add_interactions_from
    (
        const Quadratic<IndexType, FloatType> &quadratic
    )
    {
        for(auto &it : quadratic)
        {
            add_interaction(it.first.first, it.first.second, it.second);
        }
    }

    

    /**
     * @brief Remove variable v and all its interactions from a binary quadratic model.
     * 
     * @param v
     */
    void remove_variable
    (
        const IndexType &v
    )
    {
        _delete_label(v);
    }

    /**
     * @brief Remove specified variables and all of their interactions from a binary quadratic model.
     * 
     * @param variables
     */
    void remove_variables_from
    (
        const std::vector<IndexType> &variables
    )
    {
        for(auto &it : variables)
        {
            remove_variable(it);
        }
    }

    /**
     * @brief Remove interaction of variables u, v from a binary quadratic model.
     * 
     * @param u
     * @param v
     */
    void remove_interaction
    (
        const IndexType &u,
        const IndexType &v      
    )
    {
        _mat(u, v) = 0;
        _delete_label(u, false);
        _delete_label(v, false);
    }

    /**
     * @brief Remove all specified interactions from the binary quadratic model.
     * 
     * @param interactions
     */
    void remove_interactions_from
    (
        const std::vector<std::pair<IndexType, IndexType>> &interactions
    )
    {
        for(auto &it : interactions)
        {
            remove_interaction(it.first, it.second);
        }
    }

    /**
     * @brief Add specified value to the offset of a binary quadratic model.
     * 
     * @param offset
     */
    void add_offset
    (
        const FloatType &offset
    )
    {
        m_offset += offset;
    }

    /**
     * @brief Set the binary quadratic model's offset to zero.
     */
    void remove_offset()
    {
        add_offset(-m_offset);
    }

    /**
     * @brief Multiply by the specified scalar all the biases and offset of a binary quadratic model.
     * 
     * @param scalar
     * @param ignored_variables
     * @param ignored_interactions
     * @param ignored_offset
     */
    void scale
    (
        const FloatType &scalar,
        const std::vector<IndexType> &ignored_variables = {},
        const std::vector<std::pair<IndexType, IndexType>> &ignored_interactions = {},
        const bool ignored_offset = false
    )
    {
        if(scalar == 0.0)
            throw std::runtime_error("scalar must not be zero");

        // scale
        _quadmat *= scalar;

        // revert scale of linear
        for(const auto &it : ignored_variables)
        {
            _mat(it) *= 1.0/scalar;
        }


        // revert scale of quadratic
        for(const auto &it : ignored_interactions)
        {
            _mat(it.first, it.second) *= 1.0/scalar;
        }

        // scaling offset
        if(!ignored_offset)
        {
            m_offset *= scalar;
        }
    }

    /**
     * @brief Normalizes the biases of the binary quadratic model such that they fall in the provided range(s), and adjusts the offset appropriately.
     * 
     * @param bias_range
     * @param use_quadratic_range
     * @param quadratic_range
     * @param ignored_variables
     * @param ignored_interactions
     * @param ignored_offset
     * 
     */
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
        FloatType lin_min = _min_linear();
        FloatType lin_max = _max_linear();
        FloatType quad_min = _min_quadratic();
        FloatType quad_max = _max_quadratic();

        std::vector<FloatType> v_scale =
        {
            lin_min / l_range.first,
            lin_max / l_range.second,
            quad_min / q_range.first,
            quad_max / q_range.second
        };

        FloatType inv_scale = *std::max_element(v_scale.begin(), v_scale.end());

        // scaling
        if(inv_scale != 0.0)
        {
            scale(1.0 / inv_scale, ignored_variables, ignored_interactions, ignored_offset);
        }
    }

    /**
     * @brief Fix the value of a variable and remove it from a binary quadratic model.
     * 
     * @param v
     * @param value
     */
    void fix_variable
    (
        const IndexType &v,
        const int32_t &value
    )
    {
        std::vector<std::pair<IndexType, IndexType>> interactions;
        const Quadratic<IndexType, FloatType>& quadratic = this->get_quadratic();
        for(const auto &it : quadratic)
        {
            if(it.first.first == v)
            {
                add_variable(it.first.second, value*it.second);
                interactions.push_back(it.first);
            }
            else if(it.first.second == v)
            {
                add_variable(it.first.first, value*it.second);
                interactions.push_back(it.first);
            }
        }
        remove_interactions_from(interactions);
        add_offset(_mat(v)*value);
        remove_variable(v);
    }

    /**
     * @brief Fix the value of the variables and remove it from a binary quadratic model.
     * 
     * @param fixed
     */
    void fix_variables
    (
        const std::vector<std::pair<IndexType, int32_t>> &fixed
    )
    {
        for(auto &it : fixed)
        {
            fix_variable(it.first, it.second);
        }
    }

    /**
     * @brief Flip variable v in a binary quadratic model.
     * 
     * @param v
     */
    void flip_variable
    (
        const IndexType &v
    )
    {

        if(m_vartype == Vartype::SPIN)
        {
            size_t i = _label_to_idx.at(v);
            _quadmat.row(i) *= -1;
            _quadmat.col(i) *= -1;
        }
        else if(m_vartype == Vartype::BINARY)
        {
            //change vartype to spin
            this->change_vartype(Vartype::SPIN);

            size_t i = _label_to_idx.at(v);
            _quadmat.row(i) *= -1;
            _quadmat.col(i) *= -1;


            //change vartype to binary
            this->change_vartype(Vartype::BINARY);
        }
    }

    /**
     * @brief Enforce u, v being the same variable in a binary quadratic model. (currently disabled)
     * 
     * @param u
     * @param v
     */
    //void contract_variables
    //(
    //    const IndexType &u,
    //    const IndexType &v
    //)
    //{

    //    auto p1 = std::make_pair(u, v);
    //    auto p2 = std::make_pair(v, u);

    //    const Quadratic<IndexType, FloatType>& quadratic = this->get_quadratic();

    //    if(quadratic.count(p1) != 0)
    //    {
    //        if(m_vartype == Vartype::BINARY)
    //        {
    //            add_variable(u, quadratic.at(p1));
    //        }
    //        else if(m_vartype == Vartype::SPIN)
    //        {
    //            add_offset(quadratic.at(p1));
    //        }
    //        remove_interaction(u, v);
    //    }
    //    if(quadratic.count(p2) != 0)
    //    {
    //        if(m_vartype == Vartype::BINARY)
    //        {
    //            add_variable(u, quadratic.at(p2));
    //        }
    //        else if(m_vartype == Vartype::SPIN)
    //        {
    //            add_offset(quadratic.at(p2));
    //        }
    //        remove_interaction(v, u);
    //    }

    //    std::vector<std::pair<IndexType, IndexType>> interactions;
    //    for(const auto &it : quadratic)
    //    {
    //        if(it.first.first == v)
    //        {
    //            std::cout << u << " " << it.first.second << std::endl;
    //            add_interaction(u, it.first.second, it.second);
    //            interactions.push_back(it.first);
    //        }
    //        else if(it.first.second == v)
    //        {
    //            std::cout << it.first.first << " " << u << std::endl;
    //            add_interaction(it.first.first, u, it.second);
    //            interactions.push_back(it.first);
    //        }
    //    }
    //    remove_interactions_from(interactions);

    //    add_variable(u, _mat(v));
    //    remove_variable(v);
    //}

    /* Transformations */

    /**
     * @brief Create a binary quadratic model with the specified vartype.
     * This function does not return any object.
     * 
     * @param vartype
     */
    void change_vartype
    (
        const Vartype &vartype
    )
    {
        if(m_vartype == Vartype::BINARY && vartype == Vartype::SPIN) // binary -> spin
        {
            _binary_to_spin();                
        }
        else if(m_vartype == Vartype::SPIN && vartype == Vartype::BINARY) // spin -> binary
        {
            _spin_to_binary();
        }
    }

    /**
     * @brief Create a binary quadratic model with the specified vartype.
     * This function generates and returns a new object.
     *
     * @param vartype
     * @param inplace if set true, the current object is converted.
     *
     * @return created object
     */
    BinaryQuadraticModel_Dense<IndexType, FloatType> change_vartype
    (
        const Vartype &vartype,
        bool inplace
    )
    {
        BinaryQuadraticModel_Dense<IndexType, FloatType> new_bqm = *this;
        if(inplace == true){
            this->change_vartype(vartype);
        }
        new_bqm.change_vartype(vartype);

        return new_bqm;
    }


    /* Methods */

    /**
     * @brief Determine the energy of the specified sample of a binary quadratic model.
     * 
     * @param sample
     * @return An energy with respect to the sample.
     */
    FloatType energy(const Sample<IndexType> &sample) const
    {
        FloatType en = m_offset;
        //initialize vector
        Vector s = Vector::Zero(_quadmat.rows());
        for(const auto& elem : sample){
            s[_label_to_idx.at(elem.first)] = elem.second;
        }
        s[_quadmat.rows()-1] = 1;
        
        return en + (s.transpose() * _quadmat * s) - 1;
    }
    
    /**
     * @brief Determine the energies of the given samples.
     * 
     * @param samples_like
     * @return A vector including energies with respect to the samples.
     */
    std::vector<FloatType> energies(const std::vector<Sample<IndexType>> &samples_like) const
    {
        std::vector<FloatType> en_vec;
        for(auto &it : samples_like)
        {
            en_vec.push_back(energy(it));
        }
        return en_vec;
    }
    
    /* Conversions */
    /**
     * @brief Convert a binary quadratic model to QUBO format.
     * 
     * @return A tuple including a quadratic bias and an offset.
     */
    std::tuple<Quadratic<IndexType, FloatType>, FloatType> to_qubo()
    {
        // change vartype to binary
        BinaryQuadraticModel_Dense bqm = change_vartype(Vartype::BINARY, false);

        const Linear<IndexType, FloatType>& linear = bqm.get_linear();
        Quadratic<IndexType, FloatType> Q = bqm.get_quadratic();
        FloatType offset = bqm.get_offset();
        for(const auto &it : linear)
        {
            Q[std::make_pair(it.first, it.first)] = it.second;
        }
        return std::make_tuple(Q, offset);
    }

    /**
     * @brief Create a binary quadratic model from a QUBO model.
     *
     * @param Q
     * @param offset
     *
     * @return Binary quadratic model with vartype set to `.Vartype.BINARY`.
     */
    static BinaryQuadraticModel_Dense from_qubo(const Quadratic<IndexType, FloatType>& Q, FloatType offset=0.0)
    {
        Linear<IndexType, FloatType> linear;
        Quadratic<IndexType, FloatType> quadratic;

        for(auto&& elem : Q){
            const auto& key = elem.first;
            const auto& value = elem.second;
            if(key.first == key.second){
                linear[key.first] = value;
            }
            else{
                quadratic[std::make_pair(key.first, key.second)] = value;
            }
        }

        return BinaryQuadraticModel_Dense<IndexType, FloatType>(linear, quadratic, offset, Vartype::BINARY);
    }

    /**
     * @brief Convert a binary quadratic model to Ising format.
     * 
     * @return A tuple including a linear bias, a quadratic bias and an offset.
     */
    std::tuple<Linear<IndexType, FloatType>, Quadratic<IndexType, FloatType>, FloatType> to_ising()
    {
        // change vartype to spin
        BinaryQuadraticModel_Dense bqm = change_vartype(Vartype::SPIN, false);

        const Linear<IndexType, FloatType>& linear = bqm.get_linear();
        const Quadratic<IndexType, FloatType>& quadratic = bqm.get_quadratic();
        FloatType offset = bqm.get_offset();
        return std::make_tuple(linear, quadratic, offset);
    }

    /**
     * @brief Create a binary quadratic model from an Ising problem.
     *
     * @param linear
     * @param quadratic
     * @param offset
     *
     * @return Binary quadratic model with vartype set to `.Vartype.SPIN`.
     */
    static BinaryQuadraticModel_Dense from_ising(const Linear<IndexType, FloatType>& linear, const Quadratic<IndexType, FloatType>& quadratic, FloatType offset=0.0)
    {
        return BinaryQuadraticModel_Dense<IndexType, FloatType>(linear, quadratic, offset, Vartype::SPIN);
    }


    /**
     * @brief generate interaction matrix with given list of indices
     * The generated matrix will be the following triangular matrix:
     * \f[
     * \begin{pmatrix}
     * J_{0,0} & J_{0,1} & \cdots & J_{0,N-1} & h_{0}\\
     * 0 & J_{1,1} & \cdots & J_{1,N-1} & h_{1}\\
     * \vdots & \vdots & \vdots & \vdots & \vdots \\
     * 0 & 0 & \cdots & J_{N-1,N-1} & h_{N-1}\\
     * 0 & 0 & \cdots & 0 & 1 \\
     * \end{pmatrix}
     * \f]
     *
     * @param indices
     *
     * @return corresponding interaction matrix (Eigen)
     */
    Matrix interaction_matrix() const {
        return this->_quadmat;
    }
//
//
//    using json = nlohmann::json;
//
//    /**
//     * @brief Convert the binary quadratic model to a serializable object
//     * user_bytes is assume to be set to False
//     *
//     * @return An object that can be serialized (nlohmann::json)
//     */
//    json to_serializable() const
//    {
//        std::string schema_version = "3.0.0";
//        //all variables are contained in the keys of m_linear.
//        //All we need is to traverse the keys of m_linear.
//        /*
//         * output sample
//         * >>> bqm = dimod.BinaryQuadraticModel({'c': -1, 'd': 1}, {('a', 'd'): 2, ('b', 'e'): 5, ('a', 'c'): 3}, 0.0, dimod.BINARY)
//         *
//         * >>> bqm.to_serializable()
//         * {'type': 'BinaryQuadraticModel', 'version': {'bqm_schema': '3.0.0'}, 'use_bytes': False, 'index_type': 'uint16', 'bias_type': 'float32', 'num_variables': 5, 'num_interactions': 3, 'variable_labels': ['a', 'b', 'c', 'd', 'e'], 'variable_type': 'BINARY', 'offset': 0.0, 'info': {}, 'linear_biases': [0.0, 0.0, -1.0, 1.0, 0.0], 'quadratic_biases': [3.0, 2.0, 5.0], 'quadratic_head': [0, 0, 1], 'quadratic_tail': [2, 3, 4]}
//         */
//
//        //set variables (sorted)
//        std::vector<IndexType> variables;
//        variables.reserve(m_linear.size());
//        for(auto&& elem : m_linear)
//        {
//            variables.push_back(elem.first);
//        }
//
//        std::sort(variables.begin(), variables.end());
//
//        size_t num_variables = variables.size();
//        
//        //set sorted linear biases
//        std::vector<FloatType> l_bias;
//        for(auto&& key : variables)
//        {
//            l_bias.push_back(m_linear.at(key));
//        }
//
//        //set quadratic head, tail and biases
//        std::vector<size_t> q_head, q_tail;
//        std::vector<FloatType> q_bias;
//        for(auto&& elem : m_quadratic)
//        {
//            auto it_head = std::find(variables.begin(), variables.end(), elem.first.first);
//            auto it_tail = std::find(variables.begin(), variables.end(), elem.first.second);
//            size_t idx_head = std::distance(variables.begin(), it_head);
//            size_t idx_tail = std::distance(variables.begin(), it_tail);
//            q_head.push_back(idx_head);
//            q_tail.push_back(idx_tail);
//            q_bias.push_back(elem.second);
//        }
//
//        size_t num_interactions =  m_quadratic.size();
//
//        //set index_dtype
//        std::string index_dtype = num_variables <= 65536UL ? "uint16" : "uint32";
//
//        //set bias_type
//        std::string bias_type;
//        if(typeid(m_offset) == typeid(float))
//        {
//            bias_type = "float32";
//        }
//        else if(typeid(m_offset) == typeid(double))
//        {
//            bias_type = "float64";
//        }
//        else
//        {
//            std::cerr << "FloatType must be float or double." << std::endl;
//        }
//
//        //set variable type
//        std::string variable_type;
//        if(m_vartype == Vartype::SPIN)
//        {
//            variable_type = "SPIN";
//        }
//        else if(m_vartype == Vartype::BINARY)
//        {
//            variable_type = "BINARY";
//        }
//        else
//        {
//            std::cerr << "Variable type must be SPIN or BINARY." << std::endl;
//        }
//
//        json output;
//        output["type"] = "BinaryQuadraticModel";
//        output["version"] = {{"bqm_schema", "3.0.0"}}
//        output["variable_labels"] = variables;
//        output["use_bytes"] = false;
//        output["index_type"] = index_dtype;
//        output["bias_type"] = bias_type;
//        output["num_variables"] = num_variables;
//        output["num_interactions"] = num_interactions;
//        output["variable_labels"] = variables;
//        output["variable_type"] = variable_type;
//        output["offset"] = m_offset;
//        output["info"] = m_info;
//        output["linear_biases"] = l_bias;
//        output["quadratic_biases"] = q_bias;
//        output["quadratic_head"] = q_head;
//        output["quadratic_tail"] = q_tail;
//
//        return output;
//    }
//
//    /**
//     * @brief Create a BinaryQuadraticModel instance from a serializable object.
//     * 
//     * @tparam IndexType_serial
//     * @tparam FloatType_serial
//     * @param input
//     * @return BinaryQuadraticModel<IndexType_serial, FloatType_serial> 
//     */
//    template <typename IndexType_serial = IndexType, typename FloatType_serial = FloatType>
//    static BinaryQuadraticModel_Dense<IndexType_serial, FloatType_serial> from_serializable(const json &input)
//    {
//        //extract type and version
//        std::string type = input["type"];
//        if(type != "BinaryQuadraticModel")
//        {
//            throw std::runtime_error("Type must be \"BinaryQuadraticModel\".\n");
//        }
//        std::string version = input["version"]["bqm_schema"];
//        if(version != "3.0.0")
//        {
//            throw std::runtime_error("bqm_schema must be 3.0.0.\n");
//        }
//
//        //extract variable_type
//        Vartype vartype;
//        std::string variable_type = input["variable_type"];
//        if(variable_type == "SPIN")
//        {
//            vartype = Vartype::SPIN;
//        }
//        else if(variable_type == "BINARY")
//        {
//            vartype = Vartype::BINARY;
//        }
//        else
//        {
//            throw std::runtime_error("variable_type must be SPIN or BINARY.");
//        }
//
//        //extract linear biases
//        Linear<IndexType_serial, FloatType_serial> linear;
//        std::vector<IndexType_serial> variables = input["variable_labels"];
//        std::vector<FloatType_serial> l_bias = input["linear_biases"];
//        for(size_t i = 0; i < variables.size(); ++i)
//        {
//            insert_or_assign(linear, variables[i], l_bias[i]);
//        }
//
//        //extract quadratic biases
//        Quadratic<IndexType_serial, FloatType_serial> quadratic;
//        std::vector<size_t> q_head = input["quadratic_head"];
//        std::vector<size_t> q_tail = input["quadratic_tail"];
//        std::vector<FloatType_serial> q_bias = input["quadratic_biases"];
//        for(size_t i = 0; i < q_head.size(); ++i)
//        {
//            insert_or_assign(quadratic, std::make_pair(variables[q_head[i]], variables[q_tail[i]]), q_bias[i]);
//        }
//
//        //extract offset and info
//        FloatType_serial offset = input["offset"];
//        std::string info = (input["info"].empty())?"":input["info"];
//
//        //construct a BinaryQuadraticModel instance
//        BinaryQuadraticModel_Dense<IndexType_serial, FloatType_serial> bqm(linear, quadratic, offset, vartype, info);
//        return bqm;
//    }

};
}
#endif
