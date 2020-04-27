//    Copyright 2019 Jij Inc.

//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at

//        http://www.apache.org/licenses/LICENSE-2.0

//    Unless required by applicable law or agreed to in writing, software 
//    distributed under the License is distributed on an "AS IS" BASIS, 
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#include <pybind11_json/pybind11_json.hpp>
#include <nlohmann/json.hpp>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <binary_quadratic_model.hpp>

namespace py = pybind11;

using namespace py::literals;
using namespace cimod;

PYBIND11_MODULE(cxxcimod, m){
    m.doc() = "C++ library for a binary quadratic model";

    /**********************************************************
    //BinaryQuadraticModel
     **********************************************************/

    //By default indices of cimod.BinaryQuadraticModel are integers.
    using IndexType = size_t;
    using FloatType = double;

    py::enum_<Vartype>(m, "Vartype")
        .value("SPIN", Vartype::SPIN)
        .value("BINARY", Vartype::BINARY)
        .value("NONE", Vartype::NONE)
        .export_values();

    using BQM = BinaryQuadraticModel<IndexType, FloatType>;

    py::class_<BQM>(m, "BinaryQuadraticModel")
        .def(py::init<Linear<IndexType, FloatType>, Quadratic<IndexType, FloatType>, FloatType, Vartype, std::string>(), "linear"_a, "quadratic"_a, "offset"_a, "vartype"_a, "info"_a="")
        .def("length", &BQM::length)
        .def("contains", &BQM::contains, "v"_a)
        .def("get_linear", &BQM::get_linear)
        .def("get_quadratic", &BQM::get_quadratic)
        .def("get_adjacency", &BQM::get_adjacency)
        .def("get_offset", &BQM::get_offset)
        .def("get_vartype", &BQM::get_vartype)
        .def("get_info", &BQM::get_info)
        .def("print", &BQM::print)
        .def("empty", &BQM::empty)
        .def("add_variable", &BQM::add_variable, "v"_a, "bias"_a, "vartype"_a=Vartype::NONE)
        .def("add_variables_from", &BQM::add_variables_from, "linear"_a, "vartype"_a)
        .def("add_interaction", &BQM::add_interaction, "u"_a, "v"_a, "bias"_a, "vartype"_a=Vartype::NONE)
        .def("add_interactions_from", &BQM::add_interactions_from, "quadratic"_a, "vartype"_a)
        .def("remove_variable", &BQM::remove_variable, "v"_a)
        .def("remove_variables_from", &BQM::remove_variables_from, "variables"_a)
        .def("remove_interaction", &BQM::remove_interaction, "u"_a, "v"_a)
        .def("remove_interactions_from", &BQM::remove_interactions_from, "interactions"_a)
        .def("add_offset", &BQM::add_offset, "offset"_a)
        .def("remove_offset", &BQM::remove_offset)
        .def("scale", &BQM::scale, "scalar"_a, "ignored_variables"_a=std::vector<IndexType>(), "ignored_interactions"_a=std::vector<std::pair<IndexType, IndexType>>(), "ignored_offset"_a=false)
        .def("normalize", &BQM::normalize, "bias_range"_a=std::pair<FloatType, FloatType>(1.0, 1.0), "use_quadratic_range"_a=false, "quadratic_range"_a=std::pair<FloatType, FloatType>(1.0, 1.0), "ignored_variables"_a=std::vector<IndexType>(), "ignored_interactions"_a=std::vector<std::pair<IndexType, IndexType>>(), "ignored_offset"_a=false)
        .def("fix_variable", &BQM::fix_variable, "v"_a, "value"_a)
        .def("fix_variables", &BQM::fix_variables, "fixed"_a)
        .def("flip_variable", &BQM::flip_variable, "v"_a)
        .def("update", &BQM::update, "bqm"_a, "ignore_info"_a=true)
        .def("contract_variables", &BQM::contract_variables, "u"_a, "v"_a)
        .def("change_vartype", &BQM::change_vartype, "vartype"_a, "implace"_a=true)
        .def_static("spin_to_binary", &BQM::spin_to_binary, "linear"_a, "quadratic"_a, "offset"_a)
        .def_static("binary_to_spin", &BQM::binary_to_spin, "linear"_a, "quadratic"_a, "offset"_a)
        .def("energy", &BQM::energy, "sample"_a)
        .def("energies", &BQM::energy, "samples_like"_a)
        .def("to_qubo", &BQM::to_qubo)
        .def("to_ising", &BQM::to_ising)
        .def_static("from_qubo", &BQM::from_qubo, "Q"_a, "offset"_a=0.0)
        .def_static("from_ising", &BQM::from_ising, "h"_a, "J"_a, "offset"_a=0.0)
        .def("interaction_matrix", &BQM::interaction_matrix, "indices"_a)
        //.def("to_serialiable", &BQM::to_serializable)
        //.def_static("from_serialiable", &BQM::from_serializable, "input"_a);
        .def("to_serializable", [](const BQM& self){return static_cast<py::object>(self.to_serializable());})
        .def_static("from_serializable", [](const py::object& input){return BQM::from_serializable(static_cast<nlohmann::json>(input));}, "input"_a)
        .def_static("_generate_indices_dict", [](const py::dict& linear, const py::dict& quadratic){
                py::set index_set;

                for(auto it = linear.begin(); it != linear.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    index_set.add(kv[0]);
                }

                for(auto it = quadratic.begin(); it != quadratic.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    py::tuple k = kv[0];
                    index_set.add(k[0]);
                    index_set.add(k[1]);
                }

                py::list indices;
                for(auto it = index_set.begin(); it != index_set.end(); ++it){
                    indices.append(*it);
                }

                py::dict num_to_ind;
                py::dict ind_to_num;

                for(size_t i=0; i<indices.size(); i++){
                    num_to_ind[py::int_(i)] = indices[i];
                    ind_to_num[indices[i]] = i;
                }

                return py::make_tuple(indices, num_to_ind, ind_to_num);
        }, "linear"_a, "quadratic"_a)
        .def_static("_Q_to_linear_quadratic", [](const py::dict& Q){
                py::dict linear;
                py::dict quadratic;
                for(auto it = Q.begin(); it != Q.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    py::tuple uv = kv[0];
                    py::object bias = kv[1];
                    py::object u = uv[0];
                    py::object v = uv[1];

                    if(u.is(v)){
                        //linear
                        linear[uv[0]] = bias;
                    }
                    else{
                        quadratic[uv] = bias;
                    }
                }

                return py::make_tuple(linear, quadratic);
                })
        .def("_conv_linear", [](const BQM& self, const py::dict& dic, const py::dict& conv_dict){
                py::dict ret;
                for(auto it = dic.begin(); it != dic.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    ret[conv_dict[kv[0]]] = kv[1];
                }

                return ret;
                }, "dic"_a, "conv_dict"_a)
        .def("_conv_quadratic", [](const BQM& self, const py::dict& dic, const py::dict& conv_dict){
                py::dict ret;
                for(auto it = dic.begin(); it != dic.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    py::tuple k = kv[0];
                    ret[py::make_tuple(conv_dict[k[0]], conv_dict[k[1]])] = kv[1];
                }

                return ret;
                }, "dic"_a, "conv_dict"_a)
        .def("_conv_adjacency", [](const BQM& self, const py::dict& dic, const py::dict& conv_dict){
                py::dict ret;
                for(auto it = dic.begin(); it != dic.end(); ++it){
                    py::tuple kv = py::cast(*it);
                    py::object index = kv[0];
                    py::dict adj_dic = kv[1];

                    py::dict converted_adj;
                    for(auto adj_it = adj_dic.begin(); adj_it != adj_dic.end(); ++adj_it){
                        py::tuple adj_kv = py::cast(*adj_it);
                        converted_adj[conv_dict[adj_kv[0]]] = adj_kv[1];
                    }

                    ret[conv_dict[index]] = converted_adj;
                }

                return ret;
                }, "dic"_a, "conv_dict"_a);

}





























