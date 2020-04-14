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


#include <binary_quadratic_model.hpp>
#include <pybind11/pybind11.h>
#include <type_traits>
#include <pybind11_json/pybind11_json.hpp>

namespace py = pybind11;

using namespace py::literals;
using namespace cimod;

PYBIND11_MODULE(cimod, m){
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

    py::class_<BinaryQuadraticModel<IndexType, FloatType>>(m, "BinaryQuadraticModel")
        .def(py::init<Linear<IndexType, FloatType>, Quadratic<IndexType, FloatType>, FloatType, Var>())

}





























