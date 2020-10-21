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

#include "impl.hpp"
#include <string>

namespace py = pybind11;

using namespace py::literals;
using namespace cimod;

PYBIND11_MODULE(cxxcimod, m){
    m.doc() = "C++ library for a binary quadratic model";

    /**********************************************************
    //BinaryQuadraticModel
     **********************************************************/

    py::enum_<Vartype>(m, "Vartype")
        .value("SPIN", Vartype::SPIN)
        .value("BINARY", Vartype::BINARY)
        .value("NONE", Vartype::NONE)
        .export_values();

    declare_BQM<size_t, double>(m, "BinaryQuadraticModel");
    declare_BQM<std::string, double>(m, "BinaryQuadraticModel_str");
    declare_BQM<std::tuple<size_t, size_t>, double>(m, "BinaryQuadraticModel_tuple2");
    declare_BQM<std::tuple<size_t, size_t, size_t>, double>(m, "BinaryQuadraticModel_tuple3");
    declare_BQM<std::tuple<size_t, size_t, size_t, size_t>, double>(m, "BinaryQuadraticModel_tuple4");

}





























