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

#include "main.hpp"

#include <string>

namespace py = pybind11;

using namespace py::literals;
using namespace cimod;

PYBIND11_MODULE( cxxcimod, m ) {
  m.doc() = "C++ library for a binary quadratic model";

  /**********************************************************
  //BinaryQuadraticModel and BinaryPolynomialModel
   **********************************************************/

  py::enum_<Vartype>( m, "Vartype" )
      .value( "SPIN", Vartype::SPIN )
      .value( "BINARY", Vartype::BINARY )
      .value( "NONE", Vartype::NONE )
      .export_values();

  declare_BQM<int64_t, double, cimod::Dense>( m, "BinaryQuadraticModel_Dense" );
  declare_BQM<std::string, double, cimod::Dense>( m, "BinaryQuadraticModel_str_Dense" );
  declare_BQM<std::tuple<size_t, size_t>, double, cimod::Dense>( m, "BinaryQuadraticModel_tuple2_Dense" );
  declare_BQM<std::tuple<size_t, size_t, size_t>, double, cimod::Dense>( m, "BinaryQuadraticModel_tuple3_Dense" );
  declare_BQM<std::tuple<size_t, size_t, size_t, size_t>, double, cimod::Dense>( m, "BinaryQuadraticModel_tuple4_Dense" );

  declare_BQM<int64_t, double, cimod::Sparse>( m, "BinaryQuadraticModel_Sparse" );
  declare_BQM<std::string, double, cimod::Sparse>( m, "BinaryQuadraticModel_str_Sparse" );
  declare_BQM<std::tuple<size_t, size_t>, double, cimod::Sparse>( m, "BinaryQuadraticModel_tuple2_Sparse" );
  declare_BQM<std::tuple<size_t, size_t, size_t>, double, cimod::Sparse>( m, "BinaryQuadraticModel_tuple3_Sparse" );
  declare_BQM<std::tuple<size_t, size_t, size_t, size_t>, double, cimod::Sparse>( m, "BinaryQuadraticModel_tuple4_Sparse" );

  declare_BQM<int64_t, double, cimod::Dict>( m, "BinaryQuadraticModel_Dict" );
  declare_BQM<std::string, double, cimod::Dict>( m, "BinaryQuadraticModel_str_Dict" );
  declare_BQM<std::tuple<size_t, size_t>, double, cimod::Dict>( m, "BinaryQuadraticModel_tuple2_Dict" );
  declare_BQM<std::tuple<size_t, size_t, size_t>, double, cimod::Dict>( m, "BinaryQuadraticModel_tuple3_Dict" );
  declare_BQM<std::tuple<size_t, size_t, size_t, size_t>, double, cimod::Dict>( m, "BinaryQuadraticModel_tuple4_Dict" );

  declare_BPM<int64_t, double>( m, "BinaryPolynomialModel" );
  declare_BPM<std::string, double>( m, "BinaryPolynomialModel_str" );
  declare_BPM<std::tuple<int64_t, int64_t>, double>( m, "BinaryPolynomialModel_tuple2" );
  declare_BPM<std::tuple<int64_t, int64_t, int64_t>, double>( m, "BinaryPolynomialModel_tuple3" );
  declare_BPM<std::tuple<int64_t, int64_t, int64_t, int64_t>, double>( m, "BinaryPolynomialModel_tuple4" );
}
