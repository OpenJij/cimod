#include "../src/binary_quadratic_model.hpp"
#include "../src/binary_polynomial_model.hpp"

using namespace cimod;

int main() {
   
   // Set linear biases and quadratic biases
   Polynominal<uint32_t, double> polynominal {
      // linear biases
      {{1}, 1.0}    , {{2}, 2.0}    , {{3}, 3.0},
      // quadratic biases
      {{1, 2}, 12.0}, {{1, 3}, 13.0}, {{1, 4}, 14.0},
      {{2, 3}, 23.0}, {{2, 4}, 24.0},
      {{3, 4}, 34.0},
      {{1, 2, 4}   , 124.0},
      {{1, 2, 2}   , 122.0}, // This is not allowed
      {{1, 2, 3, 4}, 1234.0}
   };
   
   // Set offset
   double offset = 0.0;
   
   // Set variable type
   Vartype vartype = Vartype::BINARY;
   
   // Create a BinaryPolynominalModel instance
   BinaryPolynomialModel<uint32_t, double> bpm(polynominal, offset, vartype);
   
   //Test remove_variable
   bpm.remove_variable(3);
   
   // Print informations of bpm
   bpm.print();
   
   // Set linear biases and quadratic biases
   Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0} /*, {4, 4.0}*/ };
   
   Quadratic<uint32_t, double> quadratic {
      {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
      {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
      {std::make_pair(3, 4), 34.0},
      {std::make_pair(2, 2), 22.0} // This is not allowed
   };
   
   // Create a BinaryQuadraticModel instance
   BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);
   
   //Test remove_variable
   bqm.remove_variable(3);
   
   //Test fix_variables
   //bqm.fix_variable(1, 2);
   
   // Print informations of bqm
   bqm.print();
   
   //Test function of energy
   Sample<uint32_t> sample = {{1, 1}, {2, 0}, {3, 1}, {4, 0}};
   printf("energy: %lf\n", bqm.energy(sample));
   
   return 0;
}
