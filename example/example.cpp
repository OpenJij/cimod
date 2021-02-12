#include "../src/binary_quadratic_model.hpp"
#include "../src/binary_polynomial_model.hpp"

using namespace cimod;

int main() {
   
   
   // Set linear biases and quadratic biases
   Polynomial<uint32_t, double> polynomial {
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

   
   // Set variable type
   Vartype vartype = Vartype::BINARY;
   
   //Test function of energy
   Sample<uint32_t> sample = {{1, 1}, {2, 0}, {3, 1}, {4, 0}};
   
   // Create a BinaryPolynominalModel instance
   BinaryPolynomialModel<uint32_t, double> bpm(polynomial, vartype);
   
   //Test remove_variable
   //bpm.remove_variable(3);
   
   //bpm.normalize();
   
   // Print informations of bpm
   bpm.print();
   
   printf("energy(polynomial): %lf\n", bpm.energy(sample));
   
   // Set linear biases and quadratic biases
   Linear<std::string, double> linear{ {"4", 4.0}, {"3", 3.0}, {"2", 2.0}, {"1", 1.0} };
   Quadratic<std::string, double> quadratic
   {
      {std::make_pair("4", "3"), 34.0},
      {std::make_pair("3", "2"), 23.0}, {std::make_pair("4", "2"), 24.0},
       {std::make_pair("1", "2"), 12.0}, {std::make_pair("1", "3"), 13.0}, {std::make_pair("1", "4"), 14.0}
   };
   
   // Set offset
   double offset = 0.0;
   
   // Create a BinaryQuadraticModel instance
   BinaryQuadraticModel<std::string, double> bqm(linear, quadratic, offset, vartype);
   bqm.print();
   //printf("energy=%lf\n", bqm.energy(sample));
   printf("/*****************************/\n");

   //Test remove_variable
   //bqm.remove_variable(3);
   
   //Test fix_variables
   //bqm.fix_variable(1, 2.0);
   
   //Test flip_variable
   //bqm.flip_variable(4);
   
   //Test contract_variables
   //bqm.contract_variables(1, 2);
   
   //bqm.normalize();
   
   // Print informations of bqm
   //bqm.print();
   
   
   //printf("energy(quadratic): %lf\n", bqm.energy(sample));
   printf("off_set: %lf\n", bqm.get_offset());
   return 0;
}
