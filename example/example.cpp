#include "../src/binary_quadratic_model.hpp"

using namespace cimod;

int main()
{
    // Set linear biases and quadratic biases
    Linear<uint32_t, double> linear{ {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0} };
    Quadratic<uint32_t, double> quadratic
    {
        {std::make_pair(1, 2), 12.0}, {std::make_pair(1, 3), 13.0}, {std::make_pair(1, 4), 14.0},
        {std::make_pair(2, 3), 23.0}, {std::make_pair(2, 4), 24.0},
        {std::make_pair(3, 4), 34.0}
    };

    // Set offset
    double offset = 0.0;

    // Set variable type
    Vartype vartype = Vartype::BINARY;

    // Create a BinaryQuadraticModel instance
    BinaryQuadraticModel<uint32_t, double> bqm(linear, quadratic, offset, vartype);

    // Print informations of bqm
    bqm.print();

    return 0;
}