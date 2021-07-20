# cimod : C++ header-only library for a binary quadratic model 

* [Documents](https://openjij.github.io/cimod/html/index.html)

# How to use

You should only include a header `src/binary_quadratic_model.hpp` in your project.

## Example

### C++

```cpp
#include "src/binary_quadratic_model.hpp"

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
BinaryQuadraticModel<uint32_t, double, cimod::Dense> bqm(linear, quadratic, offset, vartype);

//linear terms -> bqm.get_linear()
//quadratic terms -> bqm.get_quadratic()

return 0;
}
```

### Python

```python
import cimod
import dimod

# Set linear biases and quadratic biases
linear = {1:1.0, 2:2.0, 3:3.0, 4:4.0}
quadratic = {(1,2):12.0, (1,3):13.0, (1,4):14.0, (2,3):23.0, (2,4):24.0, (3,4):34.0}

# Set offset
offset = 0.0

# Set variable type
vartype = dimod.BINARY

# Create a BinaryQuadraticModel instance
bqm = cimod.BinaryQuadraticModel(linear, quadratic, offset, vartype)

print(bqm.linear)
print(bqm.quadratic)

```

## Install 

### via this directory

```sh
$ python -m pip install .
```

### via pip

```sh
$ pip install jij-cimod 
```

## Benchmark

### Benchmark code

```python
import dimod
import cimod
import time

fil = open("benchmark", "w")
fil.write("N t_dimod t_cimod\n")

def benchmark(N, test_fw):
    linear = {}
    quadratic = {}

    spin = {}

    # interactions

    for i in range(N):
        spin[i] = 1

    for elem in range(N):
        linear[elem] = 2.0*elem;

    for i in range(N):
        for j in range(i+1, N):
            if i != j:
                quadratic[(i,j)] = (i+j)/(N)

    t1 = time.time()

    # initialize
    a = test_fw.BinaryQuadraticModel(linear, quadratic, 0, test_fw.BINARY)
    a.change_vartype(test_fw.SPIN)

    # calculate energy for 50 times.
    for _ in range(50):
        print(a.energy(spin))

    t2 = time.time()

    return t2-t1

d_arr = []
c_arr = []

for N in [25, 50, 100, 200, 300, 400, 600, 800,1000, 1600, 2000, 3200, 5000]:
    print("N {}".format(N))
    d = benchmark(N, dimod)
    c = benchmark(N, cimod)
    print("{} {} {}".format(N, d, c))
    fil.write("{} {} {}\n".format(N, d, c))
```

### Software versions

|       Package        | Version |
| -------------------- | ------- |
| [cimod](https://github.com/OpenJij/cimod)     | 1.0.3   |
| [dimod](https://github.com/dwavesystems/dimod)| 0.9.2   |

### Result
![benchmark](https://github.com/OpenJij/cimod/blob/image_store/figure.png)



