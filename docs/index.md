# cimod : C++ header-only library for a binary quadratic model

[![PyPI version shields.io](https://img.shields.io/pypi/v/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI implementation](https://img.shields.io/pypi/implementation/jij-cimod.svg)](https://pypi.python.org/pypi/ji-cimod/)
[![PyPI format](https://img.shields.io/pypi/format/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI license](https://img.shields.io/pypi/l/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI download month](https://img.shields.io/pypi/dm/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)

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

### For Users

```sh
# Binary package (recommended)
$ pip install jij-cimod

# From source  
$ pip install --no-binary=jij-cimod jij-cimod 

# Latest development version
$ pip install git+https://github.com/OpenJij/cimod.git
```

### For Developers

```sh
# Clone repository
$ git clone https://github.com/OpenJij/cimod.git
$ cd cimod

# Create virtual environment
$ python -m venv .venv
$ source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Install with development dependencies
$ pip install -e ".[dev]"

# Verify installation
$ python -c "import cimod; print('cimod installed successfully')"
$ pytest tests/ -v --tb=short
```

## Development

### Dependency Groups

The project uses [PEP 621](https://peps.python.org/pep-0621/) optional dependencies in `pyproject.toml`:

| Group | Purpose | Install Command |
|-------|---------|-----------------|
| **dev** | Complete development environment | `pip install -e ".[dev]"` |
| **test** | Testing tools (pytest, coverage) | `pip install -e ".[test]"` |
| **docs** | Documentation generation | `pip install -e ".[docs]"` |
| **format** | Code formatting (ruff only) | `pip install -e ".[format]"` |

**Note**: The `dev` group excludes `docs` dependencies for faster installation and C++ build avoidance.  
For full functionality: `pip install -e ".[dev,docs]"`

### System Requirements
- **Python**: 3.9-3.13
- **C++**: C++17 compatible compiler  
- **CMake**: 3.20+ (for C++ development)

## Testing

### Python

```sh
$ python -m venv .venv
$ source .venv/bin/activate
$ pip install -e ".[test]"
$ pytest tests/ -v --cov=cimod --cov-report=html 
$ python -m coverage html
```

### C++

```sh
$ mkdir build 
$ cmake -DCMAKE_BUILD_TYPE=Debug -S . -B build
$ cmake --build build --parallel
$ cd build
$ ./tests/cimod_test
# Alternatively Use CTest 
$ ctest --extra-verbose --parallel --schedule-random
```

**Requirements**: CMake > 3.22, C++17

## Code Quality

### Unified Tooling with Ruff

```sh
# Install development dependencies (includes ruff)
$ pip install -e ".[dev]"

# Check and fix all issues
$ ruff check .              # Lint check
$ ruff format .             # Format code  
$ ruff check . --fix        # Auto-fix issues

# All-in-one check (recommended)
$ ruff check . && ruff format --check .
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
