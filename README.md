# cimod : C++ header-only library for a binary quadratic model

[![PyPI version shields.io](https://img.shields.io/pypi/v/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI implementation](https://img.shields.io/pypi/implementation/jij-cimod.svg)](https://pypi.python.org/pypi/ji-cimod/)
[![PyPI format](https://img.shields.io/pypi/format/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI license](https://img.shields.io/pypi/l/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![PyPI download month](https://img.shields.io/pypi/dm/jij-cimod.svg)](https://pypi.python.org/pypi/jij-cimod/)
[![Downloads](https://pepy.tech/badge/jij-cimod)](https://pepy.tech/project/jij-cimod)

[![Test](https://github.com/OpenJij/cimod/actions/workflows/ci-test.yml/badge.svg)](https://github.com/OpenJij/cimod/actions/workflows/ci-test.yml)
[![Build&Upload](https://github.com/OpenJij/cimod/actions/workflows/build_and_upload.yaml/badge.svg)](https://github.com/OpenJij/cimod/actions/workflows/build_and_upload.yaml)
[![CodeQL](https://github.com/OpenJij/cimod/actions/workflows/codeql.yml/badge.svg)](https://github.com/OpenJij/cimod/actions/workflows/codeql.yml)
[![Build Documentation](https://github.com/OpenJij/cimod/actions/workflows/buid-doc.yml/badge.svg)](https://github.com/OpenJij/cimod/actions/workflows/buid-doc.yml)
[![pages-build-deployment](https://github.com/OpenJij/cimod/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/OpenJij/cimod/actions/workflows/pages/pages-build-deployment)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/55990ff022864098a2413c0cc4ab8299)](https://www.codacy.com/gh/OpenJij/cimod/dashboard?utm_source=github.com&utm_medium=referral&utm_content=OpenJij/cimod&utm_campaign=Badge_Grade)
[![Maintainability](https://api.codeclimate.com/v1/badges/59876c82cc2200ef1dfa/maintainability)](https://codeclimate.com/github/OpenJij/cimod/maintainability)
[![codecov](https://codecov.io/gh/OpenJij/cimod/branch/master/graph/badge.svg?token=BE45W9FJHA)](https://codecov.io/gh/OpenJij/cimod)

## Coverage Graph

| **Sunburst**                                                                                                                                                     | **Grid**                                                                                                                                                     | **Icicle**                                                                                                                                                     |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <a href="https://codecov.io/gh/OpenJij/cimod"><img src="https://codecov.io/gh/OpenJij/cimod/branch/main/graphs/sunburst.svg?token=BE45W9FJHA" width="100%"/></a> | <a href="https://codecov.io/gh/OpenJij/cimod"><img src="https://codecov.io/gh/OpenJij/cimod/branch/main/graphs/tree.svg?token=BE45W9FJHA" width="100%"/></a> | <a href="https://codecov.io/gh/OpenJij/cimod"><img src="https://codecov.io/gh/OpenJij/cimod/branch/main/graphs/icicle.svg?token=BE45W9FJHA" width="100%"/></a> |

- [Documents](https://openjij.github.io/Cimod-Documentation/)
- [Python Documents](https://openjij.github.io/cimod/)

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

## For Contributor

Use `pre-commit` for auto chech before git commit.
`.pre-commit-config.yaml`

```
# pipx install pre-commit 
# or 
# pip install pre-commit
pre-commit install
```

## Install

### via this directory

```sh
$ python -m pip install -vvv .
```

### via pip

```sh
# Binary
$ pip install jij-cimod
# From Source 
$ pip install --no-binary=jij-cimod jij-cimod 
```

## Test

### Python

```sh
$ python -m venv .venv
$ pip install pip-tools 
$ pip-compile setup.cfg
$ pip-compile dev-requirements.in
$ pip-sync requirements.txt dev-requirements.txt
$ source .venv/bin/activate
$ export CMAKE_BUILD_TYPE=Debug
$ python setup.py --force-cmake install --build-type Debug -G Ninja
$ python setup.py --build-type Debug test 
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

Needs: CMake > 3.22, C++17

- Format

```sh
$ pip-compile format-requirements.in
$ pip-sync format-requirements.txt
```

```sh
$ python -m isort 
$ python -m black 
```

- Aggressive Format

```sh
$ python -m isort --force-single-line-imports --verbose ./cimod
$ python -m autoflake --in-place --recursive --remove-all-unused-imports --ignore-init-module-imports --remove-unused-variables ./cimod
$ python -m autopep8 --in-place --aggressive --aggressive  --recursive ./cimod
$ python -m isort ./cimod
$ python -m black ./cimod
```

- Lint

```sh
$ pip-compile setup.cfg
$ pip-compile dev-requirements.in
$ pip-compile lint-requirements.in
$ pip-sync requirements.txt dev-requirements.txt lint-requirements.txt
```

```sh
$ python -m flake8
$ python -m mypy
$ python -m pyright
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

| Package                                        | Version |
| ---------------------------------------------- | ------- |
| [cimod](https://github.com/OpenJij/cimod)      | 1.0.3   |
| [dimod](https://github.com/dwavesystems/dimod) | 0.9.2   |

### Result

![benchmark](https://github.com/OpenJij/cimod/blob/image_store/figure.png)

### Licences

Copyright 2022 Jij Inc.

Licensed under the Apache License, Version 2.0 (the "License");\
you may not use this file except in compliance with the License.\
You may obtain a copy of the License at

```
 http://www.apache.org/licenses/LICENSE-2.0  
```

Unless required by applicable law or agreed to in writing, software\
distributed under the License is distributed on an "AS IS" BASIS,\
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\
See the License for the specific language governing permissions and\
limitations under the License.
