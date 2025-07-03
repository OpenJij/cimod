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
[![codecov](https://codecov.io/gh/OpenJij/cimod/branch/master/graph/badge.svg?token=BE45W9FJHA)](https://codecov.io/gh/OpenJij/cimod)


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

This project uses [uv](https://docs.astral.sh/uv/) for dependency management.

```sh
# Clone repository
$ git clone https://github.com/OpenJij/cimod.git
$ cd cimod

# Install uv (choose one method)
$ curl -LsSf https://astral.sh/uv/install.sh | sh  # macOS/Linux
# or: brew install uv                              # Homebrew
# or: pip install uv                               # fallback option

# Install with development dependencies (exact versions)
$ uv sync --locked --group dev

# Verify installation
$ uv run python -c "import cimod; print('cimod installed successfully')"
$ uv run pytest tests/ -v --tb=short
```

## Development

### Dependency Groups

The project uses [PEP 735](https://peps.python.org/pep-0735/) dependency groups in `pyproject.toml`:

| Group | Purpose | Command |
|-------|---------|---------|
| **dev** | Development environment (build + test + format) | `uv sync --group dev` |
| **test** | Testing tools (pytest, coverage) | `uv sync --group test` |
| **docs** | Documentation generation | `uv sync --group docs` |
| **format** | Code formatting (ruff only) | `uv sync --group format` |
| **all** | Complete environment (dev + docs) | `uv sync --group all` |

**Lock file usage**:
- **For exact reproduction** (CI/CD, verification): `uv sync --locked --group dev`
- **For development** (may update dependencies): `uv sync --group dev`
- **To update lock file**: `uv lock` or `uv lock --upgrade`

Dependencies are locked in `uv.lock` for reproducible builds across environments.

### System Requirements
- **Python**: 3.9-3.13
- **C++**: C++17 compatible compiler  
- **CMake**: 3.20+ (for C++ development)

## Testing

### Python Tests

```sh
# Install test dependencies (exact versions)
$ uv sync --locked --group test

# Basic test run
$ uv run pytest tests/ -v

# With coverage report
$ uv run pytest tests/ -v --cov=cimod --cov-report=html
$ uv run python -m coverage html
```

### C++ Tests

```sh
# Build C++ tests (independent of Python environment)
$ mkdir build 
$ cmake -DCMAKE_BUILD_TYPE=Debug -S . -B build
$ cmake --build build --parallel

# Run C++ tests
$ cd build
$ ./tests/cimod_test
# Alternatively Use CTest 
$ ctest --extra-verbose --parallel --schedule-random
```

**Requirements**: CMake > 3.22, C++17

## Code Quality

### Unified Tooling with Ruff

```sh
# Install format dependencies (exact versions)
$ uv sync --locked --group format

# Check and fix all issues
$ uv run ruff check .              # Lint check
$ uv run ruff format .             # Format code  
$ uv run ruff check . --fix        # Auto-fix issues

# All-in-one check (recommended)
$ uv run ruff check . && uv run ruff format --check .
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

## Notes
* As explained in https://github.com/OpenJij/cimod/issues/48, specifying self-loop index (e.g. `{(2, 2): 5}`) in the `quadratic` argument in `BinaryQuadraticModel` is not allowed.

### Licences

Copyright 2025 Jij Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

```
 http://www.apache.org/licenses/LICENSE-2.0  
```

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.