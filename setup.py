from skbuild.cmaker import CMaker
from skbuild import setup

setup(
  setup_requires = [ 
    "setuptools_scm[toml]",
    "pytest-runner",
    "cmake",
    "pybind11",
    "scikit-build",
    "setuptools",
  ],
  packages = [ 
    "cimod", 
    "cimod.model", 
    "cimod.model.legacy", 
    "cimod.utils",
  ],
  zip_safe = False,
)
