try:
    import typing 
except ImportError:
    from typing_extensions import *
from cimod.__version import __version__
from cimod.vartype import SPIN, BINARY, Vartype
from cimod.model.binary_quadratic_model import make_BinaryQuadraticModel, make_BinaryQuadraticModel_from_JSON, BinaryQuadraticModel
from cimod.model.binary_polynomial_model import make_BinaryPolynomialModel, make_BinaryPolynomialModel_from_JSON, BinaryPolynomialModel
