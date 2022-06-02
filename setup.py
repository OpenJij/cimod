import sys

from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

setup_requires=[ 
  'pybind11', 
]

# Require pytest-runner only when running tests.
if any(arg in sys.argv for arg in ('pytest', 'test')):
	setup_requires.append('pytest-runner')

# Add CMake as a build requirement if cmake is not installed or is too low a version.
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion('3.20'):
        setup_requires.append('cmake')
except SKBuildError:
    setup_requires.append('cmake')  

setup(
	setup_requires=setup_requires,
	packages=[  
 		'cimod', 
 		'cimod.model', 
 		'cimod.model.legacy',
 		'cimod.utils',
         ],
	include_package_data=True,
	zip_safe=False,
)
