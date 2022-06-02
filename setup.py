from skbuild import setup

setup(
    cmake_source_dir = "cimod/", 
    cmake_with_sdist = True, 
    cmake_languages = ('C', 'CXX',), 
    cmake_minimum_required_version = '2.20', 
    
)
