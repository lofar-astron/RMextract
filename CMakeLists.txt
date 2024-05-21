cmake_minimum_required(VERSION 3.18)
project(RMextract C CXX Fortran)

find_package(PythonExtensions REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(F2PY REQUIRED)
find_package(Boost COMPONENTS python numpy REQUIRED)

# if(SKBUILD)
#   message(STATUS "The project is built using scikit-build")
#   set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIR} ${Python_NumPy_INCLUDE_DIRS})
# else()
#   set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
# endif()

# -- Python module: RMextract.EMM._EMM_Model -- ##
add_library(_EMM_Model MODULE
    RMextract/EMM/EMM_Model.cc
    RMextract/EMM/EMM_Model_wrap.cc
    RMextract/EMM/GeomagnetismLibrary.c)
python_extension_module(_EMM_Model)
install(TARGETS _EMM_Model LIBRARY DESTINATION RMextract/EMM)

# -- Python module: RMextract.pyiri._iri -- ##
add_library(_iri MODULE
    RMextract/pyiri/iri.pyf
    RMextract/pyiri/cira.for
    RMextract/pyiri/igrf.for
    RMextract/pyiri/iridreg.for
    RMextract/pyiri/iriflip.for
    RMextract/pyiri/irifun.for
    RMextract/pyiri/irisub.for
    RMextract/pyiri/iritec.for
    RMextract/pyiri/iriget.for)
python_extension_module(_iri)
install(TARGETS _iri LIBRARY DESTINATION RMextract/pyiri)

# # -- Python module: RMextract.pyiriplas._iriplas -- ##
add_library(_iriplas MODULE
    RMextract/pyiriplas/Iris2017.for
    RMextract/pyiriplas/igrf.for
    RMextract/pyiriplas/indx2017.for
    RMextract/pyiriplas/irif2019.for
    RMextract/pyiriplas/iriplas.pyf
    RMextract/pyiriplas/iriplas_main.for)
python_extension_module(_iriplas)
install(TARGETS _iriplas LIBRARY DESTINATION RMextract/pyiriplas)
