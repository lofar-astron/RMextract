cmake_minimum_required(VERSION 3.18)
project(RMextract C CXX Fortran)

find_package(PythonExtensions REQUIRED)
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(F2PY REQUIRED)

## -- Python module: RMextract.EMM._EMM_Model -- ##
add_library(_EMM_Model MODULE
    RMextract/EMM/EMM_Model.cc
    RMextract/EMM/EMM_Model_wrap.cc
    RMextract/EMM/GeomagnetismLibrary.c)
python_extension_module(_EMM_Model)
install(TARGETS _EMM_Model LIBRARY DESTINATION RMextract/EMM)

## -- Python module: RMextract.pyiri._iri -- ##
set(F2PY_MODULE_NAME "_iri")
set(F2PY_MODULE_FILE ${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(F2PY_MODULE_SRC
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/iri.pyf
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/cira.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/igrf.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/iridreg.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/iriflip.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/irifun.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/irisub.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/iritec.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiri/iriget.for)
add_custom_target(${F2PY_MODULE_NAME} ALL
    DEPENDS ${F2PY_MODULE_FILE}
)
add_custom_command(
    OUTPUT ${F2PY_MODULE_FILE}
    COMMAND ${F2PY_EXECUTABLE} -m ${F2PY_MODULE_NAME} -c ${F2PY_MODULE_SRC}
    COMMENT "[F2PY] Building python extension module ${F2PY_MODULE_NAME}"
)
install(FILES ${F2PY_MODULE_FILE} DESTINATION RMextract/pyiri)

## -- Python module: RMextract.pyiriplas._iriplas -- ##
set(F2PY_MODULE_NAME "_iriplas")
set(F2PY_MODULE_FILE ${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_NAME}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(F2PY_MODULE_SRC
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/Iris2017.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/igrf.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/indx2017.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/irif2019.for
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/iriplas.pyf
    ${CMAKE_SOURCE_DIR}/RMextract/pyiriplas/iriplas_main.for)
add_custom_target(${F2PY_MODULE_NAME} ALL
    DEPENDS ${F2PY_MODULE_FILE}
)
add_custom_command(
    OUTPUT ${F2PY_MODULE_FILE}
    COMMAND ${F2PY_EXECUTABLE} -m ${F2PY_MODULE_NAME} -c ${F2PY_MODULE_SRC}
    COMMENT "[F2PY] Building python extension module ${F2PY_MODULE_NAME}"
)
install(FILES ${F2PY_MODULE_FILE} DESTINATION RMextract/pyiriplas)
