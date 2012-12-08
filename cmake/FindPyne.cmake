# - Find PyNE libraries
# This module finds the libraries corresponding to the PyNE library, using the Python
# interpreter.
# This code sets the following variables:
#
#  PYNE_LIBS_FOUND            - have the PyNE libs been found
#  PYNE_PREFIX                - path to the PyNE installation
#  PYNE_LIBRARIES             - path to the PyNE libs dir
#  PYNE_INCLUDE_DIRS          - path to where PyNE header files are

# Use the Python interpreter to find the libs.
if(Pyne_FIND_REQUIRED)
    find_package(PythonInterp REQUIRED)
else()
    find_package(PythonInterp)
endif()

if(NOT PYTHONINTERP_FOUND)
    set(PYNE_LIBS_FOUND FALSE)
    return()
endif()

execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import os; from pyne import pyne_config as pc; print(os.path.split(pc.__file__)[0]); print(pc.lib); print(pc.includes)"
    RESULT_VARIABLE _PYNE_SUCCESS
    OUTPUT_VARIABLE _PYNE_VALUES
    ERROR_VARIABLE _PYNE_ERROR_VALUE
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _PYNE_SUCCESS MATCHES 0)
    if(Pyne_FIND_REQUIRED)
        message(FATAL_ERROR
            "Pyne config failure:\n${_PYTHON_ERROR_VALUE}")
    endif()
    set(PYNE_LIBS_FOUND FALSE)
    return()
endif()

# Convert the process output into a list
string(REGEX REPLACE ";" "\\\\;" _PYTHON_VALUES ${_PYTHON_VALUES})
string(REGEX REPLACE "\n" ";" _PYTHON_VALUES ${_PYTHON_VALUES})
list(GET _PYTHON_VALUES 0 _PYTHON_VERSION_LIST)
list(GET _PYTHON_VALUES 1 PYTHON_PREFIX)
list(GET _PYTHON_VALUES 2 PYTHON_INCLUDE_DIR)
list(GET _PYTHON_VALUES 3 PYTHON_SITE_PACKAGES)
list(GET _PYTHON_VALUES 4 PYTHON_MODULE_EXTENSION)
list(GET _PYTHON_VALUES 5 PYTHON_IS_DEBUG)
list(GET _PYTHON_VALUES 6 PYTHON_SIZEOF_VOID_P)

# Make sure all directory separators are '/'
string(REGEX REPLACE "\\\\" "/" PYTHON_PREFIX ${PYTHON_PREFIX})
string(REGEX REPLACE "\\\\" "/" PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR})
string(REGEX REPLACE "\\\\" "/" PYTHON_SITE_PACKAGES ${PYTHON_SITE_PACKAGES})

# TODO: All the nuances of CPython debug builds have not been dealt with/tested.
if(PYTHON_IS_DEBUG)
    set(PYTHON_MODULE_EXTENSION "_d${PYTHON_MODULE_EXTENSION}")
endif()

MARK_AS_ADVANCED(
  PYTHON_LIBRARY
  PYTHON_INCLUDE_DIR
)

# We use PYTHON_INCLUDE_DIR, PYTHON_LIBRARY and PYTHON_DEBUG_LIBRARY for the
# cache entries because they are meant to specify the location of a single
# library. We now set the variables listed by the documentation for this
# module.
SET(PYTHON_INCLUDE_DIRS "${PYTHON_INCLUDE_DIR}")
SET(PYTHON_LIBRARIES "${PYTHON_LIBRARY}")
SET(PYTHON_DEBUG_LIBRARIES "${PYTHON_DEBUG_LIBRARY}")


# Don't know how to get to this directory, just doing something simple :P
#INCLUDE(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(PythonLibs DEFAULT_MSG PYTHON_LIBRARIES PYTHON_INCLUDE_DIRS)
find_package_message(PYNE
    "Found PythonLibs: ${PYTHON_LIBRARY}"
    "${PYTHON_EXECUTABLE}${PYTHON_VERSION}")
