# This project name is cycstub
PROJECT(bright_cyclus)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Direct any out-of-source builds to this directory
SET( CYCAMORE_SOURCE_DIR ${CMAKE_SOURCE_DIR} )

# Direct any binary installation paths to this directory
SET( CYCAMORE_BINARY_DIR ${CMAKE_BINARY_DIR} )

# This makes all the libraries build as SHARED
SET(BUILD_SHARED_LIBS true)

# Setup build locations.
IF(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CYCAMORE_BINARY_DIR}/bin)
endif()
IF(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CYCAMORE_BINARY_DIR}/lib)
endif()
IF(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CYCAMORE_BINARY_DIR}/lib)
ENDIF()

SET(CYCAMORE_EXECUTABLE_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/cycstub/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/cycstub/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/cycstub/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Tell CMake where the modules are
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_DIR}/share/cmake-2.8/Modules")
SET( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CYCAMORE_SOURCE_DIR}/CMake )

# Include macros
INCLUDE(ModuleMacros)

# Debian installs useful LibXML2 files to /usr/include/libxml2/libxml
# KDHFLAG Necessary?
FIND_PACKAGE( LibXml2 REQUIRED )
SET(CYCAMORE_INCLUDE_DIR ${CYCAMORE_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})
ADD_DEFINITIONS (${LIBXML2_DEFINITIONS})
SET(LIBS ${LIBS} ${LIBXML2_LIBRARIES})

# Find Sqlite3
# KDHFLAG Necessary?
FIND_PACKAGE( Sqlite3 REQUIRED )
SET(CYCAMORE_INCLUDE_DIR ${CYCAMORE_INCLUDE_DIR} ${SQLITE3_INCLUDE_DIR})
SET(LIBS ${LIBS} ${SQLITE3_LIBRARIES})

# Find cycluscore
FIND_PACKAGE( cycluscore REQUIRED )
SET(CYCAMORE_INCLUDE_DIR ${CYCAMORE_INCLUDE_DIR} ${CYCLUS_CORE_INCLUDE_DIR})
SET(LIBS ${LIBS} ${CYCLUS_CORE_LIBRARIES})

# Include the boost header files and the program_options library
# KDHFLAG Necessary?
SET(Boost_USE_STATIC_LIBS       OFF)
SET(Boost_USE_STATIC_RUNTIME    OFF)
FIND_PACKAGE( Boost COMPONENTS program_options filesystem system REQUIRED)
SET(CYCAMORE_INCLUDE_DIR ${CYCAMORE_INCLUDE_DIR} ${BOOST_INCLUDE_DIR})
SET(LIBS ${LIBS} ${Boost_PROGRAM_OPTIONS_LIBRARY})
SET(LIBS ${LIBS} ${Boost_SYSTEM_LIBRARY})
SET(LIBS ${LIBS} ${Boost_FILESYSTEM_LIBRARY})

# include the model directories
SET(CYCAMORE_INCLUDE_DIR ${CYCAMORE_INCLUDE_DIR} Testing )

# Find the Pyne installation
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/../cmake)
find_package(Pyne REQUIRED)
link_directories(${PYNE_LIBS_DIR})
include_directories(${PYNE_INCLUDE_DIR})


# include all the directories we just found
INCLUDE_DIRECTORIES( ${CYCAMORE_INCLUDE_DIR} )
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../cpp")


# ------------------------- Add the Models -----------------------------------
ADD_SUBDIRECTORY(Models)
