cmake_minimum_required (VERSION 3.1)
project (so3g)

include(local.cmake OPTIONAL)

# Default to Release because we want that -O3.  This is what spt3g_software does too.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# For this to be found, make sure the spt3g build directory can be
# searched; i.e. -DCMAKE_PREFIX_PATH=/path/to/spt3g_software/build
find_package(Spt3g REQUIRED)

find_package(PythonInterp 3)
find_package(PythonLibs 3)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
  message(WARNING "OpenMP not being linked -- this may affect performance.")
endif()

# Determine the location of site-packages.
execute_process ( COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

# Numpy include directory?
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
                 OUTPUT_VARIABLE NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR} )
include_directories(${NUMPY_INCLUDE_DIR})

#
# Define the so3g build target.  This is a shared library.
#

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY so3g)
add_library(so3g SHARED
  src/main.cxx
  src/test.cxx
  src/hkagg.cxx
  src/G3Ndarray.cxx
  src/G3WCS.cxx
  src/G3Ndmap.cxx
  src/G3IndexedReader.cxx
  src/Intervals.cxx
  src/Butterworth.cxx
  src/Ranges.cxx
  src/Rebundler.cxx
  src/Projection.cxx
  src/exceptions.cxx
)

# We could disable the lib prefix on the output library... but let's not.
#set_target_properties(so3g PROPERTIES PREFIX "")

# Make a list of .py files for the library.
file(GLOB MY_PYTHONS
  "${CMAKE_CURRENT_SOURCE_DIR}/python/*.py")
file(GLOB MY_PYTHONS_HK
  "${CMAKE_CURRENT_SOURCE_DIR}/python/hk/*.py")
file(GLOB MY_PYTHONS_PROJ
  "${CMAKE_CURRENT_SOURCE_DIR}/python/proj/*.py")
file(GLOB MY_PYTHONS_SMURF
  "${CMAKE_CURRENT_SOURCE_DIR}/python/smurf/*.py")

# Provide list of libs to link against.
target_link_libraries(so3g spt3g::core)

# This custom target generates _version.h, in the build tree.  That is all.
add_custom_target(so3g-version
                  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/version_h.py
                    SO3G_VERSION_STRING ${CMAKE_CURRENT_BINARY_DIR}/_version.h
                  SOURCES version_h.py
                   )

add_dependencies(so3g so3g-version)

# Define the install rules.

if(DEFINED PYTHON_INSTALL_DEST)
  get_filename_component(INSTALL_DEST ${PYTHON_INSTALL_DEST}/so3g
                         ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR})
  message("local.cmake has specified the install dir: ${INSTALL_DEST}")
else()
  set(INSTALL_DEST ${PYTHON_SITE_PACKAGES}/so3g)
endif()

install(TARGETS so3g
        DESTINATION ${INSTALL_DEST})

install(FILES ${MY_PYTHONS}
        DESTINATION ${INSTALL_DEST})
install(FILES ${MY_PYTHONS_HK}
        DESTINATION ${INSTALL_DEST}/hk)
install(FILES ${MY_PYTHONS_PROJ}
        DESTINATION ${INSTALL_DEST}/proj)
install(FILES ${MY_PYTHONS_SMURF}
        DESTINATION ${INSTALL_DEST}/smurf)

add_custom_target(prep-readthedocs
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/docs/extract_docstrings.py
          --prep-rtd --source-branch=master
  )
