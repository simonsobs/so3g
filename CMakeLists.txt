cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX C)

include(local.cmake OPTIONAL)

# Warn if the user invokes CMake directly
if (NOT SKBUILD)
  message(WARNING "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

# # Default to Release because we want that -O3.
# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE "Release" CACHE STRING
#     "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel" FORCE)
# endif(NOT CMAKE_BUILD_TYPE)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# We are building libraries that will eventually be linked into shared
# modules.  All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(
  Python 3.9
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
)

find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(libso3g STABLE_ABI
  src/main.cxx
  src/hkagg.cxx
  src/exceptions.cxx
  # src/Intervals.cxx
  # src/Butterworth.cxx
  # src/Ranges.cxx
  # src/Projection.cxx
  # src/so_linterp.cxx
  # src/array_ops.cxx
)

# Install directive for scikit-build-core
install(TARGETS libso3g LIBRARY DESTINATION so3g)

# Numpy include directory
execute_process(COMMAND ${Python_EXECUTABLE} -c
  "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR} )
include_directories(${NUMPY_INCLUDE_DIR})

# External Dependencies

# OpenMP

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else()
  message(WARNING "OpenMP not being linked -- this may affect performance.")
endif()

# GSL

find_package(GSL)
target_include_directories(libso3g PRIVATE ${GSL_INCLUDE_DIR})
target_link_libraries(libso3g PUBLIC ${GSL_LIBRARIES})

# You probably want to select openblas, so pass -DBLA_VENDOR=OpenBLAS
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message("-- BLAS found: ${BLAS_LIBRARIES}")
  target_link_libraries(libso3g PUBLIC ${BLAS_LIBRARIES})

  # The BLAS library may or may not include the cblas_* bindings.
  # This variable set is needed by check_function_exists; starting in
  # cmake v3.18 you can say BLAS::BLAS instead of the lib path...
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  check_function_exists(cblas_sgemm CBLAS_OK)
  if(${CBLAS_OK})
    message("-- cblas bindings are included in the BLAS library")
  else()
    message("-- cblas bindings not found in BLAS; adding cblas.")
    target_link_libraries(libso3g PUBLIC cblas)
  endif()

  # On MacOS with clang linking to the Accelerate framework, the cblas
  # headers are not always found.  Handle this case.  Also note that the
  # Accelerate framework has documented numerical problems- consider using
  # a better BLAS/LAPACK implementation.
  if(BLAS_Accelerate_LIBRARY)
    target_include_directories(libso3g PRIVATE ${BLAS_Accelerate_LIBRARY}/Versions/A/Frameworks/vecLib.framework/Headers)
  endif()
endif(BLAS_FOUND)

# This custom target generates _version.h, in the build tree.  That is all.
add_custom_target(so3g-version
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/version_h.py
  SO3G_VERSION_STRING ${CMAKE_CURRENT_BINARY_DIR}/_version.h
  SOURCES version_h.py
)
add_dependencies(libso3g so3g-version)

# # Make a list of .py files for the library.
# file(GLOB MY_PYTHONS
#   "${CMAKE_CURRENT_SOURCE_DIR}/python/*.py")
# file(GLOB MY_PYTHONS_HK
#   "${CMAKE_CURRENT_SOURCE_DIR}/python/hk/*.py")
# file(GLOB MY_PYTHONS_PROJ
#   "${CMAKE_CURRENT_SOURCE_DIR}/python/proj/*.py")
# file(GLOB MY_PYTHONS_SMURF
#   "${CMAKE_CURRENT_SOURCE_DIR}/python/smurf/*.py")

# # Define the install rules.

# if(DEFINED PYTHON_INSTALL_DEST)
#   get_filename_component(INSTALL_DEST ${PYTHON_INSTALL_DEST}/so3g
#                          ABSOLUTE BASE_DIR ${CMAKE_BINARY_DIR})
#   message("local.cmake has specified the install dir: ${INSTALL_DEST}")
# else()
#   set(INSTALL_DEST ${PYTHON_SITE_PACKAGES}/so3g)
# endif()

# install(TARGETS so3g
#         DESTINATION ${INSTALL_DEST})

# install(FILES ${MY_PYTHONS}
#         DESTINATION ${INSTALL_DEST})
# install(FILES ${MY_PYTHONS_HK}
#         DESTINATION ${INSTALL_DEST}/hk)
# install(FILES ${MY_PYTHONS_PROJ}
#         DESTINATION ${INSTALL_DEST}/proj)
# install(FILES ${MY_PYTHONS_SMURF}
#         DESTINATION ${INSTALL_DEST}/smurf)

# # To add a prefix, pass CMAKE_INSTALL_PREFIX.
# install(PROGRAMS scripts/so-hk-tool DESTINATION bin)

add_custom_target(prep-readthedocs
  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/docs/extract_docstrings.py
          --prep-rtd --source-branch=master
  )
