cmake_minimum_required (VERSION 3.1)
project (so3g)

include(local.cmake OPTIONAL)

# cmake policies -- best to keep these in sync with spt3g!
if(POLICY CMP0060) # Suppress cmake stripping full paths from libraries in some cases
  cmake_policy(SET CMP0060 NEW)
endif()
cmake_policy(SET CMP0012 NEW) # Allow use of true in boolean expressions
if(POLICY CMP0042) # Enable RPATH on OSX
   cmake_policy(SET CMP0042 NEW)
endif()

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
find_package(FLAC)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
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
  src/G3SuperTimestream.cxx
  src/so_linterp.cxx
  src/exceptions.cxx
  src/array_ops.cxx
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
target_link_libraries(so3g spt3g::core FLAC)

# You probably want to select openblas, so pass -DBLA_VENDOR=OpenBLAS
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message("-- BLAS found: ${BLAS_LIBRARIES}")
  target_link_libraries(so3g ${BLAS_LIBRARIES})

  # The BLAS library may or may not include the cblas_ bindings.
  # This variable set is needed by check_function_exists; starting in
  # cmake v3.18 you can say BLAS::BLAS instead of the lib path...
  set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES})
  check_function_exists(cblas_sgemm CBLAS_OK)
  if(${CBLAS_OK})
    message("-- cblas bindings are included in the BLAS library")
  else()
    message("-- cblas bindings not found in BLAS; adding cblas.")
    target_link_libraries(so3g cblas)
  endif()

  # On MacOS with clang linking to the Accelerate framework, the cblas
  # headers are not always found.  Handle this case.
  if(BLAS_Accelerate_LIBRARY)
    #string(REGEX REPLACE "^(.*)/System/Library/Frameworks/Accelerate.framework" "\\1" FRAMEWORK_INCLUDE "${BLAS_Accelerate_LIBRARY}")
    #target_compile_options(so3g PRIVATE -F ${FRAMEWORK_INCLUDE})
    target_include_directories(so3g PRIVATE ${BLAS_Accelerate_LIBRARY}/Versions/A/Frameworks/vecLib.framework/Headers)
  endif()
endif(BLAS_FOUND)

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
