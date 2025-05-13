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

find_package(OpenMP)
find_package(FLAC)
find_package(GSL)
find_package(Ceres)
find_package(nanobind CONFIG REQUIRED)

set(libso3g_SOURCES
    src/main.cxx
    src/hkagg.cxx
    src/Intervals.cxx
    src/Butterworth.cxx
    src/Ranges.cxx
    src/Projection.cxx
    src/so_linterp.cxx
    src/exceptions.cxx
    src/array_ops.cxx
    src/fitting_ops.cxx
)

# Documentation for the options here:
# https://nanobind.readthedocs.io/en/latest/api_cmake.html
nanobind_add_module(libso3g NOMINSIZE NB_STATIC ${libso3g_SOURCES})

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

# Add external dependencies to our include / link commands

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else()
    message(WARNING "OpenMP not being linked -- this may affect performance.")
endif()

target_link_libraries(libso3g PUBLIC Ceres::ceres Eigen3::Eigen)

target_include_directories(libso3g PRIVATE ${FLAC_INCLUDE_DIR})

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
        target_include_directories(libso3g PRIVATE
            ${BLAS_Accelerate_LIBRARY}/Versions/A/Frameworks/vecLib.framework/Headers
        )
    endif()
endif(BLAS_FOUND)

# This custom target generates _version.h, in the build tree.  That is all.
add_custom_target(so3g-version
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/version_h.py
    SO3G_VERSION_STRING ${CMAKE_CURRENT_BINARY_DIR}/_version.h
    SOURCES version_h.py
)
add_dependencies(libso3g so3g-version)

add_custom_target(prep-readthedocs
    COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/docs/extract_docstrings.py
            --prep-rtd --source-branch=master
)
