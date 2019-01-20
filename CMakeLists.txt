cmake_minimum_required (VERSION 2.8)
project (so3g)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")

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

set(Boost_PYTHON_TYPE python)
set(Boost_NUMPY_TYPE numpy3)
if(${PYTHON_LIBRARIES} MATCHES "libpython3.*\\.so")
	# Handle possible naming of the Boost python library as
	# libboost_python3.so
	set(_Boost_PYTHON3_HEADERS "boost/python.hpp")
	find_package(Boost COMPONENTS python3)
	if (${Boost_PYTHON3_FOUND})
		set(Boost_PYTHON_TYPE python3)
	endif()
else()
        message("FYI: ${PYTHON_LIBRARIES}")
endif()
find_package(Boost COMPONENTS system iostreams filesystem ${Boost_NUMPY_TYPE}
                   ${Boost_PYTHON_TYPE} REQUIRED)

include_directories(${SPT3G_INCLUDES})
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR} )


#
# Define the so3g build target.  This is a shared library.
#

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY so3g)
add_library(so3g SHARED
  src/main.cxx
  src/test.cxx
  src/hkagg.cxx
  src/G3IndexedReader.cxx
  src/Intervals.cxx
  src/Butterworth.cxx
)

# We could disable the lib prefix on the output library... but let's not.
#set_target_properties(so3g PROPERTIES PREFIX "")

# Make a list of .py files for the library.
file(GLOB MY_PYTHONS "${CMAKE_CURRENT_SOURCE_DIR}/python/*.py")

# Provide list of libs to link against.
target_link_libraries(so3g pthread ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${SPT3G_LIBRARIES})

# This custom target generates _version.h, in the build tree.  That is all.
add_custom_target(so3g-version
                  COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/version_h.py
                    SO3G_VERSION_STRING ${CMAKE_CURRENT_BINARY_DIR}/_version.h
                  SOURCES version_h.py
                   )

add_dependencies(so3g so3g-version)
                   

#
# Define the install rules.
#

include(local.cmake OPTIONAL)
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
