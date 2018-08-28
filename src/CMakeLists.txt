set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
#find_package(GLFW REQUIRED)
#find_package(OpenGL REQUIRED)
#find_package(GLEW REQUIRED)
#find_package(FFTW REQUIRED)
#find_package(X11 REQUIRED)

find_package(Spt3g REQUIRED)

find_package(PythonInterp)
find_package(PythonLibs)

set(Boost_PYTHON_TYPE python)
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
find_package(Boost COMPONENTS system iostreams filesystem ${Boost_PYTHON_TYPE} REQUIRED)

include_directories(${SPT3G_INCLUDES})
include_directories(${Boost_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} ${PYTHON_INCLUDE_PATH})
include_directories(${CMAKE_SOURCE_DIR}/include)

# Remove the lib* prefix, since this is a python lib??
##set_target_properties(pylib PROPERTIES SUFFIX .so)

#set(CMAKE_SHARED_LIBRARY_PREFIX "")
add_library(so3g SHARED
#  test.cxx
  G3IndexedReader.cxx
  SampleFlags.cxx
)

set_target_properties(so3g PROPERTIES PREFIX "")

macro(link_python_dir)
        execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
	if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR}/python ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PROJECT})
	else(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
		execute_process(COMMAND ln -fsn ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${PROJECT})
	endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/python)
endmacro(link_python_dir)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spt3g)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

target_link_libraries(so3g pthread ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${SPT3G_LIBRARIES})
link_python_dir()
