# - Find the Eigen library
# 
# EIGEN_FOUND
# EIGEN_INCLUDE_DIRS
# EIGEN_LIBRARIES
# EIGEN_LIBRARY_DIRS

if (EIGEN_INCLUDE_DIR)
# Already in cache, be silent
set (EIGEN_FIND_QUIETLY TRUE)
endif (EIGEN_INCLUDE_DIR)

find_path(EIGEN_INCLUDE_DIR "Eigen/Core"
    HINTS ENV EIGEN_DIR
    PATH_SUFFIXES eigen3
)

add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES 
    INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_INCLUDE_DIR}")

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Eigen3 DEFAULT_MSG EIGEN_INCLUDE_DIR)

mark_as_advanced(EIGEN_INCLUDE_DIR)
