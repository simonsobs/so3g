# - Find the Ceres Solver library
#
# CERES_FOUND
# CERES_INCLUDE_DIRS
# CERES_LIBRARIES
# CERES_LIBRARY_DIRS

# Look for the Ceres package
find_path(CERES_INCLUDE_DIR NAMES ceres/ceres.h HINTS ENV CERES_DIR PATH_SUFFIXES include)
find_library(CERES_LIBRARY NAMES ceres HINTS ENV CERES_DIR PATH_SUFFIXES lib)

# Get Dependencies
find_package(Eigen3 REQUIRED)
find_package(Glog REQUIRED)
find_package(Gflags REQUIRED)

# Create the imported Ceres target
add_library(Ceres::Ceres UNKNOWN IMPORTED)
set_target_properties(Ceres::Ceres PROPERTIES
    IMPORTED_LOCATION "${CERES_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CERES_INCLUDE_DIR};${EIGEN3_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${GLOG_LIBRARIES};${GFLAGS_LIBRARIES};${CERES_LIBRARY}"
)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Ceres DEFAULT_MSG CERES_LIBRARY CERES_INCLUDE_DIR)

# Set the results so they can be used by the project
mark_as_advanced(CERES_INCLUDE_DIR CERES_LIBRARY)