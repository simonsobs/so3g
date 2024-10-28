# Try to find the Ceres Solver library

# Look for the Ceres package
find_path(CERES_INCLUDE_DIR ceres/ceres.h
    PATH_SUFFIXES ceres
)

find_library(CERES_LIBRARY ceres)

# Check if Eigen is needed and available
find_package(Eigen3 REQUIRED)

# Check if Ceres is found
if (CERES_INCLUDE_DIR AND CERES_LIBRARY AND TARGET Eigen3::Eigen)
    # Ceres and Eigen were found
    set(CERES_FOUND TRUE)
    set(CERES_LIBRARIES ${CERES_LIBRARY})
    set(CERES_INCLUDE_DIRS ${CERES_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})

    # Optionally, find other dependencies (like gflags and glog, which Ceres uses)
    find_package(Glog REQUIRED)
    find_package(Gflags REQUIRED)

    set(CERES_DEPENDENCIES ${GLOG_LIBRARIES} ${GFLAGS_LIBRARIES} Eigen3::Eigen)
else()
    # Ceres or Eigen was not found
    set(CERES_FOUND FALSE)
endif()

# Set the results so they can be used by the project
mark_as_advanced(CERES_INCLUDE_DIR CERES_LIBRARY)

# Provide an interface for usage
if (CERES_FOUND)
    message(STATUS "Found Ceres: ${CERES_INCLUDE_DIR}")
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
else()
    message(WARNING "Could not find Ceres Solver or Eigen3")
endif()