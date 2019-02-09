# Find the core and cereal modules from spt3g.  This relies on
# environment variables
#
#   SPT3G_SOFTWARE_PATH -- path to the source tree
#   SPT3G_SOFTWARE_BUILD_PATH -- path to the build directory
#
# One way to set these variables is with env-shell.sh in spt3g build
# tree.

if (SPT3G_INCLUDES)
  set (SPT3G_FIND_QUIETLY TRUE)
endif (SPT3G_INCLUDES)

find_path(SPT3G_CORE_INCLUDE G3Pipeline.h PATHS
  $ENV{SPT3G_SOFTWARE_PATH}/core/include/core/)
find_path(SPT3G_CEREAL_INCLUDE cereal/cereal.hpp PATHS
  $ENV{SPT3G_SOFTWARE_PATH}/core/include)

find_library(SPT3G_CORE_LIBRARY core.so PATHS
  $ENV{SPT3G_SOFTWARE_BUILD_PATH}/spt3g/)

if(NOT SPT3G_CORE_LIBRARY)
  find_library(SPT3G_CORE_LIBRARY core.dylib PATHS $ENV{SPT3G_SOFTWARE_BUILD_PATH}/spt3g/)
endif(NOT SPT3G_CORE_LIBRARY)


include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (SPT3G DEFAULT_MSG 
  SPT3G_CORE_INCLUDE SPT3G_CORE_LIBRARY )

mark_as_advanced(SPT3G_INCLUDES SPT3G_LIBRARIES)

set(SPT3G_INCLUDES ${SPT3G_CORE_INCLUDE} ${SPT3G_CEREAL_INCLUDE})
set(SPT3G_LIBRARIES ${SPT3G_CORE_LIBRARY})

message(STATUS ${SPT3G_INCLUDES} ${SPT3G_LIBRARIES})
