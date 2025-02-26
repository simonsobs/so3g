# - Find FFTW
# Find the native FFTW includes and library
#
# FFTW_INCLUDES - where to find fftw3.h
# FFTW_LIBRARIES - List of libraries when using FFTW.
# FFTW_FOUND - True if FFTW found.
if (FFTW_INCLUDES)
# Already in cache, be silent
set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)
find_path (FFTW_INCLUDES fftw3.h HINTS ENV FFTW_INC)
find_library (FFTW_LIBRARIES NAMES fftw3 HINTS ENV FFTW_DIR)

#handles finding omp implementations
find_library (FFTW_OMP_LIBRARY NAMES fftw3_omp HINTS ENV FFTW_DIR)

include (FindPackageHandleStandardArgs)
set(FPHSA_NAME_MISMATCHED TRUE)
find_package_handle_standard_args(FFTW_OMP DEFAULT_MSG FFTW_OMP_LIBRARY)
mark_as_advanced(FFTW_OMP_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)
mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)

# FFTWF_LIBRARIES - List of libraries when using FFTWF.
# FFTWF_FOUND - True if FFTWF is found.
# if (FFTW_INCLUDES)
  # Already in cache, be silent
#  set (FFTWF_FIND_QUIETLY TRUE)
# endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3f.h HINTS ENV FFTW_INC)
find_library (FFTWF_LIBRARIES NAMES fftw3f HINTS ENV FFTW_DIR)

# Handles finding omp implementations
find_library (FFTWF_OMP_LIBRARY NAMES fftw3f_omp HINTS ENV FFTW_DIR)

include (FindPackageHandleStandardArgs)
set(FPHSA_NAME_MISMATCHED TRUE)
find_package_handle_standard_args(FFTWF_OMP DEFAULT_MSG FFTWF_OMP_LIBRARY)
mark_as_advanced(FFTWF_OMP_LIBRARY)

# Handle the QUIETLY and REQUIRED arguments and set FFTWF_FOUND to TRUE if
# all listed variables are TRUE
find_package_handle_standard_args (FFTWF DEFAULT_MSG FFTWF_LIBRARIES FFTW_INCLUDES)
mark_as_advanced (FFTWF_LIBRARIES FFTW_INCLUDES)