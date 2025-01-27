# these directories require OOPS to build
if ( oops_FOUND )
  message( STATUS "Building tools that link to OOPS/JEDI" )
  add_subdirectory( ioda-extract )
  add_subdirectory( ioda-stats )
endif()
if ( BUILD_REGRID_STATES )
  add_subdirectory( regridStates )
endif()
