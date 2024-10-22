# these directories require OOPS to build
if ( oops_FOUND )
  message( STATUS "Building tools that link to OOPS/JEDI" )
  add_subdirectory( ioda-extract )
  add_subdirectory( ioda-stats )
endif()
add_subdirectory( regridStates )
add_subdirectory( scf-proc )
add_subdirectory( land-add-incr )
