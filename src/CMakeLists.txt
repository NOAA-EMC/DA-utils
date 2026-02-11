# these directories require OOPS to build
if ( oops_FOUND )
  message( STATUS "Building tools that link to OOPS/JEDI" )
  add_subdirectory( ioda-extract )
  add_subdirectory( ioda-stats )
  add_subdirectory( ioda-summary )
  add_subdirectory( ioda-dump )
  add_subdirectory( gsi_varbc )
  if ( fv3jedi_FOUND )
    add_subdirectory( add_increment )
  endif()
endif()

# ----------------------------------------------------------------------
# copy check_rsrd.py and check_exprsrd.py into build/bin/
# ----------------------------------------------------------------------
add_custom_target(copy_restriction ALL
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/ioda-restrict/check_rsrd.py
                ${CMAKE_BINARY_DIR}/bin/check_rsrd.py
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/ioda-restrict/check_exprsrd.py
                ${CMAKE_BINARY_DIR}/bin/check_exprsrd.py
)
