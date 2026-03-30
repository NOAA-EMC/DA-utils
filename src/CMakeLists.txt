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
# copy unified restriction script into build/bin/
# ----------------------------------------------------------------------
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/bin/ioda_restriction_filter.py
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${CMAKE_CURRENT_SOURCE_DIR}/ioda-restrict/ioda_restriction_filter.py
            ${CMAKE_BINARY_DIR}/bin/ioda_restriction_filter.py
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/ioda-restrict/ioda_restriction_filter.py
)

add_custom_target(copy_restriction ALL
    DEPENDS ${CMAKE_BINARY_DIR}/bin/ioda_restriction_filter.py
)

