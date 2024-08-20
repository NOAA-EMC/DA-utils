# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# Create Data directory for test input config and symlink all files
list( APPEND utils_test_input
  testinput/meanioda.yaml
)

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testrun)

CREATE_SYMLINK( ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${utils_test_input} )

# tests that require OOPS/JEDI
if ( oops_FOUND )
  # Test example IODA utility that computes the mean of a variable
  ecbuild_add_test( TARGET  test_dautils_ioda_example
                    COMMAND ${CMAKE_BINARY_DIR}/bin/meanioda.x
                    ARGS    "testinput/meanioda.yaml"
                    LIBS    da-utils)
endif()
