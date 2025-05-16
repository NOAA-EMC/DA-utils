ecbuild_add_executable( TARGET ${CMAKE_PROJECT_NAME}_ioda_diag-stats.x
                        SOURCES ioda_diag-stats.cc )

target_compile_features( ${CMAKE_PROJECT_NAME}_ioda_diag-stats.x PUBLIC cxx_std_17)
target_link_libraries( ${CMAKE_PROJECT_NAME}_ioda_diag-stats.x PUBLIC NetCDF::NetCDF_CXX oops ioda)