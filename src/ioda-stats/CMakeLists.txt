ecbuild_add_executable( TARGET ioda-stats.x
                        SOURCES iodastats.cc )

target_compile_features( ioda-stats.x PUBLIC cxx_std_17)
target_link_libraries( ioda-stats.x PUBLIC NetCDF::NetCDF_CXX oops ioda)
