ecbuild_add_executable( TARGET iodadump.x
                        SOURCES iodadump.cc )

target_compile_features( iodadump.x PUBLIC cxx_std_17)
target_link_libraries( iodadump.x PUBLIC oops ioda)