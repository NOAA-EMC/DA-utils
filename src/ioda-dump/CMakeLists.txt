ecbuild_add_executable( TARGET ioda-dump.x
                        SOURCES ioda-dump.cc )

target_compile_features( ioda-dump.x PUBLIC cxx_std_17)
target_link_libraries( ioda-dump.x PUBLIC oops ioda)
