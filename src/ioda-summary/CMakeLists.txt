ecbuild_add_executable( TARGET ioda-summary.x
                        SOURCES ioda-summary.cc )

target_compile_features( ioda-summary.x PUBLIC cxx_std_17)
target_link_libraries( ioda-summary.x PUBLIC oops ioda)
