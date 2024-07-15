ecbuild_add_executable( TARGET meanioda.x
                        SOURCES meanioda.cc )

target_compile_features( meanioda.x PUBLIC cxx_std_17)
target_link_libraries( meanioda.x PUBLIC oops ioda)

