if(CMAKE_PROJECT_NAME)
    set(target_name "${CMAKE_PROJECT_NAME}_ioda_stats.x")
else()
    set(target_name "ioda_stats.x")
endif()

ecbuild_add_executable( TARGET ${target_name}
                        SOURCES ioda_stats.cc )

target_compile_features(${target_name} PUBLIC cxx_std_17)
target_link_libraries(${target_name} PUBLIC NetCDF::NetCDF_CXX oops ioda)
