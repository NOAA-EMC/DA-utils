if(CMAKE_PROJECT_NAME)
    set(serial_target_name "${CMAKE_PROJECT_NAME}_fv3jedi_add_increment.x")
    set(ensemble_target_name "${CMAKE_PROJECT_NAME}_fv3jedi_ensemble_add_increment.x")
else()
    set(serial_target_name "fv3jedi_add_increment.x")
    set(ensemble_target_name "fv3jedi_ensemble_add_increment.x")
endif()

# AddIncrement
ecbuild_add_executable( TARGET ${serial_target_name}
                        SOURCES fv3jedi_add_increment.cc add_increment.h)
target_compile_features( ${serial_target_name} PUBLIC cxx_std_17)
target_link_libraries( ${serial_target_name} PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)

# Ensemble AddIncrement
ecbuild_add_executable( TARGET ${ensemble_target_name}
                        SOURCES fv3jedi_ensemble_add_increment.cc add_increment.h)
target_compile_features( ${ensemble_target_name} PUBLIC cxx_std_17)
target_link_libraries( ${ensemble_target_name} PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)
