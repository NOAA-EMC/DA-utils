# AddIncrement
ecbuild_add_executable( TARGET fv3jedi_add_increment.x
                        SOURCES fv3jedi_add_increment.cc add_increment.h)
target_compile_features( fv3jedi_add_increment.x PUBLIC cxx_std_17)
target_link_libraries( fv3jedi_add_increment.x PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)

# Ensemble AddIncrement
ecbuild_add_executable( TARGET fv3jedi_ensemble_add_increment.x
                        SOURCES fv3jedi_ensemble_add_increment.cc add_increment.h)
target_compile_features( fv3jedi_ensemble_add_increment.x PUBLIC cxx_std_17)
target_link_libraries( fv3jedi_ensemble_add_increment.x PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)
