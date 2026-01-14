# AddIncrement
ecbuild_add_executable( TARGET gdas_fv3jedi_add_increment.x
                        SOURCES gdas_fv3jedi_add_increment.cc gdas_add_increment.h)
target_compile_features( gdas_fv3jedi_add_increment.x PUBLIC cxx_std_17)
target_link_libraries( gdas_fv3jedi_add_increment.x PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)

# Ensemble AddIncrement
ecbuild_add_executable( TARGET gdas_fv3jedi_ensemble_add_increment.x
                        SOURCES gdas_fv3jedi_ensemble_add_increment.cc gdas_add_increment.h)
target_compile_features( gdas_fv3jedi_ensemble_add_increment.x PUBLIC cxx_std_17)
target_link_libraries( gdas_fv3jedi_ensemble_add_increment.x PUBLIC NetCDF::NetCDF_CXX oops fv3jedi)
