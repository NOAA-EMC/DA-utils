ecbuild_add_executable( TARGET regridStates.x
                        SOURCES regridStates.F90
                                fv3_grid.F90
                                utils.F90)

target_link_libraries(
  regridStates.x
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

