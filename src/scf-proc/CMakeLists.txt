ecbuild_add_executable( TARGET calcIMS.x
	SOURCES driver_fIMS.f90
	IMSaggregate_mod.f90)

target_link_libraries(
	calcIMS.x
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)
