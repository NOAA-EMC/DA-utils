ecbuild_add_executable( TARGET apply_incr.exe
	SOURCES apply_incr_noahmp_snow.f90
	NoahMPdisag_module.f90)

target_link_libraries(
	apply_incr.exe	
  PUBLIC
  esmf
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)
