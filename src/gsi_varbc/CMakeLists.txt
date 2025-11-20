
ecbuild_add_executable( TARGET  satbias2ioda.x
                        SOURCES SatBiasConverter.cpp GsiSatBiasReader.cpp GsiSatBiasReader.h
                        LIBS    Eigen3::Eigen ${oops_LIBRARIES} ioda_engines )

ecbuild_add_executable( TARGET  acftbias2ioda.x
                        SOURCES GsiAircraftBiasConverter.cc GsiAircraftBiasReader.cc GsiAircraftBiasReader.h
                        LIBS    Eigen3::Eigen ${oops_LIBRARIES} ioda_engines )

