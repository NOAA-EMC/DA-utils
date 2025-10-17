#pragma once

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

namespace dautils {
  class CalcIodaStats {
    public :
      // -----------------------------------------------------------------------------
      CalcIodaStats(const eckit::Configuration & fullConfig,
                    const eckit::mpi::Comm & comm = oops::mpi::world())
        : config_(fullConfig), comm_(comm)
        {}
    public:
      void run();
    private:
      const eckit::Configuration & config_;
      const eckit::mpi::Comm & comm_;
  };
  void convertLongitudes(std::vector<float>& longitudes);
} // namespace dautils
