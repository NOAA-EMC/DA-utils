#pragma once

namespace dautils {
  class CalcIodaStats {
    public :
      // -----------------------------------------------------------------------------
      CalcIodaStats(const eckit::Configuration & fullConfig,
                    const eckit::mpi::Comm & comm = oops::mpi::world())
        {}
    public:
      void run();
  };

} // namespace dautils
