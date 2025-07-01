#pragma once

namespace dautils {
  class IodaStatsDriver {
    public:
     IodaStatsDriver(const eckit::Configuration & fullConfig, const eckit::mpi::Comm & comm)
       : fullConfig_(fullConfig), comm_(comm) {}
     void run();
  }
}
