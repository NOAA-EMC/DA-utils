#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

#include "calc_iodastats.h"

namespace dautils {
  class IodaStats : public oops::Application {
    public:
      // -----------------------------------------------------------------------------
      explicit IodaStats(const eckit::mpi::Comm & comm = oops::mpi::world())
        : Application(comm) {}
      static const std::string classname() {return "dautils::IodaStats";}
      // -----------------------------------------------------------------------------
      int execute(const eckit::Configuration & fullConfig) const {
        // Initialize the application
        // Pass the full configuration to the calculation class
        dautils::CalcIodaStats calc(fullConfig, this->getComm());
        calc.run();        
        return 0;
      }

    private:
      std::string appname() const {
        return "dautils::IodaStats";
      }
    // -----------------------------------------------------------------------------
  };
}  // namespace dautils
