#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"

#include "ioda_stats_driver.h"

namespace dautils {
  class IodaStats : public oops::Application {
    public:
      // -----------------------------------------------------------------------------
      explicit IodaStats(const eckit::mpi::Comm & comm = oops::mpi::world())
        : Application(comm), fillVal_(util::missingValue<float>()) {
        oceans_["Atlantic"] = 1;
        oceans_["Pacific"] = 2;
        oceans_["Indian"] = 3;
        oceans_["Arctic"] = 4;
        oceans_["Southern"] = 5;
        }
      static const std::string classname() {return "dautils::IodaStats";}
      // -----------------------------------------------------------------------------
      int execute(const eckit::Configuration & fullConfig) const {
      // Initialize the application
      // Pass the full configuration to the driver class
      dautils::IodaStatsDriver driver(fullConfig, this->getComm());
      driver.run();
      return 0;
    }

   private:
      // -----------------------------------------------------------------------------
      std::string appname() const {
        return "dautils::IodaStats";
    }
  };
}  // namespace dautils