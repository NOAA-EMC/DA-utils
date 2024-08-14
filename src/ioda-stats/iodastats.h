#pragma once

#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Group.h"
#include "ioda/ObsDataIoParameters.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/PostProcessor.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/TimeWindow.h"

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
      int execute(const eckit::Configuration & fullConfig, bool /*validate*/) const {
        // define the time window
        const eckit::LocalConfiguration timeWindowConf(fullConfig, "time window");
        const util::TimeWindow timeWindow(timeWindowConf);

        // get the list of obs spaces to process
        std::vector<eckit::LocalConfiguration> obsSpaces;
        fullConfig.get("obs spaces", obsSpaces);

        // for now, just do this serially, eventually, make it parallelized
        for (int i = 0; i < obsSpaces.size(); i++) {
          // get the configuration for this obs space
          auto obsSpace = obsSpaces[i];
          eckit::LocalConfiguration obsConfig(obsSpace, "obs space");

          // open the IODA file
          std::string obsFile;
          obsConfig.get("obsdatain.engine.obsfile", obsFile);
          oops::Log::info() << "IODA-Stats: Processing " << obsFile << std::endl;
          ioda::ObsSpace ospace(obsConfig, getComm(), timeWindow, getComm());
          const size_t nlocs = ospace.nlocs();
          oops::Log::info() << obsFile << ": nlocs =" << nlocs << std::endl;

          // get the list of variables (and channels if applicable) to process
          std::vector<std::string> variables;
          
        }
      }
  };
}  // namespace dautils
