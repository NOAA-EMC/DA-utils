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

#include "./statfile.h"

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
          std::vector<int> channels;
          obsSpace.get("variables", variables);
          if (obsSpace.has("channels")) {
            obsSpace.get("channels", channels);
          }

          // channels only works if there is one variable, so need to check this
          if (variables.size() > 1 && !channels.empty()) {
            throw eckit::Exception("Cannot use channels with multiple variables.");
          }

          // get the list of groups and statistics to process/compute
          std::vector<std::string> groups;
          std::vector<std::string> stats;
          std::vector<std::string> qcgroups;
          obsSpace.get("groups to process", groups);
          obsSpace.get("qc groups", qcgroups);
          obsSpace.get("statistics to compute", stats);

          // assert that the QC groups list is the same size as groups
          assert(groups.size() == qcgroups.size());

          // initialize netCDF output file for writing
          std::string outfile;
          obsSpace.get("output file", outfile);
          StatFile statfile;
          statfile.initializeNcfile(outfile, timeWindow, variables, channels, groups, stats);
          //statfile.writeTest(outfile, timeWindow);

          // Loop over variables
          for (int var = 0; var < variables.size(); var++) {
            // loop over groups
            for (int g = 0; g < groups.size(); g++) {
              oops::Log::info() << obsFile << ": Now processing "
                                << groups[g] << "/" << variables[var] << std::endl;
              std::vector<float> buffer(nlocs);
              // we have to process differently if there are channels
              if (channels.empty()) {
                // read the full variable
                ospace.get_db(groups[g], variables[var], buffer);
              } else {
                // give the list of channels to read
                ospace.get_db(groups[g], variables[var], buffer, channels);
              }
              // To-Do, bin by region/basin/etc.
              // loop over stats
              for (int s = 0; s < stats.size(); s++) {
                oops::Log::info() << "Now computing " << stats[s] << std::endl;

              }
            }
          }
        }
      }
    // -----------------------------------------------------------------------------
    // Data members
    std::map<std::string, int> oceans_;
    double fillVal_;
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "dautils::IodaExample";
    }
    // -----------------------------------------------------------------------------
  };
}  // namespace dautils
