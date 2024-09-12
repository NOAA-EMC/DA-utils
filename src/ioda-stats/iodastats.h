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

#include "./calcstats.h"
#include "./statfile.h"

namespace dautils {
  class IodaStats : public oops::Application {
    public:
      // -----------------------------------------------------------------------------
      explicit IodaStats(const eckit::mpi::Comm & comm = oops::mpi::world())
        : Application(comm), fillVal_(util::missingValue<float>()) {
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

          // get the lists of everything to process/compute
          std::vector<std::string> groups;
          std::vector<std::string> stats;
          std::vector<std::string> qcgroups;
          std::vector<eckit::LocalConfiguration> domains;
          
          obsSpace.get("groups to process", groups);
          obsSpace.get("qc groups", qcgroups);
          obsSpace.get("statistics to compute", stats);

          if (obsSpace.has("domains to proceses")) {
            obsSpace.get("domains to process", domains);
          }

          // optional regular binning configuration
          std::vector<std::vector<float>> lats1;
          std::vector<std::vector<float>> lons1;
          std::vector<std::vector<float>> lats2;
          std::vector<std::vector<float>> lons2;
          std::vector<std::vector<float>> lats_cen;
          std::vector<std::vector<float>> lons_cen;
          float bin_res;
          bool do_binning = false;
          if (obsSpace.has("binning resolution in degrees")) {
            obsSpace.get("binning resolution in degrees", bin_res);
            do_binning = true;
          }

          if (do_binning) {
            // create the bins, note this will assume global model coverage
            int nlons;
            int nlats;
            nlons = static_cast<int>(360.0 / bin_res);
            nlats = static_cast<int>(180.0 / bin_res);
          }
          
          // loop over all domains and get their definitions
          std::vector<std::string> domainNames;
          std::vector<std::string> domainMaskVar1;
          std::vector<std::string> domainMaskVar2;
          std::vector<std::string> domainMaskVar3;
          std::vector<std::vector<float>> domainMaskVals1;
          std::vector<std::vector<float>> domainMaskVals2;
          std::vector<std::vector<float>> domainMaskVals3;
          for (int idom = 0; idom < domains.size(); idom++ ) {
            auto domain = domains[idom];
            eckit::LocalConfiguration domainConf(domain, "domain");
            std::string domainname;
            std::string maskvar1, maskvar2, maskvar3;
            std::vector<float> maskvals1, maskvals2, maskvals3;
            domainConf.get("name", domainname);
            if (domainConf.has("first mask variable")) {
              domainConf.get("first mask variable", maskvar1);
              domainConf.get("first mask range", maskvals1);
            }
            if (domainConf.has("second mask variable")) {
              domainConf.get("second mask variable", maskvar2);
              domainConf.get("second mask range", maskvals2);
            }
            if (domainConf.has("third mask variable")) {
              domainConf.get("third mask variable", maskvar3);
              domainConf.get("third mask range", maskvals3);
            }
            domainNames.push_back(domainname);
            domainMaskVar1.push_back(maskvar1);
            domainMaskVar2.push_back(maskvar2);
            domainMaskVar3.push_back(maskvar3);
            domainMaskVals1.push_back(maskvals1);
            domainMaskVals2.push_back(maskvals2);
            domainMaskVals3.push_back(maskvals3);
          }

          // assert that the QC groups list is the same size as groups
          assert(groups.size() == qcgroups.size());

          // initialize netCDF output file for writing
          std::string outfile;
          obsSpace.get("output file", outfile);
          StatFile statfile;
          statfile.initializeNcfile(outfile, timeWindow, variables, channels, groups, stats, domainNames);

          // loop over domains, compute the masks for each
          std::vector<std::vector<int>> mask(domains.size()+1, std::vector<int>(nlocs, 0));
          for (int idom = 0; idom < domains.size(); idom++ ) {
            // compute mask with function 3 times, one for each possible mask
            ObsStats obstatmask;
            std::vector<float> maskvalues(nlocs);
            if (!domainMaskVar1[idom].empty()) {
              ospace.get_db("MetaData", domainMaskVar1[idom], maskvalues);
              mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals1[idom][0], domainMaskVals1[idom][1], mask[idom]);
            }
            if (!domainMaskVar2[idom].empty()) {
              ospace.get_db("MetaData", domainMaskVar2[idom], maskvalues);
              mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals2[idom][0], domainMaskVals2[idom][1], mask[idom]);
            }
            if (!domainMaskVar3[idom].empty()) {
              ospace.get_db("MetaData", domainMaskVar3[idom], maskvalues);
              mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals3[idom][0], domainMaskVals3[idom][1], mask[idom]);
            }
          }

          // loop over variables
          for (int var = 0; var < variables.size(); var++) {
            // loop over groups
            for (int g = 0; g < groups.size(); g++) {
              oops::Log::info() << obsFile << ": Now processing "
                                << groups[g] << "/" << variables[var] << std::endl;
              std::vector<float> buffer(nlocs);
              std::vector<int> qcflag(nlocs);
              // we have to process differently if there are channels
              if (channels.empty()) {
                // read the full variable
                ospace.get_db(groups[g], variables[var], buffer);
                // get the QC group
                ospace.get_db(qcgroups[g], variables[var], qcflag);
              } else {
                // give the list of channels to read
                ospace.get_db(groups[g], variables[var], buffer, channels);
                // get the QC group
                ospace.get_db(qcgroups[g], variables[var], qcflag, channels);

              }
              // loop over domains
              for (int idom = 0; idom < domains.size()+1; idom++ ) {
                if (idom < domains.size()) {
                  oops::Log::info() << "Processing domain: " << domainNames[idom] << std::endl;
                  oops::Log::info() << domainMaskVar1[idom] << "=" <<  domainMaskVals1[idom] << ";"
                           << domainMaskVar2[idom] << "=" <<  domainMaskVals2[idom] << ";"
                           << domainMaskVar3[idom] << "=" <<  domainMaskVals3[idom] << ";"
                           << std::endl;
                }
                // loop over stats
                ObsStats obstat;
                for (int s = 0; s < stats.size(); s++) {
                  // Maybe eventually set this up as a factory but for now just do it
                  // with this old school if/else if way
                  std::vector<int> intstat;
                  std::vector<float> floatstat;
                  if (stats[s] == "count") {
                    intstat = obstat.getObsCount(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "Count:" << intstat << std::endl; 
                  } else if (stats[s] == "mean") {
                    floatstat = obstat.getMean(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "Mean:" << floatstat << std::endl; 
                  } else if (stats[s] == "RMS") {
                    floatstat = obstat.getRMS(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "RMS:" << floatstat << std::endl;
                  } else {
                    oops::Log::info() << stats[s] << " not supported. Skipping." << std::endl;
                  }
                  if (stats[s] == "count") {
                    statfile.write(outfile, groups[g], variables[var],
                                   stats[s], idom, intstat);
                  } else {
                    statfile.write(outfile, groups[g], variables[var],
                                   stats[s], idom, floatstat);
                  }
                }
              }
            }
          }
        }
        return 0;
      }

    // -----------------------------------------------------------------------------
    // Data members
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
