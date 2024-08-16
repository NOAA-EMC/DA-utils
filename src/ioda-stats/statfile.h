#pragma once

#include <netcdf>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"

namespace dautils {
  class StatFile {
    public:

    int initializeNcfile(const std::string filename, const util::TimeWindow timeWindow,
                      std::vector<std::string> variables, std::vector<int> channels,
                      std::vector<std::string> groups, std::vector<std::string> stats) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::replace);
      oops::Log::info() << "Opening " << filename << " for writing..." << std::endl;
      // create an unlimited time dimension
      netCDF::NcDim tDim = ncFile.addDim("analysisCycle");
      // if channel is not empty, create a channel dimension
      if (!channels.empty()) {
        ncFile.addDim("Channel", channels.size());
      }
      // create validTime variable
      netCDF::NcVar time = ncFile.addVar("validTime", netCDF::ncString, tDim);
      // put the analysis time in the file
      util::DateTime analysisTime = timeWindow.midpoint();
      std::vector<size_t> idxout;
      idxout.push_back(0);
      time.putVar(idxout, analysisTime.toString());

      // loop over group, then variables, then stats to create /group/var/stat in file
      for (int g = 0; g < groups.size(); g++) {
        // create group group
        netCDF::NcGroup group = ncFile.addGroup(groups[g]);
        // loop over variables
        for (int var = 0; var < variables.size(); var++) {
          // create variable group
          netCDF::NcGroup group2 = group.addGroup(variables[var]);
          // loop over statistics to write out
          for (int s = 0; s < stats.size(); s++) {
            netCDF::NcVar varout;
            if (stats[s] == "count") {
              varout = group2.addVar(stats[s], netCDF::ncInt, tDim);
            } else {
              varout = group2.addVar(stats[s], netCDF::ncFloat, tDim);
            }
          }
        }
      }

      return 0;
    };
    // int writeTest(const std::string filename, const util::TimeWindow timeWindow) {
    //   netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
    //   netCDF::NcVar time = ncFile.getVar("validTime");
    //   util::DateTime analysisTime = timeWindow.midpoint();
    //   std::vector<size_t> idxout;
    //   idxout.push_back(0);
    //   oops::Log::info() << "stuff " << analysisTime << std::endl;
    //   time.putVar(idxout, analysisTime.toString());


    //   return 0;
    // }
  };
}  // namespace dautils