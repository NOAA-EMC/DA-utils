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
                      std::vector<std::string> groups, std::vector<std::string> stats,
                      std::vector<std::string> domainNames) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::replace);
      oops::Log::info() << "Opening " << filename << " for writing..." << std::endl;
      // create an unlimited time dimension
      netCDF::NcDim tDim = ncFile.addDim("analysisCycle");
      // create domain dimension
      int ndomains = domainNames.size() + 1;
      netCDF::NcDim dDim = ncFile.addDim("Domain", ndomains);
      // vector of dimensions
      std::vector<netCDF::NcDim> dimVector;
      dimVector.push_back(tDim);
      dimVector.push_back(dDim);
      // if channel is not empty, create a channel dimension
      netCDF::NcDim cDim;
      if (!channels.empty()) {
        cDim = ncFile.addDim("Channel", channels.size());
        dimVector.push_back(cDim);
      }
      // create validTime variable
      netCDF::NcVar time = ncFile.addVar("validTime", netCDF::ncString, tDim);
      // put the analysis time in the file
      util::DateTime analysisTime = timeWindow.midpoint();
      std::vector<size_t> idxout;
      idxout.push_back(0);
      time.putVar(idxout, analysisTime.toString());

      // create domain variable
      netCDF::NcVar domain = ncFile.addVar("statisticDomain", netCDF::ncString, dDim);
      domainNames.insert(domainNames.begin(), "Global");
      domain.putVar(domainNames.data());

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
              varout = group2.addVar(stats[s], netCDF::ncInt, dimVector);
            } else {
              varout = group2.addVar(stats[s], netCDF::ncFloat, dimVector);
            }
          }
        }
      }
      oops::Log::info() << "Output file " << filename << " has been created." << std::endl;
      return 0;
    };

    // TODO make these into an overloaded method at some point
    int writeInt(const std::string filename, const std::string group, const std::string variable,
                 const std::string stat, const std::vector<int> intvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup outgroup1 = ncFile.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      idxout.push_back(0);
      outvar.putVar(idxout, intvals[0]);
      return 0;
    };

    int writeFloat(const std::string filename, const std::string group, const std::string variable,
                 const std::string stat, const std::vector<float> floatvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup outgroup1 = ncFile.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      idxout.push_back(0);
      outvar.putVar(idxout, floatvals[0]);
      return 0;
    };
  };
}  // namespace dautils
