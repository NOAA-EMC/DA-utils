#pragma once

#include <netcdf>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "oops/util/missingValues.h"

namespace dautils {
  class StatFile {
    public:
    float fillVal_ = util::missingValue<float>();
    int initializeNcfile(const std::string filename, const util::TimeWindow timeWindow,
                      std::vector<std::string> variables, std::vector<int> channels,
                      std::vector<std::string> groups, std::vector<std::string> stats,
                      std::vector<std::string> domainNames,
                      int nbins_x, int nbins_y,
                      std::vector<std::string> bins_z,
                      std::vector<float> bin_lons = std::vector<float>(),
                      std::vector<float> bin_lats = std::vector<float>()) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::replace);
      oops::Log::info() << "Opening " << filename << " for writing..." << std::endl;
      // create an unlimited time dimension
      netCDF::NcDim tDim = ncFile.addDim("analysisCycle");
      // create domain dimension
      int ndomains = domainNames.size() + 1;
      netCDF::NcDim dDim = ncFile.addDim("Domain", ndomains);
      // vector of dimensions
      std::vector<netCDF::NcDim> domainDimVector, binningDimVector;
      domainDimVector.push_back(tDim);
      domainDimVector.push_back(dDim);
      // if channel is not empty, create a channel dimension
      netCDF::NcDim cDim;
      if (!channels.empty()) {
        cDim = ncFile.addDim("Channel", channels.size());
        domainDimVector.push_back(cDim);
      }
      // if nbins_x or nbins_y are nonzero, make them dimensions too
      netCDF::NcDim xDim, yDim, zDim;
      if (nbins_x > 0 || nbins_y > 0) {
        binningDimVector.push_back(tDim);
        if (!channels.empty()) {
          binningDimVector.push_back(cDim);
        } else {
          if (bins_z.size() > 0) {
            zDim = ncFile.addDim("binsZDim", bins_z.size());
            binningDimVector.push_back(zDim);
          }
        }
      }

      if (nbins_y > 0) {
        yDim = ncFile.addDim("binsYDim", nbins_y);
        binningDimVector.push_back(yDim);
      }
      if (nbins_x > 0) {
        xDim = ncFile.addDim("binsXDim", nbins_x);
        binningDimVector.push_back(xDim);
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
      domainNames.push_back("Global");
      for (int idom = 0; idom < ndomains; idom++) {
        std::vector<size_t> idxdom;
        idxdom.push_back(idom);
        domain.putVar(idxdom, domainNames[idom]);
      }

      // create vertical bin variable
      if (channels.empty() && bins_z.size() > 0 && (nbins_x > 0 || nbins_y > 0)) {
        netCDF::NcVar zbins = ncFile.addVar("verticalBin", netCDF::ncString, zDim);
        for (int ibin = 0; ibin < bins_z.size(); ibin++) {
          std::vector<size_t> idxbin;
          idxbin.push_back(ibin);
          zbins.putVar(idxbin, bins_z[ibin]);
        }
      }

      // loop over group, then variables, then stats to create byDomains/group/var/stat in file
      netCDF::NcGroup domaingroup = ncFile.addGroup("byDomains");
      for (int g = 0; g < groups.size(); g++) {
        // create group group
        netCDF::NcGroup group = domaingroup.addGroup(groups[g]);
        // loop over variables
        for (int var = 0; var < variables.size(); var++) {
          // create variable group
          netCDF::NcGroup group2 = group.addGroup(variables[var]);
          // loop over statistics to write out
          for (int s = 0; s < stats.size(); s++) {
            netCDF::NcVar varout;
            if (stats[s] == "count") {
              varout = group2.addVar(stats[s], netCDF::ncInt, domainDimVector);
            } else {
              varout = group2.addVar(stats[s], netCDF::ncFloat, domainDimVector);
              varout.putAtt("_FillValue", netCDF::ncFloat, fillVal_);
            }
          }
        }
      }

      // loop over group, then variables, then stats to create griddedBins/group/var/stat in file
      if (nbins_x > 0 || nbins_y > 0) {
        netCDF::NcGroup bingroup = ncFile.addGroup("griddedBins");
        
        // Write latitude and longitude coordinate arrays if provided
        // Both nbins_x and nbins_y must be > 0 for 2D gridded bins
        if (!bin_lats.empty() && !bin_lons.empty() && nbins_x > 0 && nbins_y > 0) {
          // Create 2D dimensions for lat/lon arrays
          std::vector<netCDF::NcDim> coordDims;
          coordDims.push_back(yDim);
          coordDims.push_back(xDim);
          
          // Create and write latitude variable
          netCDF::NcVar latVar = bingroup.addVar("latitude", netCDF::ncFloat, coordDims);
          latVar.putAtt("units", "degrees_north");
          latVar.putAtt("long_name", "latitude of bin centers");
          
          // Create and write longitude variable
          netCDF::NcVar lonVar = bingroup.addVar("longitude", netCDF::ncFloat, coordDims);
          lonVar.putAtt("units", "degrees_east");
          lonVar.putAtt("long_name", "longitude of bin centers");
          
          // Write the data
          std::vector<size_t> start = {0, 0};
          std::vector<size_t> count = {static_cast<size_t>(nbins_y), static_cast<size_t>(nbins_x)};
          latVar.putVar(start, count, bin_lats.data());
          lonVar.putVar(start, count, bin_lons.data());
        }
        
        for (int g = 0; g < groups.size(); g++) {
          // create group group
          netCDF::NcGroup group = bingroup.addGroup(groups[g]);
          // loop over variables
          for (int var = 0; var < variables.size(); var++) {
            // create variable group
            netCDF::NcGroup group2 = group.addGroup(variables[var]);
            // loop over statistics to write out
            for (int s = 0; s < stats.size(); s++) {
              netCDF::NcVar varout;
              if (stats[s] == "count") {
                varout = group2.addVar(stats[s], netCDF::ncInt, binningDimVector);
              } else {
                varout = group2.addVar(stats[s], netCDF::ncFloat, binningDimVector);
                varout.putAtt("_FillValue", netCDF::ncFloat, fillVal_);
              }
            }
          }
        }
      }

      oops::Log::info() << "Output file " << filename << " has been created." << std::endl;
      return 0;
    };

    // Overloaded write methods
    int writeByDomains(const std::string filename, const std::string group, const std::string variable,
              const std::string stat, const int idom, const std::vector<int> intvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup domaingroup = ncFile.getGroup("byDomains");
      netCDF::NcGroup outgroup1 = domaingroup.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      idxout.push_back(0);
      idxout.push_back(idom);
      if (intvals.size() > 1) {
        idxout.push_back(0);
        std::vector<size_t> countout;
        countout.push_back(1);
        countout.push_back(1);
        countout.push_back(intvals.size());
        outvar.putVar(idxout, countout, intvals.data());
      } else {
        outvar.putVar(idxout, intvals[0]);
      }
      return 0;
    };

    int writeByDomains(const std::string filename, const std::string group, const std::string variable,
              const std::string stat, const int idom, const std::vector<float> floatvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup domaingroup = ncFile.getGroup("byDomains");
      netCDF::NcGroup outgroup1 = domaingroup.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      idxout.push_back(0);
      idxout.push_back(idom);
      if (floatvals.size() > 1) {
        idxout.push_back(0);
        std::vector<size_t> countout;
        countout.push_back(1);
        countout.push_back(1);
        countout.push_back(floatvals.size());
        outvar.putVar(idxout, countout, floatvals.data());
      } else {
        outvar.putVar(idxout, floatvals[0]);
      }
      return 0;
    };

    int writeByBins(const std::string filename, const std::string group, const std::string variable,
              const std::string stat, const int ibin, const int ny, const int nx,
              const std::vector<std::vector<int>> intvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup bingroup = ncFile.getGroup("griddedBins");
      netCDF::NcGroup outgroup1 = bingroup.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      std::vector<size_t> countout;
      int idxtmp = 0;
      idxout.push_back(0);
      idxout.push_back(ibin);
      idxout.push_back(0);
      idxout.push_back(0);
      countout.push_back(1);
      countout.push_back(1);
      countout.push_back(1);
      countout.push_back(nx);
      for (auto i: intvals) {
        idxout[2] = idxtmp;
        outvar.putVar(idxout, countout, i.data());
        ++idxtmp;
      }
      return 0;
    };

    int writeByBins(const std::string filename, const std::string group, const std::string variable,
              const std::string stat, const int ibin, const int ny, const int nx,
              const std::vector<std::vector<float>> floatvals) {
      netCDF::NcFile ncFile(filename, netCDF::NcFile::write);
      netCDF::NcGroup bingroup = ncFile.getGroup("griddedBins");
      netCDF::NcGroup outgroup1 = bingroup.getGroup(group);
      netCDF::NcGroup outgroup2 = outgroup1.getGroup(variable);
      netCDF::NcVar outvar = outgroup2.getVar(stat);
      std::vector<size_t> idxout;
      std::vector<size_t> countout;
      int idxtmp = 0;
      idxout.push_back(0);
      idxout.push_back(ibin);
      idxout.push_back(0);
      idxout.push_back(0);
      countout.push_back(1);
      countout.push_back(1);
      countout.push_back(1);
      countout.push_back(nx);
      for (auto i: floatvals) {
        idxout[2] = idxtmp;
        outvar.putVar(idxout, countout, i.data());
        ++idxtmp;
      }
      return 0;
    };

  };
}  // namespace dautils
