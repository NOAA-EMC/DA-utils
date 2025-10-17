#pragma once

#include <netcdf>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "oops/util/missingValues.h"

namespace dautils {
  class StatNcFile {
    public:
      float fillVal_ = util::missingValue<float>();
      std::vector<std::string> use_categories = {"assimilated", "monitored", "rejected"};
      
      int initializeNcfile(const std::string filename, const util::TimeWindow timeWindow,
                          std::vector<std::string> variables, std::vector<int> channels,
                          std::vector<std::string> groups, std::vector<std::string> stats,
                          std::vector<std::string> domainNames,
                          int nbins_x, int nbins_y,
                          std::vector<std::string> bins_z,
                          std::vector<float> bin_lons = std::vector<float>(),
                          std::vector<float> bin_lats = std::vector<float>());

      // Overloaded write methods
      int writeByDomains(const std::string filename, const std::string group, const std::string variable,
                        const std::string stat, const int idom, const std::vector<int> intvals);

      int writeByDomains(const std::string filename, const std::string group, const std::string variable,
                        const std::string stat, const int idom, const std::vector<float> floatvals);

      int writeByBins(const std::string filename, const std::string group, const std::string variable,
                     const std::string stat, const int ibin, const int ny, const int nx,
                     const std::vector<std::vector<int>> intvals);

      int writeByBins(const std::string filename, const std::string group, const std::string variable,
                     const std::string stat, const int ibin, const int ny, const int nx,
                     const std::vector<std::vector<float>> floatvals);
  };
}  // namespace dautils
