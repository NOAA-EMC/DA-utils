#pragma once

#include <netcdf>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "oops/util/missingValues.h"

namespace dautils {
  class IodaStatsFile {
    public:
      float fillVal_ = util::missingValue<float>();

      int initializeNcfile(const std::string obspace, const std::string filename,
                           const util::TimeWindow & timeWindow,
                           std::vector<std::string> variables, std::vector<int> channels,
                           std::vector<std::string> groups, std::vector<std::string> stats,
                           std::vector<std::string> domainNames,
                           int nbins_x, int nbins_y,
                           std::vector<std::string> bins_z);
  };
}  // namespace dautils
