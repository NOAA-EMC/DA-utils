#pragma once

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/TimeWindow.h"
#include "oops/util/missingValues.h"

namespace dautils {
  class StatTxtFile {
    public:
      float fillVal_ = util::missingValue<float>();

      int initializeTxtFile(const std::string &filename, const util::TimeWindow &timeWindow,
                            const std::string &obsSpaceName, bool hasChannels = false, std::vector<std::string> zBins = {});
      int writeData(const std::vector<float> &data);
      int closeFile();

    private:
      std::ofstream txtFile_;
  };
}
