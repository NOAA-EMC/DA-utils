#pragma once

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace dautils {
  class ObsStats {
    public:
    float fillVal_ = util::missingValue<float>();
    std::vector<int> getObsCount(const std::vector<float> &data,
                                 const std::vector<int> &qcvals,
                                 const std::vector<int> &channels,
                                 const std::vector<int> &mask) {
      std::vector<int> counts;
      if (channels.empty()) {
        int count(0);
        for (size_t i = 0; i < data.size(); ++i) {
          if (data[i] != fillVal_ && qcvals[i] == 0 && mask[i] == 0) {
            count += 1;
          }
        }
        counts.push_back(count);
      } else {
        for (int ch = 0; ch < channels.size(); ch++) {
          oops::Log::info() << ch << std::endl;
        }
      }
      return counts;
    }
    // -----------------------------------------------------------------------------
    std::vector<float> getMean(const std::vector<float> &data,
                               const std::vector<int> &qcvals,
                               const std::vector<int> &channels,
                               const std::vector<int> &mask) {
      std::vector<float> means;
      if (channels.empty()) {
        int count(0);
        float mean(0.0);
        for (size_t i = 0; i < data.size(); ++i) {
          if (data[i] != fillVal_ && qcvals[i] == 0 && mask[i] == 0) {
            count += 1;
            mean += data[i];
          }
        }
        if (count > 0) {
          mean = mean / count;
        }
        means.push_back(mean);
      } else {
        for (int ch = 0; ch < channels.size(); ch++) {
          oops::Log::info() << ch << std::endl;
        }
      }
    //   if (channels.empty()) {
    //     int count(0);
    //     for (size_t i = 0; i < data.size(); ++i) {
    //       if (data[i] != fillVal_ && qcvals[i] == 0) {
    //         count += 1;
    //       }
    //     }
    //     counts.push_back(count);
    //   } else {
    //     for (int ch = 0; ch < channels.size(); ch++) {
    //       oops::Log::info() << ch << std::endl;
    //     }
    //   }
      return means;
    }
    // -----------------------------------------------------------------------------
    // update mask provided values and a min/max value to mask outside of
    std::vector<int> update_mask(std::vector<float> maskvalues, float minval, float maxval, const std::vector<int>& inputMask) {
      std::vector<int> updatedMask = inputMask;
      for (int i = 0; i < maskvalues.size(); i++) {
        if (maskvalues[i] < minval || maskvalues[i] > maxval) {
          updatedMask[i] = 1;
        }
      }
      return updatedMask;
    }
    // -----------------------------------------------------------------------------
  };
}  // namespace dautils
