#include "calcstats.h"

#include <cmath>

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace dautils {

std::vector<std::vector<int>> getObsCount(const std::vector<float> &data,
                                          const std::vector<int> &qcvals,
                                          const std::vector<int> &channels,
                                          const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<std::vector<int>> counts;
  if (channels.empty()) {
    int count_assim(0), count_mon(0), count_rej(0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count_assim += 1;
      } else if (data[i] != fillVal && qcvals[i] == 1 && mask[i] == 0) {
        count_mon += 1;
      } else if (data[i] != fillVal && qcvals[i] > 1 && mask[i] == 0) {
        count_rej += 1;
      }
    }
    counts.push_back({count_assim, count_mon, count_rej});
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count_assim(0), count_mon(0), count_rej(0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count_assim += 1;
        } else if (data[ii] != fillVal && qcvals[ii] == 1 && mask[i] == 0) {
          count_mon += 1;
        } else if (data[ii] != fillVal && qcvals[ii] > 1 && mask[i] == 0) {
          count_rej += 1;
        }
      }
      counts.push_back({count_assim, count_mon, count_rej});
    }
  }
  return counts;
}

// -----------------------------------------------------------------------------
std::vector<std::vector<float>> getMean(const std::vector<float> &data,
                                        const std::vector<int> &qcvals,
                                        const std::vector<int> &channels,
                                        const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<std::vector<float>> means;
  if (channels.empty()) {
    int count_assim(0), count_mon(0), count_rej(0);
    float mean_assim(0.0), mean_mon(0.0), mean_rej(0.0);
    double sum_assim(0.0), sum_mon(0.0), sum_rej(0.0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count_assim += 1;
        sum_assim += data[i];
      } else if (data[i] != fillVal && qcvals[i] == 1 && mask[i] == 0) {
        count_mon += 1;
        sum_mon += data[i];
      } else if (data[i] != fillVal && qcvals[i] > 1 && mask[i] == 0) {
        count_rej += 1;
        sum_rej += data[i];
      }
    }
    if (count_assim > 0) {
      mean_assim = sum_assim / count_assim;
    } else {
      mean_assim = fillVal;
    }
    if (count_mon > 0) {
      mean_mon = sum_mon / count_mon;
    } else {
      mean_mon = fillVal;
    }
    if (count_rej > 0) {
      mean_rej = sum_rej / count_rej;
    } else {
      mean_rej = fillVal;
    }
    means.push_back({mean_assim, mean_mon, mean_rej});
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count_assim(0), count_mon(0), count_rej(0);
      float mean_assim(0.0), mean_mon(0.0), mean_rej(0.0);
      double sum_assim(0.0), sum_mon(0.0), sum_rej(0.0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count_assim += 1;
          sum_assim += data[ii];
        } else if (data[ii] != fillVal && qcvals[ii] == 1 && mask[i] == 0) {
          count_mon += 1;
          sum_mon += data[ii];
        } else if (data[ii] != fillVal && qcvals[ii] > 1 && mask[i] == 0) {
          count_rej += 1;
          sum_rej += data[ii];
        }
      }
      if (count_assim > 0) {
        mean_assim = sum_assim / count_assim;
      } else {
        mean_assim = fillVal;
      }
      if (count_mon > 0) {
        mean_mon = sum_mon / count_mon;
      } else {
        mean_mon = fillVal;
      }
      if (count_rej > 0) {
        mean_rej = sum_rej / count_rej;
      } else {
        mean_rej = fillVal;
      }
      means.push_back({mean_assim, mean_mon, mean_rej});
    }
  }
  return means;
}

// -----------------------------------------------------------------------------
std::vector<std::vector<float>> getRMS(const std::vector<float> &data,
                                       const std::vector<int> &qcvals,
                                       const std::vector<int> &channels,
                                       const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<std::vector<float>> rmsvals;
  if (channels.empty()) {
    int count_assim(0), count_mon(0), count_rej(0);
    float rms_assim(0.0), rms_mon(0.0), rms_rej(0.0);
    double sum_assim(0.0), sum_mon(0.0), sum_rej(0.0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count_assim += 1;
        sum_assim += pow(data[i], 2);
      } else if (data[i] != fillVal && qcvals[i] == 1 && mask[i] == 0) {
        count_mon += 1;
        sum_mon += pow(data[i], 2);
      } else if (data[i] != fillVal && qcvals[i] > 1 && mask[i] == 0) {
        count_rej += 1;
        sum_rej += pow(data[i], 2);
      }
    }
    if (count_assim > 0) {
      rms_assim = sqrt(sum_assim / count_assim);
    } else {
      rms_assim = fillVal;
    }
    if (count_mon > 0) {
      rms_mon = sqrt(sum_mon / count_mon);
    } else {
      rms_mon = fillVal;
    }
    if (count_rej > 0) {
      rms_rej = sqrt(sum_rej / count_rej);
    } else {
      rms_rej = fillVal;
    }
    rmsvals.push_back({rms_assim, rms_mon, rms_rej});
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count_assim(0), count_mon(0), count_rej(0);
      float rms_assim(0.0), rms_mon(0.0), rms_rej(0.0);
      double sum_assim(0.0), sum_mon(0.0), sum_rej(0.0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count_assim += 1;
          sum_assim += pow(data[ii], 2);
        } else if (data[ii] != fillVal && qcvals[ii] == 1 && mask[i] == 0) {
          count_mon += 1;
          sum_mon += pow(data[ii], 2);
        } else if (data[ii] != fillVal && qcvals[ii] > 1 && mask[i] == 0) {
          count_rej += 1;
          sum_rej += pow(data[ii], 2);
        }
      }
      if (count_assim > 0) {
        rms_assim = sqrt(sum_assim / count_assim);
      } else {
        rms_assim = fillVal;
      }
      if (count_mon > 0) {
        rms_mon = sqrt(sum_mon / count_mon);
      } else {
        rms_mon = fillVal;
      }
      if (count_rej > 0) {
        rms_rej = sqrt(sum_rej / count_rej);
      } else {
        rms_rej = fillVal;
      }
      rmsvals.push_back({rms_assim, rms_mon, rms_rej});
    }
  }
  return rmsvals;
}

// -----------------------------------------------------------------------------
// update mask provided values and a min/max value to mask outside of
std::vector<int> update_mask(std::vector<float> maskvalues, 
                             float minval, 
                             float maxval, 
                             const std::vector<int>& inputMask) {
  std::vector<int> updatedMask = inputMask;
  for (size_t i = 0; i < maskvalues.size(); i++) {
    if (maskvalues[i] < minval || maskvalues[i] >= maxval) {
      updatedMask[i] = 1;
    }
  }
  return updatedMask;
}

}  // namespace dautils