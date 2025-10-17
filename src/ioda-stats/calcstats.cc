#include "calcstats.h"

#include <cmath>

namespace dautils {

std::vector<int> getObsCount(const std::vector<float> &data,
                             const std::vector<int> &qcvals,
                             const std::vector<int> &channels,
                             const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<int> counts;
  if (channels.empty()) {
    int count(0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count += 1;
      }
    }
    counts.push_back(count);
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count(0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count += 1;
        }
      }
      counts.push_back(count);
    }
  }
  return counts;
}

// -----------------------------------------------------------------------------
std::vector<float> getMean(const std::vector<float> &data,
                           const std::vector<int> &qcvals,
                           const std::vector<int> &channels,
                           const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<float> means;
  if (channels.empty()) {
    int count(0);
    float mean(0.0);
    double sum(0.0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count += 1;
        sum += data[i];
      }
    }
    if (count > 0) {
      mean = sum / count;
    } else {
      mean = fillVal;
    }
    means.push_back(mean);
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count(0);
      float mean(0.0);
      double sum(0.0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count += 1;
          sum += data[ii];
        }
      }
      if (count > 0) {
        mean = sum / count;
      } else {
        mean = fillVal;
      }
      means.push_back(mean);
    }
  }
  return means;
}

// -----------------------------------------------------------------------------
std::vector<float> getRMS(const std::vector<float> &data,
                          const std::vector<int> &qcvals,
                          const std::vector<int> &channels,
                          const std::vector<int> &mask) {
  float fillVal = util::missingValue<float>();
  std::vector<float> rmsvals;
  if (channels.empty()) {
    int count(0);
    float rms(0.0);
    double sum(0.0);
    for (size_t i = 0; i < data.size(); ++i) {
      if (data[i] != fillVal && qcvals[i] == 0 && mask[i] == 0) {
        count += 1;
        sum += pow(data[i], 2);
      }
    }
    if (count > 0) {
      rms = sqrt(sum / count);
    } else {
      rms = fillVal;
    }
    rmsvals.push_back(rms);
  } else {
    int _nlocs = data.size() / channels.size();
    for (int ch = 0; ch < channels.size(); ch++) {
      int count(0);
      float rms(0.0);
      double sum(0.0);
      for (size_t i = 0; i < _nlocs; ++i) {
        int ii = ch + (i * channels.size());
        if (data[ii] != fillVal && qcvals[ii] == 0 && mask[i] == 0) {
          count += 1;
          sum += pow(data[ii], 2);
        }
      }
      if (count > 0) {
        rms = sqrt(sum / count);
      } else {
        rms = fillVal;
      }
      rmsvals.push_back(rms);
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
  for (int i = 0; i < maskvalues.size(); i++) {
    if (maskvalues[i] < minval || maskvalues[i] >= maxval) {
      updatedMask[i] = 1;
    }
  }
  return updatedMask;
}

}  // namespace dautils