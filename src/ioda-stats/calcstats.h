#pragma once

#include <vector>

namespace dautils {
  
  // Statistical utility functions for observation data processing

  std::vector<std::vector<int>> getObsCount(const std::vector<float> &data,
                                            const std::vector<int> &qcvals,
                                            const std::vector<int> &channels,
                                            const std::vector<int> &mask);

  std::vector<std::vector<float>> getMean(const std::vector<float> &data,
                                          const std::vector<int> &qcvals,
                                          const std::vector<int> &channels,
                                          const std::vector<int> &mask);

  std::vector<std::vector<float>> getRMS(const std::vector<float> &data,
                                         const std::vector<int> &qcvals,
                                         const std::vector<int> &channels,
                                         const std::vector<int> &mask);

  std::vector<int> update_mask(std::vector<float> maskvalues, 
                               float minval, 
                               float maxval, 
                               const std::vector<int>& inputMask);
                               
}  // namespace dautils
