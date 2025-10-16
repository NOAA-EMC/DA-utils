#include "stat_txtfile.h"

namespace dautils {

    int StatTxtFile::initializeTxtFile(const std::string &filename, const util::TimeWindow &timeWindow,
                                    const std::string &obsSpaceName, bool hasChannels, std::vector<std::string> zBins) {
        txtFile_.open(filename);
        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error opening file: " << filename << std::endl;
            return -1;
        }
        txtFile_ << "==================================================" << std::endl;
        txtFile_ << "Observation Space: " << obsSpaceName << std::endl;
        txtFile_ << "Analysis Time: " << timeWindow.midpoint().toString() << std::endl;
        txtFile_ << "--------------------------------------------------" << std::endl;
        if (hasChannels) {
            // use a different header if channels are present
            txtFile_ << std::left << std::setw(25) << "Obs Space"
                     << std::left << std::setw(25) << "Variable"
                     << std::left << std::setw(5) << "Chan"
                     << std::left << std::setw(5) << "Loop"
                     << std::left << std::setw(12) << "Use"
                     << std::endl;
        } else {
            // use a header that allows for vertical bins
            txtFile_ << std::left << std::setw(25) << "Obs Space"
                     << std::left << std::setw(25) << "Variable"
                     << std::left << std::setw(5) << "Loop"
                     << std::left << std::setw(12) << "Use";
            for (const auto& bin : zBins) {
            txtFile_ << std::left << std::setw(10) << bin;
            }
            txtFile_ << std::endl;
        }
        return 0;
    }
}