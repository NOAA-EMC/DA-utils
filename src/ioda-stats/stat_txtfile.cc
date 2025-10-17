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
                     << std::left << std::setw(15) << "Group"
                     << std::left << std::setw(5) << "Chan"
                     << std::left << std::setw(12) << "Use"
                     << std::left << std::setw(8) << "Stat"
                     << std::endl;
            txtFile_ << std::setfill('-') << std::setw(90) << "-" << std::endl;
            txtFile_ << std::setfill(' '); // reset fill character
        } else {
            // use a header that allows for vertical bins
            txtFile_ << std::left << std::setw(25) << "Obs Space"
                     << std::left << std::setw(25) << "Variable"
                     << std::left << std::setw(15) << "Group"
                     << std::left << std::setw(12) << "Use"
                     << std::left << std::setw(8) << "Stat"
                     << std::left << std::setw(3) << " | "
                     << std::right << std::setw(8) << "All Bins";
            for (const auto& bin : zBins) {
            txtFile_ << std::left << std::setw(3) << " | "
                     << std::right << std::setw(8) << bin;
            }
            txtFile_ << std::endl;
            txtFile_ << std::setfill('-') << std::setw(96 + zBins.size() * 11) << "-" << std::endl;
            txtFile_ << std::setfill(' '); // reset fill character
        }
        return 0;
    }

    int StatTxtFile::writeTxtStat(const std::string &obsSpaceName, const std::string &variable,
                                  const int &ch, const std::string &group, const std::string &assim,
                                  const std::string &statname, const std::vector<std::vector<float>> &values) {
        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error writing to file" << std::endl;
            return -1;
        }
        txtFile_ << std::left << std::setw(25) << obsSpaceName
                << std::left << std::setw(25) << variable
                << std::left << std::setw(15) << group;
        if (ch >= 0) {
            txtFile_ << std::left << std::setw(5) << ch;
        }
        txtFile_ << std::left << std::setw(12) << assim
                << std::left << std::setw(8) << statname;
        if (ch < 0) {
            for (const auto& val : values) {
                if (val[0] == fillVal_) {
                    txtFile_ << std::left << std::setw(3) << " | "
                             << std::right << std::setw(8) << "NaN";
                } else {
                    txtFile_ << std::left << std::setw(3) << " | "
                             << std::right << std::setw(8) << val[0];
                }
            }
        } else {
            for (const auto& val : values[0]) {
                txtFile_ << std::right << std::setw(8) << val;
            }
        }
        txtFile_ << std::endl;
        return 0;
    }
    int StatTxtFile::writeTxtStat(const std::string &obsSpaceName, const std::string &variable,
                                  const int &ch, const std::string &group, const std::string &assim,
                                  const std::string &statname, const std::vector<std::vector<int>> &values) {
        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error writing to file" << std::endl;
            return -1;
        }
        txtFile_ << std::left << std::setw(25) << obsSpaceName
                << std::left << std::setw(25) << variable
                << std::left << std::setw(15) << group;
        if (ch >= 0) {
            txtFile_ << std::left << std::setw(5) << ch;
        }
        txtFile_ << std::left << std::setw(12) << assim
                << std::left << std::setw(8) << statname;
        if (ch < 0) {
            for (const auto& val : values) {
                std::cout << val << std::endl;
                txtFile_ << std::left << std::setw(3) << " | "
                         << std::right << std::setw(8) << val[0];
            }
        } else {
            for (const auto& val : values[0]) {
                txtFile_ << std::right << std::setw(8) << val;
            }
        }
        txtFile_ << std::endl;
        return 0;
    }
    int StatTxtFile::closeFile() {
        if (txtFile_.is_open()) {
            txtFile_.close();
        }
        return 0;
    }
}