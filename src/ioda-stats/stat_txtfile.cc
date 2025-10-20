#include "stat_txtfile.h"

namespace dautils {

    int StatTxtFile::initializeTxtFile(const std::string &filename, const util::TimeWindow &timeWindow,
                                    const std::string &obsSpaceName, const int &nlocs, bool hasChannels,
                                    std::vector<std::string> zBins) {
        txtFile_.open(filename);
        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error opening file: " << filename << std::endl;
            return -1;
        }
        txtFile_ << std::setfill('=') << std::setw(90) << "=" << std::endl;
        txtFile_ << std::setfill(' '); // reset fill character
        txtFile_ << "Observation Space: " << obsSpaceName 
                 << "     Analysis Time: " << timeWindow.midpoint().toString()
                 << "     nlocs=" << nlocs << std::endl;
        txtFile_ << std::setfill('-') << std::setw(90) << "-" << std::endl;
        txtFile_ << std::setfill(' '); // reset fill character
        if (hasChannels) {
            // use a different format if channels are present
            txtFile_ << std::left << std::setw(25) << "Obs Space"
                     << std::left << std::setw(25) << "Variable"
                     << std::left << std::setw(15) << "Group"
                     << std::left << std::setw(10) << "Chan"
                     << std::left << std::setw(8) << "Stat"
                     << std::left << std::setw(3) << "   "
                     << std::left << std::setw(11) << "assimilated"
                     << std::left << std::setw(3) << "   "
                     << std::left << std::setw(11) << "monitored"
                     << std::left << std::setw(3) << "   "
                     << std::left << std::setw(11) << "rejected"
                     << std::endl;
            txtFile_ << std::setfill('-') << std::setw(126) << "-" << std::endl;
            txtFile_ << std::setfill(' '); // reset fill character
        } else {

            // use a header that allows for vertical bins
            txtFile_ << std::left << std::setw(25) << "Obs Space"
                     << std::left << std::setw(25) << "Variable"
                     << std::left << std::setw(15) << "Group"
                     << std::left << std::setw(14) << "Use"
                     << std::left << std::setw(8) << "Stat"
                     << std::left << std::setw(3) << " | "
                     << std::right << std::setw(8) << "All Bins";
            for (const auto& bin : zBins) {
            txtFile_ << std::left << std::setw(3) << " | "
                     << std::right << std::setw(8) << bin;
            }
            txtFile_ << std::endl;
            txtFile_ << std::setfill('-') << std::setw(98 + zBins.size() * 11) << "-" << std::endl;
            txtFile_ << std::setfill(' '); // reset fill character
        }
        return 0;
    }

    int StatTxtFile::writeTxtStat(const std::string &obsSpaceName, const std::string &variable,
                                  const std::vector<int> &ch, const std::string &group, const std::string &statname,
                                  const std::vector<std::vector<std::vector<float>>> &values) {

        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error writing to file" << std::endl;
            return -1;
        }
        // determine if we have channels or not
        bool hasChans = false;
        if (ch[0] >= 0) {
            hasChans = true;
        }
        
        std::vector<std::string> useTypes = {"assimilated", "monitored", "rejected"};
        if (hasChans) {
            // for obs spaces with channels, for each statistic,
            // we have rows for each channel and columns of use type
            for (size_t ich = 0; ich < ch.size(); ++ich) {
                int chan = ch[ich];
                txtFile_ << std::left << std::setw(25) << obsSpaceName
                         << std::left << std::setw(25) << variable
                         << std::left << std::setw(15) << group
                         << std::left << std::setw(10) << chan
                         << std::left << std::setw(8) << statname;
                for (size_t iuse = 0; iuse < useTypes.size(); ++iuse) {
                    if (values[iuse][0][ich] == fillVal_) {
                        txtFile_ << std::left << std::setw(3) << " | "
                                 << std::right << std::setw(11) << "NaN";
                    } else {
                        txtFile_ << std::left << std::setw(3) << " | "
                                 << std::right << std::setw(11) << values[iuse][0][ich];
                    }
                }
                txtFile_ << std::endl;
            }
                
        } else {
            // for obs spaces without channels, for each statistic,
            // we have rows for each use type and columns of bins

            // loop over assimilated, monitored, rejected
            for (size_t iuse = 0; iuse < useTypes.size(); ++iuse) {
                std::string assim = useTypes[iuse];
                txtFile_ << std::left << std::setw(25) << obsSpaceName
                        << std::left << std::setw(25) << variable
                        << std::left << std::setw(15) << group
                        << std::left << std::setw(14) << assim
                        << std::left << std::setw(8) << statname;
                for (const auto& val : values[iuse]) {
                    if (val[0] == fillVal_) {
                        txtFile_ << std::left << std::setw(3) << " | "
                                << std::right << std::setw(8) << "NaN";
                    } else {
                        txtFile_ << std::left << std::setw(3) << " | "
                                << std::right << std::setw(8) << val[0];
                    }
                }
                txtFile_ << std::endl;             
            }
        }
        return 0;
    }
    int StatTxtFile::writeTxtStat(const std::string &obsSpaceName, const std::string &variable,
                                  const std::vector<int> &ch, const std::string &group, const std::string &statname,
                                  const std::vector<std::vector<std::vector<int>>> &values) {
        if (!txtFile_.is_open()) {
            oops::Log::error() << "Error writing to file" << std::endl;
            return -1;
        }
        // determine if we have channels or not
        bool hasChans = false;
        if (ch[0] >= 0) {
            hasChans = true;
        }
        
        std::vector<std::string> useTypes = {"assimilated", "monitored", "rejected"};
        if (hasChans) {
            // for obs spaces with channels, for each statistic,
            // we have rows for each channel and columns of use type
            for (size_t ich = 0; ich < ch.size(); ++ich) {
                int chan = ch[ich];
                txtFile_ << std::left << std::setw(25) << obsSpaceName
                         << std::left << std::setw(25) << variable
                         << std::left << std::setw(15) << group
                         << std::left << std::setw(10) << chan
                         << std::left << std::setw(8) << statname;
                for (size_t iuse = 0; iuse < useTypes.size(); ++iuse) {
                    txtFile_ << std::left << std::setw(3) << " | "
                                << std::right << std::setw(11) << values[iuse][0][ich];
                }
                txtFile_ << std::endl;
            }
                
        } else {
            // for obs spaces without channels, for each statistic,
            // we have rows for each use type and columns of bins

            // loop over assimilated, monitored, rejected
            for (size_t iuse = 0; iuse < useTypes.size(); ++iuse) {
                std::string assim = useTypes[iuse];
                txtFile_ << std::left << std::setw(25) << obsSpaceName
                        << std::left << std::setw(25) << variable
                        << std::left << std::setw(15) << group
                        << std::left << std::setw(14) << assim
                        << std::left << std::setw(8) << statname;
                for (const auto& val : values[iuse]) {
                    txtFile_ << std::left << std::setw(3) << " | "
                            << std::right << std::setw(8) << val[0];
                }
                txtFile_ << std::endl;             
            }
        }
        return 0;
    }

    int StatTxtFile::closeFile() {
        if (txtFile_.is_open()) {
            txtFile_.close();
        }
        return 0;
    }
}