#include "calc_iodastats.h"

void dautils::CalcIodaStats::run() {
    // Main driver code for calculating observation space
    // statistics from IODA files, and depending on configuration,
    // writing them to a new IODA file for concatenation and future use/plotting
    // or writing to ASCII files for quick viewing.

    // First, get the full configuration passed to this class
    const eckit::LocalConfiguration & fullConfig = this->getConfig();
    // define the time window
    const eckit::LocalConfiguration timeWindowConf(fullConfig, "time window");
    const util::TimeWindow timeWindow(timeWindowConf);
}