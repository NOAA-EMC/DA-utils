

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "ioda_stats_driver.h"

void dautils::IodaStatsDriver::run() {
    // Main driver function of the
    // IodaStats application.
    // First, define the time window
    const eckit::LocalConfiguration timeWindowConf(fullConfig, "time window");
    const util::TimeWindow timeWindow(timeWindowConf);

    // get the list of obs spaces to process
    std::vector<eckit::LocalConfiguration> obsSpaces;
    fullConfig.get("obs spaces", obsSpaces);

    // get the communicator for just me
    const eckit::mpi::Comm & mycomm = oops::mpi::myself();

    // get MPI comm size
    int nprocs = getComm().size();
}