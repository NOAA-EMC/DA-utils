#include <mpi.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Group.h"
#include "ioda/ObsDataIoParameters.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/PostProcessor.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/TimeWindow.h"

#include "calc_iodastats.h"
#include "stat_ncfile.h"
#include "stat_txtfile.h"

void dautils::CalcIodaStats::run() {
    // Main driver code for calculating observation space
    // statistics from IODA files, and depending on configuration,
    // writing them to a new IODA file for concatenation and future use/plotting
    // or writing to ASCII files for quick viewing.

    // define the time window
    const eckit::LocalConfiguration timeWindowConf(config_, "time window");
    const util::TimeWindow timeWindow(timeWindowConf);

    // get the list of obs spaces to process
    std::vector<eckit::LocalConfiguration> obsSpaces;
    config_.get("obs spaces", obsSpaces);

    // get the communicator for just me
    const eckit::mpi::Comm & mycomm = oops::mpi::myself();

    // get MPI comm size
    int nprocs = comm_.size();

    // loop over obs spaces, distributing over MPI ranks
    for (int i = comm_.rank(); i < obsSpaces.size(); i+= nprocs) {
        // Now process each observation space in turn
        // get the configuration for this obs space
        auto obsSpace = obsSpaces[i];
        eckit::LocalConfiguration obsConfig(obsSpace, "obs space");
        std::string obsSpaceName;
        obsConfig.get("name", obsSpaceName);

        // open the IODA file
        std::string obsFile;
        obsConfig.get("obsdatain.engine.obsfile", obsFile);
        oops::Log::info() << "Processing " << obsSpaceName << ":" << obsFile << std::endl;
        ioda::ObsSpace ospace(obsConfig, mycomm, timeWindow, mycomm);
        const size_t nlocs = ospace.nlocs();
        oops::Log::info() << obsSpaceName << ": nlocs =" << nlocs << std::endl;
        if (nlocs == 0) {
            oops::Log::info() << "ObsSpace is empty, skipping " << obsSpaceName << ":" << obsFile << std::endl;
            continue;
        }

        // get the list of variables (and channels if applicable) to process
        std::vector<std::string> variables;
        std::vector<int> channels;
        obsSpace.get("variables", variables);
        if (obsSpace.has("channels")) {
            obsSpace.get("channels", channels);
        }

        // channels only works if there is one variable, so need to check this
        if (variables.size() > 1 && !channels.empty()) {
            throw eckit::Exception("Cannot use channels with multiple variables.");
        }

        // get the lists of everything to process/compute
        std::vector<std::string> groups;
        std::vector<std::string> stats;
        std::vector<std::string> qcgroups;
        std::vector<std::string> errorGroups;
        std::vector<eckit::LocalConfiguration> domains;

        obsSpace.get("groups to process", groups);
        obsSpace.get("qc groups", qcgroups);
        if (obsSpace.has("error groups")) {
            obsSpace.get("error groups", errorGroups);
            if (errorGroups.size() != groups.size()) {
                throw eckit::Exception("If error groups are provided, there must be one for each group.");
            }
        } else {
            // if no error groups provided, create empty strings for each group
            errorGroups.resize(groups.size(), "");
        }
        obsSpace.get("statistics to compute", stats);
        std::vector<std::string> qccategories = obsSpace.getStringVector("qc categories", {"all"});

        obsSpace.get("domains to process", domains);
        // loop over all domains and get their definitions
        std::vector<std::string> domainNames;
        std::vector<std::string> domainMaskVar1, domainMaskVar2, domainMaskVar3;
        std::vector<std::vector<float>> domainMaskVals1, domainMaskVals2, domainMaskVals3;
        for (int idom = 0; idom < domains.size(); idom++ ) {
            eckit::LocalConfiguration domainConf(domains[idom], "domain");
            std::string domName;
            domainConf.get("name", domName);
            domainNames.push_back(domName);
            std::string maskVar1, maskVar2, maskVar3;
            std::vector<float> maskVals1, maskVals2, maskVals3;
            if (domainConf.has("first mask variable")) {
                domainConf.get("first mask variable", maskVar1);
                domainConf.get("first mask range", maskVals1);
            }
            if (domainConf.has("second mask variable")) {
                domainConf.get("second mask variable", maskVar2);
                domainConf.get("second mask range", maskVals2);
            }
            if (domainConf.has("third mask variable")) {
                domainConf.get("third mask variable", maskVar3);
                domainConf.get("third mask range", maskVals3);
            }
            domainMaskVar1.push_back(maskVar1);
            domainMaskVar2.push_back(maskVar2);
            domainMaskVar3.push_back(maskVar3);
            domainMaskVals1.push_back(maskVals1);
            domainMaskVals2.push_back(maskVals2);
            domainMaskVals3.push_back(maskVals3);
        }

        // determine if we are doing regular binning for this obs space
        int nbins_x = 0;
        int nbins_y = 0;
        std::vector<float> bin_lons_centers;
        std::vector<float> bin_lats_centers;
        std::vector<std::string> zBinNames;
        std::vector<std::string> zBinMaskVar;
        std::vector<std::vector<float>> zBinMaskVals;
        if (obsSpace.has("regular grid binning")) {
            eckit::LocalConfiguration binConfig;
            obsSpace.get("regular grid binning", binConfig);
            float binsize;
            binConfig.get("bin size in degrees", binsize);
            nbins_x = int(360.0 / binsize);
            nbins_y = int(180.0 / binsize);
            
            // Compute bin centers for lat/lon
            // Use -180 to 180 longitude range (standard convention)
            float dx = 360.0 / float(nbins_x);
            float dy = 180.0 / float(nbins_y);
            
            // Allocate space for 2D arrays (flattened in row-major order)
            bin_lons_centers.resize(nbins_y * nbins_x);
            bin_lats_centers.resize(nbins_y * nbins_x);
            
            // Compute centers
            for (int iy = 0; iy < nbins_y; iy++) {
                float lat_min = -90.0 + iy * dy;
                float lat_center = lat_min + dy / 2.0;
                for (int ix = 0; ix < nbins_x; ix++) {
                    float lon_min = -180.0 + ix * dx;
                    float lon_center = lon_min + dx / 2.0;
                    int idx = iy * nbins_x + ix;
                    bin_lats_centers[idx] = lat_center;
                    bin_lons_centers[idx] = lon_center;
                }
            }
            
            std::vector<eckit::LocalConfiguration> zBins;
            if (binConfig.has("vertical bins")) {
                binConfig.get("vertical bins", zBins);
                for (int idom = 0; idom < zBins.size(); idom++ ) {
                    auto zBin = zBins[idom];
                    eckit::LocalConfiguration binConf(zBin, "vertical bin");
                    std::string binname, maskvar;
                    std::vector<float> maskvals;
                    binConf.get("name", binname);
                    binConf.get("mask variable", maskvar);
                    binConf.get("mask range", maskvals);
                    zBinNames.push_back(binname);
                    zBinMaskVar.push_back(maskvar);
                    zBinMaskVals.push_back(maskvals);
                }
            }
        }

        // assert that the QC groups list is the same size as groups
        assert(groups.size() == qcgroups.size());

        // if the zBins are empty, create a dummy one for 'all'
        if (zBinNames.size() == 0){
            zBinNames.push_back("all");
            zBinMaskVar.push_back("latitude");
            zBinMaskVals.push_back(std::vector<float> {-90.0f, 90.0f});
        }

        // initialize netCDF output file for writing
        std::string outncfile;
        obsSpace.get("output file", outncfile);
        StatNcFile statncfile;
        statncfile.initializeNcfile(outncfile, timeWindow, variables, channels, groups,
                                    stats, domainNames, nbins_x, nbins_y, zBinNames,
                                    bin_lons_centers, bin_lats_centers);
        
        // if desired, get ranges of vertical bins for ascii reporting
        std::vector<std::string> asciiZBins;
        std::vector<float> asciiZBinRanges;
        std::string asciiZBinVar;
        if (obsSpace.has("ascii vertical bins")) {
            eckit::LocalConfiguration asciiBinConfig;
            obsSpace.get("ascii vertical bins", asciiBinConfig);
            asciiBinConfig.get("vertical bin names", asciiZBins);
            asciiBinConfig.get("vertical bin ranges", asciiZBinRanges);
            asciiBinConfig.get("vertical bin variable", asciiZBinVar);
            if (asciiZBins.size() != asciiZBinRanges.size()) {
                throw eckit::Exception("ascii vertical bin names and ranges must be the same size");
            }
        } 

        // initialize ASCII output for writing
        std::string outasciifile = obsSpaceName + "_ioda_stats.txt";
        if (obsSpace.has("output ascii file")) {
            obsSpace.get("output ascii file", outasciifile);
        }
        StatTxtFile stattxtfile;
        stattxtfile.initializeTxtFile(outasciifile, timeWindow, obsSpaceName, obsSpace.has("channels"), asciiZBins);
        
        // first let us loop over the ascii vertical bins if they are defined, and a total if not

    } // end of obs space loop
} // end of run