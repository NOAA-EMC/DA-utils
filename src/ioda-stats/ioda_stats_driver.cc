

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"

#include "ioda/Engines/EngineUtils.h"
#include "ioda/Group.h"
#include "ioda/ObsDataIoParameters.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "ioda_stats_driver.h"
#include "ioda_stats_file.h"

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

    // loop over all obs spaces
    oops::Log::info() << "IodaStats: Processing " << obsSpaces.size()
                      << " obs spaces on " << nprocs << " processors." << std::endl;
    for (int i = getComm().rank(); i < obsSpaces.size(); i+= nprocs) {
        // get the configuration for this obs space
        auto obsSpace = obsSpaces[i];
        eckit::LocalConfiguration obsConfig(obsSpace, "obs space");

        // open the IODA file
        std::string obsFile;
        std::string name;
        obsConfig.get("name", name);
        obsConfig.get("obsdatain.engine.obsfile", obsFile);
        oops::Log::info() << "IodaStats:" << name << ": Processing " << obsFile << std::endl;
        ioda::ObsSpace ospace(obsConfig, mycomm, timeWindow, mycomm);
        const size_t nlocs = ospace.nlocs();
        oops::Log::info() << "IodaStats:" << name << ": nlocs =" << nlocs << std::endl;
        if (nlocs == 0) {
            oops::Log::info() << "IodaStats:" << name << ": ObsSpace is empty, skipping " << obsFile << std::endl;
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
            // warn and continue to process other obs spaces
            oops::Log::info() << "IodaStats:" << name
                      << ": Cannot use channels with multiple variables, skipping." << std::endl;
            continue;
        }

        // get the lists of everything to process/compute
        std::vector<std::string> groups;
        std::vector<std::string> stats;
        std::vector<std::string> qcgroups;
        std::vector<eckit::LocalConfiguration> domains;
        
        obsSpace.get("groups to process", groups);
        obsSpace.get("qc groups", qcgroups);
        obsSpace.get("statistics to compute", stats);

        obsSpace.get("domains to process", domains);
        // loop over all domains and get their definitions
        std::vector<std::string> domainNames;
        std::vector<std::string> domainMaskVar1;
        std::vector<std::string> domainMaskVar2;
        std::vector<std::string> domainMaskVar3;
        std::vector<std::vector<float>> domainMaskVals1;
        std::vector<std::vector<float>> domainMaskVals2;
        std::vector<std::vector<float>> domainMaskVals3;
        for (int idom = 0; idom < domains.size(); idom++ ) {
            auto domain = domains[idom];
            eckit::LocalConfiguration domainConf(domain, "domain");
            std::string domainname;
            std::string maskvar1, maskvar2, maskvar3;
            std::vector<float> maskvals1, maskvals2, maskvals3;
            domainConf.get("name", domainname);
            if (domainConf.has("first mask variable")) {
                domainConf.get("first mask variable", maskvar1);
                domainConf.get("first mask range", maskvals1);
            }
            if (domainConf.has("second mask variable")) {
                domainConf.get("second mask variable", maskvar2);
                domainConf.get("second mask range", maskvals2);
            }
            if (domainConf.has("third mask variable")) {
                domainConf.get("third mask variable", maskvar3);
                domainConf.get("third mask range", maskvals3);
            }
            domainNames.push_back(domainname);
            domainMaskVar1.push_back(maskvar1);
            domainMaskVar2.push_back(maskvar2);
            domainMaskVar3.push_back(maskvar3);
            domainMaskVals1.push_back(maskvals1);
            domainMaskVals2.push_back(maskvals2);
            domainMaskVals3.push_back(maskvals3);
        }

        // Determine if we are doing regular binning for this obs space
        int nbins_x = 0;
        int nbins_y = 0;
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
        assert(qcgroups.size() == groups.size());

        // if the zBins are empty, create a dummy one for 'all'
        if (zBinNames.size() == 0){
            zBinNames.push_back("all");
            zBinMaskVar.push_back("latitude");
            zBinMaskVals.push_back(std::vector<float> {-90.0f, 90.0f});
        }

        // initialize netCDF output file for writing
        std::string outfile;
        obsSpace.get("output file", outfile);
        dautils::IodaStatsFile statfile;
        statfile.initializeNcfile(name, outfile, timeWindow,
                                  variables, channels, groups, stats,
                                  domainNames, nbins_x, nbins_y,
                                  zBinNames);
        
        // initialize ascii output file for writing if necessary stats are computed
        bool asciiOutput = false;
        std::string outasciifile = outfile;
        if (std::find(stats.begin(), stats.end(), "mean") != stats.end() &&
            std::find(stats.begin(), stats.end(), "count") != stats.end() &&
            std::find(stats.begin(), stats.end(), "rms") != stats.end() &&
            std::find(groups.begin(), groups.end(), "ombg") != groups.end() &&
            std::find(groups.begin(), groups.end(), "oman") != groups.end()) {
            // Place code here to initialize ascii output file
            size_t pos = outasciifile.rfind(".nc");
            if (pos != std::string::npos) {
                outasciifile.replace(pos, 3, ".txt");
            }
            statfile.initializeAsciiFile(name, outasciifile);
            asciiOutput = true;
        }

        // --------------------------------------------------------------------------
        // First, compute stats over specified domains (or global only)
        // --------------------------------------------------------------------------
        std::vector<std::vector<int>> mask(domains.size()+1, std::vector<int>(nlocs, 0));
        for (int idom = 0; idom < domains.size(); idom++ ) {
            // compute mask with function 3 times, one for each possible mask
            ObsStats obstatmask;
            std::vector<float> maskvalues(nlocs);
            if (!domainMaskVar1[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar1[idom], maskvalues);
                mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals1[idom][0], domainMaskVals1[idom][1], mask[idom]);
            }
            if (!domainMaskVar2[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar2[idom], maskvalues);
                mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals2[idom][0], domainMaskVals2[idom][1], mask[idom]);
            }
            if (!domainMaskVar3[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar3[idom], maskvalues);
                mask[idom] = obstatmask.update_mask(maskvalues, domainMaskVals3[idom][0], domainMaskVals3[idom][1], mask[idom]);
            }
        }

        // loop over variables
        for (int var = 0; var < variables.size(); var++) {
            // loop over groups
            for (int g = 0; g < groups.size(); g++) {
                std::vector<float> buffer(nlocs);
                std::vector<int> qcflag(nlocs);
                // we have to process differently if there are channels
                if (channels.empty()) {
                    // read the full variable
                    ospace.get_db(groups[g], variables[var], buffer);
                    // get the QC group
                    ospace.get_db(qcgroups[g], variables[var], qcflag);
                } else {
                    // give the list of channels to read
                    ospace.get_db(groups[g], variables[var], buffer, channels);
                    // get the QC group
                    ospace.get_db(qcgroups[g], variables[var], qcflag, channels);
                }
                // loop over domains
                for (int idom = 0; idom < domains.size()+1; idom++ ) {
                    // loop over stats
                    for (int s = 0; s < stats.size(); s++) {
                        
                    }  // loop over stats
                }  // loop over domains
            }  // loop over groups
        }  // loop over variables
    }  // loop over obs spaces
    return 0
}