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
#include "calcstats.h"

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

        // Check that the QC groups list is the same size as groups
        if (groups.size() != qcgroups.size()) {
            throw eckit::Exception("QC groups list size does not match groups list size", Here());
        }

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
        std::string asciiZBinVar = "latitude"; // this is so that below when getting the total mask it works even if no vertical binning
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
        stattxtfile.initializeTxtFile(outasciifile, timeWindow, obsSpaceName, nlocs, obsSpace.has("channels"), asciiZBins);
        
        // first let us loop over the ascii vertical bins if they are defined, and a total if not
        std::vector<std::vector<int>> ascii_binmask(asciiZBins.size() + 1, std::vector<int>(nlocs, 0));
        std::vector<float> maskvalues(nlocs);
        ospace.get_db("MetaData", asciiZBinVar, maskvalues);
        // total column regardless of vertical binning
        ascii_binmask[0] = update_mask(maskvalues, -1.0e30f, 1.0e30f, ascii_binmask[0]);
        // make the first bin max value very large to include all values above min
        asciiZBinRanges.insert(asciiZBinRanges.begin(), 1.0e30f);
        for (int izbin = 0; izbin < asciiZBins.size(); izbin++) {
            ascii_binmask[izbin+1] = update_mask(maskvalues, asciiZBinRanges[izbin+1], asciiZBinRanges[izbin], ascii_binmask[izbin+1]);
        }

        // loop over variables
        for (int var = 0; var < variables.size(); var++) {
            // loop over groups
            for (int g = 0; g < groups.size(); g++) {
                std::vector<float> buffer(nlocs);
                std::vector<int> qcflag(nlocs);
                std::vector<float> errorvals(nlocs);
                // we have to process differently if there are channels
                if (channels.empty()) {
                    // read the full variable
                    ospace.get_db(groups[g], variables[var], buffer);
                    // get the QC group
                    ospace.get_db(qcgroups[g], variables[var], qcflag);
                    // get the error group if provided
                    if (!errorGroups[g].empty()) {
                        ospace.get_db(errorGroups[g], variables[var], errorvals);
                    }
                } else {
                    // give the list of channels to read
                    ospace.get_db(groups[g], variables[var], buffer, channels);
                    // get the QC group
                    ospace.get_db(qcgroups[g], variables[var], qcflag, channels);
                    // get the error group if provided
                    if (!errorGroups[g].empty()) {
                        ospace.get_db(errorGroups[g], variables[var], errorvals, channels);
                    }
                }
                // loop over stats
                for (int s = 0; s < stats.size(); s++) {
                    std::vector<std::vector<std::vector<int>>> intstat;
                    std::vector<std::vector<std::vector<float>>> floatstat;
                    // loop over vertical bins
                    for (int izbin = 0; izbin < asciiZBins.size()+1; izbin++) {
                        if (stats[s] == "count") {
                            intstat.push_back(getObsCount(buffer, qcflag, channels, ascii_binmask[izbin]));
                        } else if (stats[s] == "mean") {
                            floatstat.push_back(getMean(buffer, qcflag, channels, ascii_binmask[izbin]));
                        }
                        else if (stats[s] == "RMS") {
                            floatstat.push_back(getRMS(buffer, qcflag, channels, ascii_binmask[izbin]));
                        } else {
                            oops::Log::info() << stats[s] << " not supported. Skipping." << std::endl;
                        }
                    }
                    // reshape the vectors for writing to the ASCII file
                    int nch = 1;
                    if (!channels.empty()) nch = channels.size();
                    std::vector<std::vector<std::vector<int>>> intstat_reshaped(3,
                        std::vector<std::vector<int>>(asciiZBins.size()+1,
                            std::vector<int>(nch, 0)));
                    std::vector<std::vector<std::vector<float>>> floatstat_reshaped(3,
                        std::vector<std::vector<float>>(asciiZBins.size()+1,
                            std::vector<float>(nch, 0.0)));
                    for (int ibin = 0; ibin < asciiZBins.size()+1; ibin++) {
                        for (int ich = 0; ich < nch; ich++) {
                            if (stats[s] == "count") {
                                intstat_reshaped[0][ibin][ich] = intstat[ibin][ich][0];
                                intstat_reshaped[1][ibin][ich] = intstat[ibin][ich][1];
                                intstat_reshaped[2][ibin][ich] = intstat[ibin][ich][2];
                            } else {
                                floatstat_reshaped[0][ibin][ich] = floatstat[ibin][ich][0];
                                floatstat_reshaped[1][ibin][ich] = floatstat[ibin][ich][1];
                                floatstat_reshaped[2][ibin][ich] = floatstat[ibin][ich][2];
                            }
                        }
                    }
                    // write to ASCII file
                    if (channels.empty()) {
                        std::vector<int> ch = {-1};
                        if (stats[s] == "count") {
                            stattxtfile.writeTxtStat(obsSpaceName, variables[var], ch, groups[g],
                                                     stats[s], intstat_reshaped);
                        } else {
                            stattxtfile.writeTxtStat(obsSpaceName, variables[var], ch, groups[g],
                                                     stats[s], floatstat_reshaped);
                        }
                    } else {
                        if (stats[s] == "count") {
                            stattxtfile.writeTxtStat(obsSpaceName, variables[var], channels, groups[g],
                                                    stats[s], intstat_reshaped);
                        } else {
                            stattxtfile.writeTxtStat(obsSpaceName, variables[var], channels, groups[g],
                                                    stats[s], floatstat_reshaped);
                        }
                    }
                } // end of stats loop
            } // end of group loop
        } // end of variable loop
        // close the ascii file
        stattxtfile.closeFile();
        oops::Log::info() << "Finished writing ASCII stats for " << obsSpaceName << std::endl;
        // Now let's process the netCDF files by domains and/or bins
        // --------------------------------------------------------------------------
        // first, compute stats over specified domains (or global only)
        // --------------------------------------------------------------------------
        // loop over domains, compute the masks for each
        std::vector<std::vector<int>> mask(domains.size()+1, std::vector<int>(nlocs, 0));
        for (int idom = 0; idom < domains.size(); idom++ ) {
            // compute mask with function 3 times, one for each possible mask
            std::vector<float> maskvalues(nlocs);
            if (!domainMaskVar1[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar1[idom], maskvalues);
                // Convert longitudes if this is a longitude mask
                if (domainMaskVar1[idom] == "longitude") {
                    convertLongitudes(maskvalues);
                }
                mask[idom] = update_mask(maskvalues, domainMaskVals1[idom][0], domainMaskVals1[idom][1], mask[idom]);
            }
            if (!domainMaskVar2[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar2[idom], maskvalues);
                // Convert longitudes if this is a longitude mask
                if (domainMaskVar2[idom] == "longitude") {
                    convertLongitudes(maskvalues);
                }
                mask[idom] = update_mask(maskvalues, domainMaskVals2[idom][0], domainMaskVals2[idom][1], mask[idom]);
            }
            if (!domainMaskVar3[idom].empty()) {
                ospace.get_db("MetaData", domainMaskVar3[idom], maskvalues);
                // Convert longitudes if this is a longitude mask
                if (domainMaskVar3[idom] == "longitude") {
                    convertLongitudes(maskvalues);
                }
                mask[idom] = update_mask(maskvalues, domainMaskVals3[idom][0], domainMaskVals3[idom][1], mask[idom]);
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
                        // Maybe eventually set this up as a factory but for now just do it
                        // with this old school if/else if way
                        std::vector<std::vector<int>> intstat;
                        std::vector<std::vector<float>> floatstat;
                        if (stats[s] == "count") {
                            intstat = getObsCount(buffer, qcflag, channels, mask[idom]);
                        } else if (stats[s] == "mean") {
                            floatstat = getMean(buffer, qcflag, channels, mask[idom]);
                        } else if (stats[s] == "RMS") {
                            floatstat = getRMS(buffer, qcflag, channels, mask[idom]);
                        } else {
                            oops::Log::info() << stats[s] << " not supported. Skipping." << std::endl;
                        }
                        int nch = 1;
                        if (!channels.empty()) nch = channels.size();
                        if (stats[s] == "count") {
                            std::vector<int> intstat_assim(nch), intstat_monit(nch), intstat_rej(nch);
                            for (size_t ich = 0; ich < nch; ich++) {
                                intstat_assim[ich] = intstat[ich][0];
                                intstat_monit[ich] = intstat[ich][1];
                                intstat_rej[ich] = intstat[ich][2];
                            }
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "assimilated_" + stats[s], idom, intstat_assim);
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "monitored_" + stats[s], idom, intstat_monit);
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "rejected_" + stats[s], idom, intstat_rej);
                        } else {
                            std::vector<float> floatstat_assim(nch), floatstat_monit(nch), floatstat_rej(nch);
                            for (size_t ich = 0; ich < nch; ich++) {
                                floatstat_assim[ich] = floatstat[ich][0];
                                floatstat_monit[ich] = floatstat[ich][1];
                                floatstat_rej[ich] = floatstat[ich][2];
                            }
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "assimilated_" + stats[s], idom, floatstat_assim);
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "monitored_" + stats[s], idom, floatstat_monit);
                            statncfile.writeByDomains(outncfile, groups[g], variables[var],
                                                    "rejected_" + stats[s], idom, floatstat_rej);
                        }
                    } // end of stats loop
                } // end of domain loop
            } // end of group loop
        } // end of variable loop
        // --------------------------------------------------------------------------
        // now, compute stats over binned regions, if applicable
        // --------------------------------------------------------------------------
        if (obsSpace.has("regular grid binning")) {
            // Read longitude data to determine if conversion is needed
            std::vector<float> lon_sample(nlocs);
            ospace.get_db("MetaData", "longitude", lon_sample);
            
            // Check if longitudes need conversion from 0-360 to -180 to 180
            float minLon = lon_sample[0];
            float maxLon = lon_sample[0];
            for (const auto& lon : lon_sample) {
                if (lon < minLon) minLon = lon;
                if (lon > maxLon) maxLon = lon;
            }
            bool needsConversion = (minLon >= 0.0 && maxLon > 180.0);
            
            // figure out if we are 0-360 or -180-180 longitudes
            eckit::LocalConfiguration binConfig;
            obsSpace.get("regular grid binning", binConfig);
            bool negLon = needsConversion;  // Use detected range instead of config
            if (binConfig.has("use negative longitudes")) {
                binConfig.get("use negative longitudes", negLon);
                // Override config if data needs conversion
                if (needsConversion) {
                    negLon = true;
                    oops::Log::info() << "Data has longitudes in 0-360 range, will convert to -180 to 180" << std::endl;
                }
            } else if (needsConversion) {
                oops::Log::info() << "Data has longitudes in 0-360 range, will convert to -180 to 180" << std::endl;
            }
            // get lat/lon ranges based on bin sizes
            std::vector<float> longitudes(nbins_x+1);
            std::vector<float> latitudes(nbins_y+1);
            float dx = 360.0 / float(nbins_x);
            if (negLon) {
                longitudes[0] = -180.0;
            } else {
                longitudes[0] = 0.0;
            }
            latitudes[0] = -90.0;
            for (int ibin = 1; ibin < nbins_x+1; ibin++ ) {
                longitudes[ibin] = longitudes[ibin-1] + dx;
            }
            for (int ibin = 1; ibin < nbins_y+1; ibin++ ) {
                latitudes[ibin] = latitudes[ibin-1] + dx;
            }

            // loop over domains, compute the masks for each
            int nzbins;
            if (channels.empty()) {
                nzbins = zBinNames.size();
            } else {
                nzbins = 1;
            }
            std::vector<std::vector<int>> binmask(nzbins * nbins_x * nbins_y, std::vector<int>(nlocs, 0));
            if (channels.empty()) { // use vertical bins
                int ibin = 0;
                for (int idom = 0; idom < nzbins; idom++ ) {
                    oops::Log::info() << "Now processing binned data for vertical bin: " << zBinNames[idom] << std::endl;
                    oops::Log::info() << "nbins_x: " << nbins_x << " nbins_y: " << nbins_y << std::endl;
                    // compute masks for the bins
                    std::vector<float> xmaskvalues(nlocs), ymaskvalues(nlocs), zmaskvalues(nlocs);
                    if (!zBinMaskVar[idom].empty()) {
                        ospace.get_db("MetaData", zBinMaskVar[idom], zmaskvalues);
                    }
                    ospace.get_db("MetaData", "latitude", ymaskvalues);
                    ospace.get_db("MetaData", "longitude", xmaskvalues);
                    // Convert longitudes from 0-360 to -180 to 180 if needed
                    convertLongitudes(xmaskvalues);
                    for (int iy=0; iy < nbins_y; iy++) {
                        for (int ix=0; ix < nbins_x; ix++) {
                            ibin = ix + (iy * nbins_x) + (idom * nbins_x * nbins_y);
                            binmask[ibin] = update_mask(zmaskvalues, zBinMaskVals[idom][0], zBinMaskVals[idom][1], binmask[ibin]);
                            binmask[ibin] = update_mask(ymaskvalues, latitudes[iy], latitudes[iy+1], binmask[ibin]);
                            binmask[ibin] = update_mask(xmaskvalues, longitudes[ix], longitudes[ix+1], binmask[ibin]);
                        }
                    }
                }
            } else { // assumes channels
                int ibin = 0;
                oops::Log::info() << "nbins_x: " << nbins_x << " nbins_y: " << nbins_y << std::endl;
                // compute masks for the bins
                std::vector<float> xmaskvalues(nlocs), ymaskvalues(nlocs);
                ospace.get_db("MetaData", "latitude", ymaskvalues);
                ospace.get_db("MetaData", "longitude", xmaskvalues);
                // Convert longitudes from 0-360 to -180 to 180 if needed
                convertLongitudes(xmaskvalues);
                for (int iy=0; iy < nbins_y; iy++) {
                    for (int ix=0; ix < nbins_x; ix++) {
                        ibin = ix + (iy * nbins_x);
                        binmask[ibin] = update_mask(ymaskvalues, latitudes[iy], latitudes[iy+1], binmask[ibin]);
                        binmask[ibin] = update_mask(xmaskvalues, longitudes[ix], longitudes[ix+1], binmask[ibin]);
                    }
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
                    // loop over stats
                    for (int s = 0; s < stats.size(); s++) {
                        if (channels.empty()) {
                            // loop over bins
                            int ibin = 0;
                            for (int idom = 0; idom < nzbins; idom++ ) {
                                std::vector<std::vector<float>> fullfloatstat_assim(nbins_y, std::vector<float>(nbins_x,0.0));
                                std::vector<std::vector<int>> fullintstat_assim(nbins_y, std::vector<int>(nbins_x,0.0));
                                std::vector<std::vector<float>> fullfloatstat_monitored(nbins_y, std::vector<float>(nbins_x,0.0));
                                std::vector<std::vector<int>> fullintstat_monitored(nbins_y, std::vector<int>(nbins_x,0.0));
                                std::vector<std::vector<float>> fullfloatstat_rejected(nbins_y, std::vector<float>(nbins_x,0.0));
                                std::vector<std::vector<int>> fullintstat_rejected(nbins_y, std::vector<int>(nbins_x,0.0));
                                
                                for (int iy=0; iy < nbins_y; iy++) {
                                    for (int ix=0; ix < nbins_x; ix++) {
                                        ibin = ix + (iy * nbins_x) + (idom * nbins_x * nbins_y);
                                        // Maybe eventually set this up as a factory but for now just do it
                                        // with this old school if/else if way
                                        std::vector<std::vector<int>> intstat;
                                        std::vector<std::vector<float>> floatstat;
                                        if (stats[s] == "count") {
                                            intstat = getObsCount(buffer, qcflag, channels, binmask[ibin]);
                                        } else if (stats[s] == "mean") {
                                            floatstat = getMean(buffer, qcflag, channels, binmask[ibin]);
                                        } else if (stats[s] == "RMS") {
                                            floatstat = getRMS(buffer, qcflag, channels, binmask[ibin]);
                                        }

                                        if (stats[s] == "count") {
                                            std::vector<int> intstat_assim, intstat_monit, intstat_rej;
                                            for (auto val : intstat) {
                                                intstat_assim.push_back(val[0]);
                                                intstat_monit.push_back(val[1]);
                                                intstat_rej.push_back(val[2]);
                                            }
                                            fullintstat_assim[iy][ix] = intstat_assim[0];
                                            fullintstat_monitored[iy][ix] = intstat_monit[0];
                                            fullintstat_rejected[iy][ix] = intstat_rej[0];
                                        } else {
                                            std::vector<float> floatstat_assim, floatstat_monit, floatstat_rej;
                                            for (auto val : floatstat) {
                                                floatstat_assim.push_back(val[0]);
                                                floatstat_monit.push_back(val[1]);
                                                floatstat_rej.push_back(val[2]);
                                            }
                                            fullfloatstat_assim[iy][ix] = floatstat_assim[0];
                                            fullfloatstat_monitored[iy][ix] = floatstat_monit[0];
                                            fullfloatstat_rejected[iy][ix] = floatstat_rej[0];
                                        }
                                   } // end of ix
                                } // end of iy
                                if (stats[s] == "count") {
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "assimilated_" + stats[s], idom, nbins_y, nbins_x, fullintstat_assim);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "monitored_" + stats[s], idom, nbins_y, nbins_x, fullintstat_monitored);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "rejected_" + stats[s], idom, nbins_y, nbins_x, fullintstat_rejected);
                                } else {
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "assimilated_" + stats[s], idom, nbins_y, nbins_x, fullfloatstat_assim);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "monitored_" + stats[s], idom, nbins_y, nbins_x, fullfloatstat_monitored);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "rejected_" + stats[s], idom, nbins_y, nbins_x, fullfloatstat_rejected);
                                }
                            }
                        } else { // variable has channels not vertical bins
                            std::vector<std::vector<std::vector<float>>> fullfloatstat_assim(channels.size(),
                              std::vector<std::vector<float>>(nbins_y,std::vector<float>(nbins_x, 0.0)));
                            std::vector<std::vector<std::vector<int>>> fullintstat_assim(channels.size(),
                              std::vector<std::vector<int>>(nbins_y,std::vector<int>(nbins_x, 0.0)));
                            std::vector<std::vector<std::vector<float>>> fullfloatstat_monitored(channels.size(),
                              std::vector<std::vector<float>>(nbins_y,std::vector<float>(nbins_x, 0.0)));
                            std::vector<std::vector<std::vector<int>>> fullintstat_monitored(channels.size(),
                              std::vector<std::vector<int>>(nbins_y,std::vector<int>(nbins_x, 0.0)));
                            std::vector<std::vector<std::vector<float>>> fullfloatstat_rejected(channels.size(),
                              std::vector<std::vector<float>>(nbins_y,std::vector<float>(nbins_x, 0.0)));
                            std::vector<std::vector<std::vector<int>>> fullintstat_rejected(channels.size(),
                              std::vector<std::vector<int>>(nbins_y,std::vector<int>(nbins_x, 0.0)));
                            int nch = 1;
                            if (!channels.empty()) nch = channels.size();
                            int ibin = 0;
                            for (int iy=0; iy < nbins_y; iy++) {
                                for (int ix=0; ix < nbins_x; ix++) {
                                    ibin = ix + (iy * nbins_x);
                                    // Maybe eventually set this up as a factory but for now just do it
                                    // with this old school if/else if way
                                    std::vector<std::vector<int>> intstat;
                                    std::vector<std::vector<float>> floatstat;
                                    if (stats[s] == "count") {
                                        intstat = getObsCount(buffer, qcflag, channels, binmask[ibin]);
                                    } else if (stats[s] == "mean") {
                                        floatstat = getMean(buffer, qcflag, channels, binmask[ibin]);
                                    } else if (stats[s] == "RMS") {
                                        floatstat = getRMS(buffer, qcflag, channels, binmask[ibin]);
                                    }
                                    // loop over channels
                                    for (int ich = 0; ich < nch; ich++ ) {
                                        if (stats[s] == "count") {
                                            fullintstat_assim[ich][iy][ix] = intstat[ich][0];
                                            fullintstat_monitored[ich][iy][ix] = intstat[ich][1];
                                            fullintstat_rejected[ich][iy][ix] = intstat[ich][2];
                                        } else {
                                            fullfloatstat_assim[ich][iy][ix] = floatstat[ich][0];
                                            fullfloatstat_monitored[ich][iy][ix] = floatstat[ich][1];
                                            fullfloatstat_rejected[ich][iy][ix] = floatstat[ich][2];
                                        }
                                    }
                                }
                            }
                            // loop over channels
                            for (int ich = 0; ich < channels.size(); ich++ ) {
                                if (stats[s] == "count") {
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "assimilated_" + stats[s], ich, nbins_y, nbins_x, fullintstat_assim[ich]);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "monitored_" + stats[s], ich, nbins_y, nbins_x, fullintstat_monitored[ich]);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "rejected_" + stats[s], ich, nbins_y, nbins_x, fullintstat_rejected[ich]);
                                } else {
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "assimilated_" + stats[s], ich, nbins_y, nbins_x, fullfloatstat_assim[ich]);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "monitored_" + stats[s], ich, nbins_y, nbins_x, fullfloatstat_monitored[ich]);
                                    statncfile.writeByBins(outncfile, groups[g], variables[var],
                                                        "rejected_" + stats[s], ich, nbins_y, nbins_x, fullfloatstat_rejected[ich]);
                                }
                            } // end of channels loop
                        } // channels or bins
                    } // end of stats loop
                } // end of group loop
            } // end of variable loop
        } // end of regular binning check
    }// end of obs space loop
} // end of run

// -----------------------------------------------------------------------------
// Convert longitudes from 0-360 to -180 to 180 range if needed
void dautils::convertLongitudes(std::vector<float>& longitudes) {
    if (longitudes.empty()) return;

    // Find min and max longitude values
    float minLon = longitudes[0];
    float maxLon = longitudes[0];
    for (const auto& lon : longitudes) {
        if (lon < minLon) minLon = lon;
        if (lon > maxLon) maxLon = lon;
    }

    // Check if longitudes are in 0-360 range
    // If min is >= 0 and max is > 180, we assume 0-360 range
    if (minLon >= 0.0 && maxLon > 180.0) {
        oops::Log::debug() << "Converting longitudes from 0-360 to -180 to 180 range" << std::endl;
        oops::Log::trace() << "Original range: [" << minLon << ", " << maxLon << "]" << std::endl;
        
        // Convert: if lon > 180, subtract 360
        for (auto& lon : longitudes) {
            if (lon > 180.0) {
                lon -= 360.0;
            }
        }
        
        // Find new min and max for logging
        minLon = longitudes[0];
        maxLon = longitudes[0];
        for (const auto& lon : longitudes) {
            if (lon < minLon) minLon = lon;
            if (lon > maxLon) maxLon = lon;
            }
        oops::Log::trace() << "Converted range: [" << minLon << ", " << maxLon << "]" << std::endl;
    }
}
