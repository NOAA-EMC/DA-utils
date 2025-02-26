#pragma once

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

#include "./calcstats.h"
#include "./statfile.h"

namespace dautils {
  class IodaStats : public oops::Application {
    public:
      // -----------------------------------------------------------------------------
      explicit IodaStats(const eckit::mpi::Comm & comm = oops::mpi::world())
        : Application(comm), fillVal_(util::missingValue<float>()) {
        oceans_["Atlantic"] = 1;
        oceans_["Pacific"] = 2;
        oceans_["Indian"] = 3;
        oceans_["Arctic"] = 4;
        oceans_["Southern"] = 5;
        }
      static const std::string classname() {return "dautils::IodaStats";}
      // -----------------------------------------------------------------------------
      int execute(const eckit::Configuration & fullConfig) const {
        // define the time window
        const eckit::LocalConfiguration timeWindowConf(fullConfig, "time window");
        const util::TimeWindow timeWindow(timeWindowConf);

        // get the list of obs spaces to process
        std::vector<eckit::LocalConfiguration> obsSpaces;
        fullConfig.get("obs spaces", obsSpaces);

        // get the communicator for just me
        const eckit::mpi::Comm & mycomm = oops::mpi::myself();

        // get MPI comm size
        int nprocs = getComm().size();
        for (int i = getComm().rank(); i < obsSpaces.size(); i+= nprocs) {
          // get the configuration for this obs space
          auto obsSpace = obsSpaces[i];
          eckit::LocalConfiguration obsConfig(obsSpace, "obs space");

          // open the IODA file
          std::string obsFile;
          obsConfig.get("obsdatain.engine.obsfile", obsFile);
          oops::Log::info() << "IODA-Stats: Processing " << obsFile << std::endl;
          ioda::ObsSpace ospace(obsConfig, mycomm, timeWindow, mycomm);
          const size_t nlocs = ospace.nlocs();
          oops::Log::info() << obsFile << ": nlocs =" << nlocs << std::endl;

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

          // determine if we are doing regular binning for this obs space
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
          assert(groups.size() == qcgroups.size());

          // initialize netCDF output file for writing
          std::string outfile;
          obsSpace.get("output file", outfile);
          StatFile statfile;
          statfile.initializeNcfile(outfile, timeWindow, variables, channels, groups,
                                    stats, domainNames, nbins_x, nbins_y, zBinNames);

          // --------------------------------------------------------------------------
          // first, compute stats over specified domains (or global only)
          // --------------------------------------------------------------------------
          // loop over domains, compute the masks for each
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
              oops::Log::info() << obsFile << ": Now processing "
                                << groups[g] << "/" << variables[var] << std::endl;
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
                if (idom < domains.size()) {
                  oops::Log::info() << "Processing domain: " << domainNames[idom] << std::endl;
                  oops::Log::info() << domainMaskVar1[idom] << "=" <<  domainMaskVals1[idom] << ";"
                           << domainMaskVar2[idom] << "=" <<  domainMaskVals2[idom] << ";"
                           << domainMaskVar3[idom] << "=" <<  domainMaskVals3[idom] << ";"
                           << std::endl;
                }
                // loop over stats
                ObsStats obstat;
                for (int s = 0; s < stats.size(); s++) {
                  // Maybe eventually set this up as a factory but for now just do it
                  // with this old school if/else if way
                  std::vector<int> intstat;
                  std::vector<float> floatstat;
                  if (stats[s] == "count") {
                    intstat = obstat.getObsCount(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "Count:" << intstat << std::endl; 
                  } else if (stats[s] == "mean") {
                    floatstat = obstat.getMean(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "Mean:" << floatstat << std::endl; 
                  } else if (stats[s] == "RMS") {
                    floatstat = obstat.getRMS(buffer, qcflag, channels, mask[idom]);
                    oops::Log::info() << "RMS:" << floatstat << std::endl;
                  } else {
                    oops::Log::info() << stats[s] << " not supported. Skipping." << std::endl;
                  }
                  if (stats[s] == "count") {
                    statfile.writeByDomains(outfile, groups[g], variables[var],
                                            stats[s], idom, intstat);
                  } else {
                    statfile.writeByDomains(outfile, groups[g], variables[var],
                                            stats[s], idom, floatstat);
                  }
                }
              }
            }
          }
          // --------------------------------------------------------------------------
          // now, compute stats over binned regions, if applicable
          // --------------------------------------------------------------------------
          if (obsSpace.has("regular grid binning")) {
            // get lat/lon ranges based on bin sizes
            std::vector<float> longitudes(nbins_x+1);
            std::vector<float> latitudes(nbins_y+1);
            float dx = 360.0 / float(nbins_x);
            longitudes[0] = 0.0;
            latitudes[0] = -90.0;
            for (int ibin = 1; ibin < nbins_x+1; ibin++ ) {
              longitudes[ibin] = longitudes[ibin-1] + dx;
            }
            for (int ibin = 1; ibin < nbins_y+1; ibin++ ) {
              latitudes[ibin] = latitudes[ibin-1] + dx;
            }
            // loop over domains, compute the masks for each
            oops::Log::info() << "--------------------------------------------" << std::endl;
            int nzbins;
            if (channels.empty()) {
              nzbins = zBinNames.size();
            } else {
              nzbins = 1;
            }
            std::vector<std::vector<int>> binmask(nzbins * nbins_x * nbins_y, std::vector<int>(nlocs, 0));
            if (channels.empty()) { // use vertical bins
              int ibin = 0;
              for (int idom = 0; idom < zBinNames.size(); idom++ ) {
                oops::Log::info() << "Now processing binned data for vertical bin: " << zBinNames[idom] << std::endl;
                oops::Log::info() << "nbins_x: " << nbins_x << " nbins_y: " << nbins_y << std::endl;
                // compute masks for the bins
                ObsStats obstatbinmask;
                std::vector<float> xmaskvalues(nlocs), ymaskvalues(nlocs), zmaskvalues(nlocs);
                if (!zBinMaskVar[idom].empty()) {
                  ospace.get_db("MetaData", zBinMaskVar[idom], zmaskvalues);
                }
                ospace.get_db("MetaData", "latitude", ymaskvalues);
                ospace.get_db("MetaData", "longitude", xmaskvalues);
                for (int iy=0; iy < nbins_y; iy++) {
                  for (int ix=0; ix < nbins_x; ix++) {
                    ibin = ix + (iy * nbins_x) + (idom * nbins_x * nbins_y);
                    binmask[ibin] = obstatbinmask.update_mask(zmaskvalues, zBinMaskVals[idom][0], zBinMaskVals[idom][1], binmask[ibin]);
                    binmask[ibin] = obstatbinmask.update_mask(ymaskvalues, latitudes[iy], latitudes[iy+1], binmask[ibin]);
                    binmask[ibin] = obstatbinmask.update_mask(xmaskvalues, longitudes[ix], longitudes[ix+1], binmask[ibin]);
                  }
                }
              }
            } else { // assumes channels
              int ibin = 0;
              oops::Log::info() << "nbins_x: " << nbins_x << " nbins_y: " << nbins_y << std::endl;
              // compute masks for the bins
              ObsStats obstatbinmask;
              std::vector<float> xmaskvalues(nlocs), ymaskvalues(nlocs);
              ospace.get_db("MetaData", "latitude", ymaskvalues);
              ospace.get_db("MetaData", "longitude", xmaskvalues);
              for (int iy=0; iy < nbins_y; iy++) {
                for (int ix=0; ix < nbins_x; ix++) {
                  ibin = ix + (iy * nbins_x);
                  binmask[ibin] = obstatbinmask.update_mask(ymaskvalues, latitudes[iy], latitudes[iy+1], binmask[ibin]);
                  binmask[ibin] = obstatbinmask.update_mask(xmaskvalues, longitudes[ix], longitudes[ix+1], binmask[ibin]);
                }
              }
            }
            // loop over variables
            for (int var = 0; var < variables.size(); var++) {
              // loop over groups
              for (int g = 0; g < groups.size(); g++) {
                oops::Log::info() << obsFile << ": Now processing "
                                  << groups[g] << "/" << variables[var] << std::endl;
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
                ObsStats obstat;
                for (int s = 0; s < stats.size(); s++) {
                  if (channels.empty()) {
                    // loop over bins
                    int ibin = 0;
                    for (int idom = 0; idom < zBinNames.size(); idom++ ) {
                      std::vector<std::vector<float>> fullfloatstat(nbins_y, std::vector<float>(nbins_x,0.0));
                      std::vector<std::vector<int>> fullintstat(nbins_y, std::vector<int>(nbins_x,0.0));
                      for (int iy=0; iy < nbins_y; iy++) {
                        for (int ix=0; ix < nbins_x; ix++) {
                          ibin = ix + (iy * nbins_x) + (idom * nbins_x * nbins_y);
                          // Maybe eventually set this up as a factory but for now just do it
                          // with this old school if/else if way
                          std::vector<int> intstat;
                          std::vector<float> floatstat;
                          if (stats[s] == "count") {
                            intstat = obstat.getObsCount(buffer, qcflag, channels, binmask[ibin]);
                          } else if (stats[s] == "mean") {
                            floatstat = obstat.getMean(buffer, qcflag, channels, binmask[ibin]);
                          } else if (stats[s] == "RMS") {
                            floatstat = obstat.getRMS(buffer, qcflag, channels, binmask[ibin]);
                          }
                          if (stats[s] == "count") {
                            fullintstat[iy][ix] = intstat[0];
                          } else {
                            fullfloatstat[iy][ix] = floatstat[0];  
                          }
                        }
                      }
                      if (stats[s] == "count") {
                        statfile.writeByBins(outfile, groups[g], variables[var],
                                            stats[s], idom, nbins_y, nbins_x, fullintstat);                  
                      } else {
                        statfile.writeByBins(outfile, groups[g], variables[var],
                                            stats[s], idom, nbins_y, nbins_x, fullfloatstat);
                      }
                    }
                  } else { // variable has channels not vertical bins
                    std::vector<std::vector<std::vector<float>>> fullfloatstat(channels.size(),
                      std::vector<std::vector<float>>(nbins_y,std::vector<float>(nbins_x, 0.0)));
                    std::vector<std::vector<std::vector<int>>> fullintstat(channels.size(),
                      std::vector<std::vector<int>>(nbins_y,std::vector<int>(nbins_x, 0.0)));
                    int ibin = 0;
                    for (int iy=0; iy < nbins_y; iy++) {
                      for (int ix=0; ix < nbins_x; ix++) {
                        ibin = ix + (iy * nbins_x);
                        // Maybe eventually set this up as a factory but for now just do it
                        // with this old school if/else if way
                        std::vector<int> intstat;
                        std::vector<float> floatstat;
                        if (stats[s] == "count") {
                          intstat = obstat.getObsCount(buffer, qcflag, channels, binmask[ibin]);
                        } else if (stats[s] == "mean") {
                          floatstat = obstat.getMean(buffer, qcflag, channels, binmask[ibin]);
                        } else if (stats[s] == "RMS") {
                          floatstat = obstat.getRMS(buffer, qcflag, channels, binmask[ibin]);
                        }
                        // loop over channels
                        for (int ich = 0; ich < channels.size(); ich++ ) {
                          if (stats[s] == "count") {
                            fullintstat[ich][iy][ix] = intstat[ich];
                          } else {
                            fullfloatstat[ich][iy][ix] = floatstat[ich];  
                          }
                        }
                      }
                    }
                    // loop over channels
                    for (int ich = 0; ich < channels.size(); ich++ ) {
                      if (stats[s] == "count") {
                        statfile.writeByBins(outfile, groups[g], variables[var],
                                            stats[s], ich, nbins_y, nbins_x, fullintstat[ich]);                  
                      } else {
                        statfile.writeByBins(outfile, groups[g], variables[var],
                                            stats[s], ich, nbins_y, nbins_x, fullfloatstat[ich]);
                      }
                    }
                  } // channels vs bins
                } // loop over stats
              } // groups
            } // variables 
          } // if binning
        } // communicator / obs space loop
        return 0;
      }

    // -----------------------------------------------------------------------------
    // Data members
    std::map<std::string, int> oceans_;
    double fillVal_;
    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------
   private:
    std::string appname() const {
      return "dautils::IodaExample";
    }
    // -----------------------------------------------------------------------------
  };
}  // namespace dautils
