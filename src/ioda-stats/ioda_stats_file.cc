int dautils::IodaStatsFile::initializeNcfile(
    const std::string obspace, const std::string filename, const util::TimeWindow & timeWindow,
    std::vector<std::string> variables, std::vector<int> channels,
    std::vector<std::string> groups, std::vector<std::string> stats,
    std::vector<std::string> domainNames,
    int nbins_x, int nbins_y,
    std::vector<std::string> bins_z) {
    
    // Initialize the netCDF file for writing obs space statistics
    netCDF::NcFile ncFile(filename, netCDF::NcFile::replace);
    oops::Log::info() << "IodaStats:" << obspace << ": Opening "
                      << filename << " for writing..." << std::endl;
    // create an unlimited time dimension
    netCDF::NcDim tDim = ncFile.addDim("analysisCycle");
    // create domain dimension
    int ndomains = domainNames.size() + 1;
    netCDF::NcDim dDim = ncFile.addDim("Domain", ndomains);
    // vector of dimensions
    std::vector<netCDF::NcDim> domainDimVector, binningDimVector;
    domainDimVector.push_back(tDim);
    domainDimVector.push_back(dDim);
    // if channel is not empty, create a channel dimension
    netCDF::NcDim cDim;
    if (!channels.empty()) {
        cDim = ncFile.addDim("Channel", channels.size());
        domainDimVector.push_back(cDim);
    }
    // if nbins_x or nbins_y are nonzero, make them dimensions too
    netCDF::NcDim xDim, yDim, zDim;
    if (nbins_x > 0 || nbins_y > 0) {
        binningDimVector.push_back(tDim);
        if (!channels.empty()) {
            binningDimVector.push_back(cDim);
        } else {
            if (bins_z.size() > 0) {
                zDim = ncFile.addDim("binsZDim", bins_z.size());
                binningDimVector.push_back(zDim);
            }
        }
    }

    if (nbins_y > 0) {
        yDim = ncFile.addDim("binsYDim", nbins_y);
        binningDimVector.push_back(yDim);
    }
    if (nbins_x > 0) {
        xDim = ncFile.addDim("binsXDim", nbins_x);
        binningDimVector.push_back(xDim);
    }

    // create validTime variable
    netCDF::NcVar time = ncFile.addVar("validTime", netCDF::ncString, tDim);
    // put the analysis time in the file
    util::DateTime analysisTime = timeWindow.midpoint();
    std::vector<size_t> idxout;
    idxout.push_back(0);
    time.putVar(idxout, analysisTime.toString());

    // create domain variable
    netCDF::NcVar domain = ncFile.addVar("statisticDomain", netCDF::ncString, dDim);
    domainNames.push_back("Global");
    for (int idom = 0; idom < ndomains; idom++) {
        std::vector<size_t> idxdom;
        idxdom.push_back(idom);
        domain.putVar(idxdom, domainNames[idom]);
    }

    // create vertical bin variable
    if (channels.empty() && bins_z.size() > 0 && (nbins_x > 0 || nbins_y > 0)) {
        netCDF::NcVar zbins = ncFile.addVar("verticalBin", netCDF::ncString, zDim);
        for (int ibin = 0; ibin < bins_z.size(); ibin++) {
            std::vector<size_t> idxbin;
            idxbin.push_back(ibin);
            zbins.putVar(idxbin, bins_z[ibin]);
        }
    }

    // loop over group, then variables, then stats to create byDomains/group/var/stat in file
    netCDF::NcGroup domaingroup = ncFile.addGroup("byDomains");
    for (int g = 0; g < groups.size(); g++) {
        // create group group
        netCDF::NcGroup group = domaingroup.addGroup(groups[g]);
        // loop over variables
        for (int var = 0; var < variables.size(); var++) {
            // create variable group
            netCDF::NcGroup group2 = group.addGroup(variables[var]);
            // loop over statistics to write out
            for (int s = 0; s < stats.size(); s++) {
                netCDF::NcVar varout;
                if (stats[s] == "count") {
                    varout = group2.addVar(stats[s], netCDF::ncInt, domainDimVector);
                } else {
                    varout = group2.addVar(stats[s], netCDF::ncFloat, domainDimVector);
                    varout.putAtt("_FillValue", netCDF::ncFloat, fillVal_);
                }
            }
        }
    }

    // loop over group, then variables, then stats to create griddedBins/group/var/stat in file
    if (nbins_x > 0 || nbins_y > 0) {
        netCDF::NcGroup bingroup = ncFile.addGroup("griddedBins");
        for (int g = 0; g < groups.size(); g++) {
            // create group group
            netCDF::NcGroup group = bingroup.addGroup(groups[g]);
            // loop over variables
            for (int var = 0; var < variables.size(); var++) {
                // create variable group
                netCDF::NcGroup group2 = group.addGroup(variables[var]);
                // loop over statistics to write out
                for (int s = 0; s < stats.size(); s++) {
                    netCDF::NcVar varout;
                    if (stats[s] == "count") {
                        varout = group2.addVar(stats[s], netCDF::ncInt, binningDimVector);
                    } else {
                        varout = group2.addVar(stats[s], netCDF::ncFloat, binningDimVector);
                        varout.putAtt("_FillValue", netCDF::ncFloat, fillVal_);
                    }
                }
            }
        }
    }

    oops::Log::info() << "IodaStats:" << obspace << ": Output file:  "
                      << filename << " has been created." << std::endl;
    return 0;
}

int dautils::IodaStatsFile::initializeAsciifile(const std::string obspace,
                                                const std::string filename) {
    // Initialize the ASCII file for writing obs space statistics
    oops::Log::info() << "IodaStats:" << obspace << ": Opening "
                      << filename << " for writing..." << std::endl;

    // Open the file for writing
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        oops::Log::error() << "IodaStats:" << obspace << ": Cannot open file " << filename << " for writing."
                           << std::endl;
        return -1;
    }

    // Write header information
    outfile << "IODA Obs Space   Variable Name               "
            << "OmB Counts   OmB Mean    OmB RMS OmA Counts   OmA Mean    OmA RMS\n";
    outfile << "---------------------------------------------"
    outfile << "-----------------------------------------------------------------\n";

    outfile.close();
    oops::Log::info() << "IodaStats:" << obspace << ": Output file:  "
                      << filename << " has been created." << std::endl;
    return 0;
}
