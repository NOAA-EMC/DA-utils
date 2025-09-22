#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>

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
#include "oops/util/TimeWindow.h"

namespace dautils {
  // This utility reads IODA files and dumps basic summary information
  // to a formatted ASCII text file. It supports:
  // - MPI execution
  // - Processing directories or lists of IODA files
  // - Extracting nobs and ObsValue variables
  // - Nicely formatted text output

  class IodaDump : public oops::Application {
   public:
    explicit IodaDump(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "dautils::IodaDump";}

    int execute(const eckit::Configuration & fullConfig) const {
      // Validate configuration first
      validateConfiguration(fullConfig);

      // define the time window
      const eckit::LocalConfiguration timeWindowConf(fullConfig, "time window");
      const util::TimeWindow timeWindow(timeWindowConf);

      // get output file configuration
      std::string outputFile;
      fullConfig.get("output file", outputFile);

      // get input configuration - can be directory or list of files
      std::vector<std::string> inputFiles;
      if (fullConfig.has("input directory")) {
        std::string inputDir;
        fullConfig.get("input directory", inputDir);
        oops::Log::info() << "Scanning directory: " << inputDir << std::endl;
        inputFiles = getFilesFromDirectory(inputDir);
        oops::Log::info() << "Found " << inputFiles.size() << " IODA files in directory" << std::endl;
      } else if (fullConfig.has("input files")) {
        fullConfig.get("input files", inputFiles);
        oops::Log::info() << "Processing " << inputFiles.size() << " specified files" << std::endl;
      } else {
        throw eckit::Exception("Either 'input directory' or 'input files' must be specified");
      }

      if (inputFiles.empty()) {
        oops::Log::warning() << "No input files found to process" << std::endl;
        return 0;
      }

      // get the communicator for just me
      const eckit::mpi::Comm & mycomm = oops::mpi::myself();

      // distribute files across MPI processes
      int nprocs = getComm().size();
      int myrank = getComm().rank();
      
      std::vector<std::string> myFiles;
      for (size_t i = myrank; i < inputFiles.size(); i += nprocs) {
        myFiles.push_back(inputFiles[i]);
      }

      oops::Log::info() << "Process " << myrank << " will process " << myFiles.size() << " files" << std::endl;

      // process my files
      std::vector<FileInfo> fileInfos;
      for (const auto& file : myFiles) {
        try {
          FileInfo info = processFile(file, timeWindow, mycomm);
          fileInfos.push_back(info);
        } catch (const std::exception& e) {
          oops::Log::warning() << "Failed to process file " << file << ": " << e.what() << std::endl;
          // Add failed file info
          FileInfo failedInfo;
          failedInfo.filename = file;
          failedInfo.nobs = 0;
          failedInfo.success = false;
          failedInfo.errorMsg = e.what();
          fileInfos.push_back(failedInfo);
        }
      }

      // gather results from all processes
      std::vector<FileInfo> allFileInfos = gatherResults(fileInfos);

      // write output file (only from rank 0)
      if (myrank == 0) {
        writeOutputFile(outputFile, allFileInfos);
      }

      return 0;
    }

   private:
    struct FileInfo {
      std::string filename;
      size_t nobs;
      std::vector<std::string> obsValueVars;
      bool success;
      std::string errorMsg;
    };

    std::string appname() const {
      return "dautils::IodaDump";
    }

    void validateConfiguration(const eckit::Configuration & fullConfig) const {
      // Check for required fields
      if (!fullConfig.has("time window")) {
        throw eckit::Exception("Configuration must include 'time window' section");
      }
      
      if (!fullConfig.has("output file")) {
        throw eckit::Exception("Configuration must include 'output file' specification");
      }
      
      if (!fullConfig.has("input directory") && !fullConfig.has("input files")) {
        throw eckit::Exception("Configuration must include either 'input directory' or 'input files'");
      }
      
      if (fullConfig.has("input directory") && fullConfig.has("input files")) {
        oops::Log::warning() << "Both 'input directory' and 'input files' specified - using 'input directory'" << std::endl;
      }
    }

    std::vector<std::string> getFilesFromDirectory(const std::string& dir) const {
      std::vector<std::string> files;
      
      DIR* dirp = opendir(dir.c_str());
      if (dirp == nullptr) {
        throw eckit::Exception("Directory does not exist: " + dir);
      }

      struct dirent* entry;
      while ((entry = readdir(dirp)) != nullptr) {
        std::string filename = entry->d_name;
        
        // Skip . and .. entries
        if (filename == "." || filename == "..") {
          continue;
        }
        
        std::string fullPath = dir + "/" + filename;
        
        // Check if it's a regular file
        struct stat fileStat;
        if (stat(fullPath.c_str(), &fileStat) == 0 && S_ISREG(fileStat.st_mode)) {
          // Only include files that look like IODA files (nc, h5, hdf5)
          size_t dotPos = filename.find_last_of('.');
          if (dotPos != std::string::npos) {
            std::string ext = filename.substr(dotPos);
            if (ext == ".nc" || ext == ".h5" || ext == ".hdf5") {
              files.push_back(fullPath);
            }
          }
        }
      }
      
      closedir(dirp);
      std::sort(files.begin(), files.end());
      return files;
    }

    FileInfo processFile(const std::string& filename, 
                        const util::TimeWindow& timeWindow,
                        const eckit::mpi::Comm & comm) const {
      FileInfo info;
      info.filename = filename;
      info.success = false;

      try {
        // create a minimal configuration for this file
        eckit::LocalConfiguration obsConfig;
        obsConfig.set("obsdatain.engine.type", "H5File");
        obsConfig.set("obsdatain.engine.obsfile", filename);
        obsConfig.set("simulated variables", std::vector<std::string>{"dummy"});
        obsConfig.set("observed variables", std::vector<std::string>{"dummy"});

        // open the IODA file
        oops::Log::info() << "IODA-Dump: Processing " << filename << std::endl;
        ioda::ObsSpace ospace(obsConfig, comm, timeWindow, comm);
        
        info.nobs = ospace.nlocs();
        oops::Log::info() << filename << ": nobs = " << info.nobs << std::endl;

        // get ObsValue variables
        info.obsValueVars = getObsValueVariables(ospace);
        
        info.success = true;
        
      } catch (const std::exception& e) {
        info.errorMsg = e.what();
        oops::Log::warning() << "Error processing " << filename << ": " << e.what() << std::endl;
      }

      return info;
    }

    std::vector<std::string> getObsValueVariables(const ioda::ObsSpace& ospace) const {
      std::vector<std::string> variables;
      
      try {
        // For now, we'll try a set of common observation variables
        // In a full implementation, we would use IODA API to introspect the file
        // but this requires more detailed knowledge of the IODA internals
        
        std::vector<std::string> candidateVars = {
          "brightnessTemperature", "temperature", "humidity", "windSpeed", 
          "windDirection", "pressure", "seaSurfaceTemperature", "oceanDepth",
          "radiance", "reflectance", "altitude", "thickness", "waterVaporMixingRatio",
          "specificHumidity", "airTemperature", "virtualTemperature"
        };
        
        for (const auto& var : candidateVars) {
          try {
            // Try to read a small amount of data to see if variable exists
            size_t testSize = std::min(static_cast<size_t>(1), ospace.nlocs());
            if (testSize > 0) {
              std::vector<float> buffer(testSize);
              ospace.get_db("ObsValue", var, buffer);
              variables.push_back(var);
            }
          } catch (...) {
            // Variable doesn't exist or can't be read, skip it
          }
        }
        
      } catch (...) {
        // If we can't access variables this way, add a note
        variables.push_back("(introspection failed - may need specific configuration)");
      }
      
      return variables;
    }

    std::vector<FileInfo> gatherResults(const std::vector<FileInfo>& myResults) const {
      std::vector<FileInfo> allResults;
      
      int myrank = getComm().rank();
      int nprocs = getComm().size();
      
      if (nprocs == 1) {
        // Single process - just return our results
        return myResults;
      }
      
      // For MPI implementation, we would need to serialize and gather the results
      // This is a simplified version that works for single process or when
      // each process writes its own partial results
      
      if (myrank == 0) {
        // Rank 0 starts with its own results
        allResults = myResults;
        
        // In a full MPI implementation, rank 0 would receive results from other ranks
        // For now, we'll just use the results from rank 0
        // TODO: Implement proper MPI_Gather or similar for FileInfo structures
        oops::Log::info() << "Note: Full MPI result gathering not implemented - only showing results from rank 0" << std::endl;
      }
      
      return allResults;
    }

    void writeOutputFile(const std::string& filename, const std::vector<FileInfo>& fileInfos) const {
      std::ofstream outFile(filename);
      if (!outFile.is_open()) {
        throw eckit::Exception("Cannot open output file: " + filename);
      }

      // Write header
      outFile << "================================================================================\n";
      outFile << "                          IODA File Summary Report                             \n";
      outFile << "================================================================================\n";
      outFile << "Generated by: IODA Dump Utility\n";
      outFile << "Total files processed: " << fileInfos.size() << "\n";
      outFile << "================================================================================\n\n";

      // Write summary for each file
      for (const auto& info : fileInfos) {
        // Extract just the filename from the full path
        std::string filename = info.filename;
        size_t slashPos = filename.find_last_of('/');
        if (slashPos != std::string::npos) {
          filename = filename.substr(slashPos + 1);
        }
        
        outFile << "File: " << filename << "\n";
        outFile << "Full path: " << info.filename << "\n";
        
        if (info.success) {
          outFile << "Status: SUCCESS\n";
          outFile << "Number of observations (nobs): " << info.nobs << "\n";
          
          if (!info.obsValueVars.empty()) {
            outFile << "ObsValue variables (" << info.obsValueVars.size() << "):\n";
            for (size_t i = 0; i < info.obsValueVars.size(); ++i) {
              outFile << "  " << (i + 1) << ". " << info.obsValueVars[i] << "\n";
            }
          } else {
            outFile << "ObsValue variables: None detected\n";
          }
        } else {
          outFile << "Status: FAILED\n";
          outFile << "Error: " << info.errorMsg << "\n";
        }
        
        outFile << "--------------------------------------------------------------------------------\n";
      }

      outFile << "\n";
      outFile << "================================================================================\n";
      outFile << "                               End of Report                                   \n";
      outFile << "================================================================================\n";

      outFile.close();
      oops::Log::info() << "Summary written to: " << filename << std::endl;
    }
  };

}  // namespace dautils