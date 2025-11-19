#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <dirent.h>
#include <sys/stat.h>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <cstring>

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

  class IodaSummary : public oops::Application {
   public:
    explicit IodaSummary(const eckit::mpi::Comm & comm = oops::mpi::world())
      : Application(comm) {}
    static const std::string classname() {return "dautils::IodaSummary";}

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
          failedInfo.nchans = 0;
          failedInfo.nrecs = 0;
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
      size_t nchans;
      size_t nrecs;
      std::vector<std::string> obsValueVars;
      bool success;
      std::string errorMsg;
    };

    std::string appname() const {
      return "dautils::IodaSummary";
    }

    std::string getCurrentTimeString() const {
      auto now = std::chrono::system_clock::now();
      auto time_t_now = std::chrono::system_clock::to_time_t(now);
      
      std::stringstream ss;
      ss << std::put_time(std::gmtime(&time_t_now), "%Y-%m-%dT%H:%M:%SZ");
      return ss.str();
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
            if (ext == ".nc" || ext == ".nc4" || ext == ".h5" || ext == ".hdf5") {
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
        obsConfig.set("name", "ioda_summary_obsspace");
        obsConfig.set("obsdatain.engine.type", "H5File");
        obsConfig.set("obsdatain.engine.obsfile", filename);
        obsConfig.set("simulated variables", std::vector<std::string>{"dummy"});
        obsConfig.set("observed variables", std::vector<std::string>{"dummy"});

        // open the IODA file
        oops::Log::info() << "IODA-Summary: Processing " << filename << std::endl;
        ioda::ObsSpace ospace(obsConfig, comm, timeWindow, comm);
        
        info.nobs = ospace.nlocs();
        oops::Log::info() << filename << ": nobs = " << info.nobs << std::endl;

        info.nchans = ospace.nchans();
        oops::Log::info() << filename << ": nchans = " << info.nchans << std::endl;

        info.nrecs = ospace.nrecs();
        oops::Log::info() << filename << ": nrecs = " << info.nrecs << std::endl;

        // get ObsValue variables
        std::vector<std::string> allVars = ospace.listVariables();

        // Identify ObsValue variables by convention (those that start with "ObsValue")
        info.obsValueVars.clear();
        for (const auto& var : allVars) {
            if (var.find("ObsValue") == 0) {
                info.obsValueVars.push_back(var);
            }
        }
        oops::Log::info() << filename << ": Found " << info.obsValueVars.size() << " ObsValue variables" << std::endl;
        info.success = true;
        
      } catch (const std::exception& e) {
        info.errorMsg = e.what();
        oops::Log::warning() << "Error processing " << filename << ": " << e.what() << std::endl;
      }

      return info;
    }

    std::vector<FileInfo> gatherResults(const std::vector<FileInfo>& myResults) const {
      std::vector<FileInfo> allResults;
      
      int myrank = getComm().rank();
      int nprocs = getComm().size();
      
      if (nprocs == 1) {
        // Single process - just return our results
        return myResults;
      }
      
      // Serialize local results into a buffer
      std::vector<char> myBuffer = serializeFileInfos(myResults);
      size_t myBufferSize = myBuffer.size();
      
      // Gather buffer sizes from all ranks
      std::vector<size_t> allBufferSizes(nprocs);
      getComm().allGather(myBufferSize, allBufferSizes.begin(), allBufferSizes.end());
      
      // Calculate total size and displacements for allGatherv
      std::vector<int> recvCounts(nprocs);
      std::vector<int> displs(nprocs);
      size_t totalSize = 0;
      for (int i = 0; i < nprocs; ++i) {
        recvCounts[i] = static_cast<int>(allBufferSizes[i]);
        displs[i] = static_cast<int>(totalSize);
        totalSize += allBufferSizes[i];
      }
      
      // Use allGatherv to gather all buffers to all ranks
      std::vector<char> allBuffers(totalSize);
      getComm().allGatherv(myBuffer.begin(), myBuffer.end(),
                          allBuffers.begin(), recvCounts.data(), displs.data());
      
      if (myrank == 0) {
        // Rank 0 deserializes all results
        allResults.clear();
        for (int rank = 0; rank < nprocs; ++rank) {
          if (recvCounts[rank] > 0) {
            std::vector<char> rankBuffer(allBuffers.begin() + displs[rank],
                                        allBuffers.begin() + displs[rank] + recvCounts[rank]);
            std::vector<FileInfo> rankResults = deserializeFileInfos(rankBuffer);
            allResults.insert(allResults.end(), rankResults.begin(), rankResults.end());
          }
        }
      }
      
      return allResults;
    }
    
    std::vector<char> serializeFileInfos(const std::vector<FileInfo>& infos) const {
      std::vector<char> buffer;
      
      // Write number of FileInfo structures
      size_t count = infos.size();
      buffer.insert(buffer.end(), reinterpret_cast<const char*>(&count), 
                   reinterpret_cast<const char*>(&count) + sizeof(size_t));
      
      for (const auto& info : infos) {
        // Serialize filename
        size_t filenameLen = info.filename.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&filenameLen),
                     reinterpret_cast<const char*>(&filenameLen) + sizeof(size_t));
        buffer.insert(buffer.end(), info.filename.begin(), info.filename.end());
        
        // Serialize numeric fields
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&info.nobs),
                     reinterpret_cast<const char*>(&info.nobs) + sizeof(size_t));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&info.nchans),
                     reinterpret_cast<const char*>(&info.nchans) + sizeof(size_t));
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&info.nrecs),
                     reinterpret_cast<const char*>(&info.nrecs) + sizeof(size_t));
        
        // Serialize obsValueVars vector
        size_t varsCount = info.obsValueVars.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&varsCount),
                     reinterpret_cast<const char*>(&varsCount) + sizeof(size_t));
        for (const auto& var : info.obsValueVars) {
          size_t varLen = var.size();
          buffer.insert(buffer.end(), reinterpret_cast<const char*>(&varLen),
                       reinterpret_cast<const char*>(&varLen) + sizeof(size_t));
          buffer.insert(buffer.end(), var.begin(), var.end());
        }
        
        // Serialize success flag
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&info.success),
                     reinterpret_cast<const char*>(&info.success) + sizeof(bool));
        
        // Serialize errorMsg
        size_t errorLen = info.errorMsg.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&errorLen),
                     reinterpret_cast<const char*>(&errorLen) + sizeof(size_t));
        buffer.insert(buffer.end(), info.errorMsg.begin(), info.errorMsg.end());
      }
      
      return buffer;
    }
    
    std::vector<FileInfo> deserializeFileInfos(const std::vector<char>& buffer) const {
      std::vector<FileInfo> infos;
      size_t pos = 0;
      
      // Read number of FileInfo structures
      size_t count;
      std::memcpy(&count, buffer.data() + pos, sizeof(size_t));
      pos += sizeof(size_t);
      
      for (size_t i = 0; i < count; ++i) {
        FileInfo info;
        
        // Deserialize filename
        size_t filenameLen;
        std::memcpy(&filenameLen, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        info.filename = std::string(buffer.begin() + pos, buffer.begin() + pos + filenameLen);
        pos += filenameLen;
        
        // Deserialize numeric fields
        std::memcpy(&info.nobs, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        std::memcpy(&info.nchans, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        std::memcpy(&info.nrecs, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        
        // Deserialize obsValueVars vector
        size_t varsCount;
        std::memcpy(&varsCount, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        for (size_t j = 0; j < varsCount; ++j) {
          size_t varLen;
          std::memcpy(&varLen, buffer.data() + pos, sizeof(size_t));
          pos += sizeof(size_t);
          std::string var(buffer.begin() + pos, buffer.begin() + pos + varLen);
          info.obsValueVars.push_back(var);
          pos += varLen;
        }
        
        // Deserialize success flag
        std::memcpy(&info.success, buffer.data() + pos, sizeof(bool));
        pos += sizeof(bool);
        
        // Deserialize errorMsg
        size_t errorLen;
        std::memcpy(&errorLen, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        info.errorMsg = std::string(buffer.begin() + pos, buffer.begin() + pos + errorLen);
        pos += errorLen;
        
        infos.push_back(info);
      }
      
      return infos;
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
      outFile << "Generated on: " << getCurrentTimeString() << "\n";
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
          outFile << "Number of records (nrecs): " << info.nrecs << "\n";
          outFile << "Number of channels (nchans): " << info.nchans << "\n";

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
