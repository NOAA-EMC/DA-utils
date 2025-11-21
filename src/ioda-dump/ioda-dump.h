#pragma once

#include <iomanip>     // for tables in the output

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

#include "ioda/Engines/ObsStore.h"
#include "ioda/Variables/Variable.h"

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
      }  else if (fullConfig.has("input files")) {
         fullConfig.get("input files", inputFiles);
         oops::Log::info() << "Processing " << inputFiles.size() << " specified files" << std::endl;
      }  else {
         throw eckit::Exception("Either 'input directory' or 'input files' must be specified");
      }
    
      // get "shared path" for mapping and query files
      std::vector<std::string> sharedPath;
      if (fullConfig.has("shared path")) {
         fullConfig.get("shared path", sharedPath);
         oops::Log::info() << "Shared Path: " << sharedPath << std::endl;
      }  else {
         throw std::runtime_error("Missing 'Shared Path' in YAML configuration");
      }

      // get "variables" and "count" from yaml if it has
      std::vector<std::string> previewVars;
      size_t previewCount = 10;  // default if it doesn't have
      size_t indexChannel = 0;  // default if it doesn't have

      if (fullConfig.has("variables")) {
         fullConfig.get("variables", previewVars);
      }
      if (fullConfig.has("count")) {
         fullConfig.get("count", previewCount);
      }
      if (fullConfig.has("channel")) {
         fullConfig.get("channel", indexChannel);
      }

      // logging
      oops::Log::info() << "Preview count set to: " << previewCount << std::endl;
      if (!previewVars.empty()) {
         oops::Log::info() << "Preview variables:" << std::endl;
         for (const std::string& varName : previewVars) {
             oops::Log::info() << "  - " << varName << std::endl;
         }
      }  else {
         oops::Log::info() << "No preview variables specified." << std::endl;
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
             FileInfo info = processFile(file, timeWindow, mycomm,
                                         previewVars, sharedPath, previewCount, indexChannel);    
             fileInfos.push_back(info);
          }  catch (const std::exception& e) {
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
         writeOutputFile(outputFile, allFileInfos, previewVars, previewCount, indexChannel);
      }

      return 0;
    }

   private:
    struct FileInfo {
      std::string filename;
      size_t nobs;
      size_t nchans;
      size_t nrecs;
      std::vector<std::string> metaDataVars;
      std::vector<std::string> obsValueVars;
      std::map<std::string, std::vector<std::string>> previewData;
      bool success;
      std::string errorMsg;
    };

    std::string appname() const {
      return "dautils::IodaDump";
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
            if (ext == ".nc" || ext == ".nc4" || ext == ".h5" || ext == ".hdf5" || ext == ".odb" || ext == ".bufr") {
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
                         const eckit::mpi::Comm & comm,
                         const std::vector<std::string>& previewVars,   
			 const std::vector<std::string>& sharedPath,
                         size_t previewCount,
                         size_t indexChannel) const {                   
      FileInfo info;
      info.filename = filename;
      info.success = false;

      try {
        // create a minimal configuration for this file
        eckit::LocalConfiguration obsConfig;

        obsConfig.set("name", "ioda_dump_obsspace");

        std::string ext = filename.substr(filename.find_last_of('.')+1);
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);     // lowercase

        // Engine Type Selection based on the files
        if (ext == "nc" || ext == "nc4" || ext == "h5" || ext == "hdf5") {
           obsConfig.set("obsdatain.engine.type", "H5File");

        }  else if (ext == "odb") {
           obsConfig.set("obsdatain.obsfile", filename);

           // Create a nested configuration for the engine
           eckit::LocalConfiguration engineConfig;
           engineConfig.set("type", "ODB");

           // Extract instrument name from filename
           std::string base = filename.substr(0, filename.find_last_of('.'));
           std::string instrument = base.substr(base.find_last_of("/\\") + 1);

           std::string yamlDir = sharedPath[0];
           std::string queryFile = yamlDir + "/iodatest_odb_" + instrument + ".yaml";

           engineConfig.set("mapping file", yamlDir + "/odb_default_name_map.yaml");
           engineConfig.set("query file", queryFile);

           // Attach the engine config under obsdatain
           obsConfig.set("obsdatain.engine", engineConfig);

           oops::Log::info() << "IODA-Dump: Engine type = ODB" << std::endl;
           oops::Log::info() << "IODA-Dump: Using query file: " << queryFile << std::endl;

        }  else if (ext == "bufr") {
           obsConfig.set("obsdatain.engine.type", "BUFRFile");
        }  else {
           throw std::runtime_error("Unsupported file extension: " + ext);
        }

        obsConfig.set("obsdatain.engine.obsfile", filename);
        obsConfig.set("simulated variables", std::vector<std::string>{"dummy"});
        obsConfig.set("observed variables", std::vector<std::string>{"dummy"});

        // open the IODA file
        oops::Log::info() << "IODA-Dump: Processing " << filename << std::endl;
        ioda::ObsSpace ospace(obsConfig, comm, timeWindow, comm);

        info.nobs = ospace.nlocs();
        oops::Log::info() << filename << ": nobs = " << info.nobs << std::endl;

        info.nchans = ospace.nchans();
        oops::Log::info() << filename << ": nchans = " << info.nchans << std::endl;

        info.nrecs = ospace.nrecs();
        oops::Log::info() << filename << ": nrecs = " << info.nrecs << std::endl;

        // get ObsValue variables
        std::vector<std::string> allVars = ospace.listVariables();

        // list Variables
        std::cout << "Variables in ospace:" << std::endl;
        for (const auto& var : allVars) {
           std::cout << "  " << var << std::endl;
        }

        // Identify MetaData/ObsValue variables by convention (those that start with "MetaData/ObsValue")
        info.metaDataVars.clear();
        for (const auto& var : allVars) {
            if (var.find("MetaData") == 0) {
                info.metaDataVars.push_back(var);
            }
            if (var.find("ObsValue") == 0) {
                info.obsValueVars.push_back(var);
            }
        }
        oops::Log::info() << filename << ": Found " << info.metaDataVars.size() << " MetaData variables" << std::endl;
        oops::Log::info() << filename << ": Found " << info.obsValueVars.size() << " ObsValue variables" << std::endl;
        info.success = true;

        // Step 1: Extract MetaData/satelliteIdentifier
        std::vector<int> satIDs(info.nobs);
        std::map<int, size_t> satIDCounts;   // satID -> count
        std::vector<std::string> satIDStrings;
        try {
            if (ospace.has("MetaData", "satelliteIdentifier")) {
                ospace.get_db("MetaData", "satelliteIdentifier", satIDs);

                // Count occurrences of each satID
                for (const auto& id : satIDs) {
                   satIDCounts[id]++;
                }
                // Prepare strings for previewData
                for (const auto& kv : satIDCounts) {
                   std::ostringstream oss;
                   oss << kv.first << "(" << kv.second << ")";
                   satIDStrings.push_back(oss.str());
                }
                // Log to console
                std::cout << "Unique satelliteIdentifier values found in " << filename << ":\n";
                for (const auto& kv : satIDCounts) {
                   std::cout << "  " << kv.first << " (" << kv.second << ")\n";
                }
                info.previewData["MetaData/satelliteIdentifier_unique_sorted"] = satIDStrings;
            }
        } catch (const std::exception& e) {
        oops::Log::warning() << "Could not read MetaData/satelliteIdentifier: " << e.what() << std::endl;
        }

        // Extract MetaData/dateTime
        std::vector<int> dateIDs(info.nobs);
        try {
           ospace.get_db("MetaData", "dateTime", dateIDs);
           std::set<int> uniqueSortedDate(dateIDs.begin(), dateIDs.end());
           std::vector<std::string> dateIDStrings;

           // Convert each timestamp to formatted UTC string
           auto formatUnixTime = [](int timestamp) -> std::string {
              std::time_t t = static_cast<std::time_t>(timestamp);
              std::tm tm = *std::gmtime(&t);  // Convert to UTC
              std::ostringstream oss;
              oss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S UTC");
              return oss.str();
           };

           for (const auto& id : uniqueSortedDate) {
               dateIDStrings.push_back(formatUnixTime(id));
           }

           //Store formatted dateTime strings for later output
           info.previewData["MetaData/dateTime_unique_sorted"] = dateIDStrings;

           // Log to console
           if (!dateIDStrings.empty()) {
              oops::Log::info() << "Start dateTime: " << dateIDStrings.front() << "\n";
              oops::Log::info() << "End dateTime:   " << dateIDStrings.back() << "\n";
           } else {
              oops::Log::info() << "No valid dateTime entries found.\n";
           }

           info.previewData["MetaData/dateTime_unique_sorted"] = dateIDStrings;
        } catch (const std::exception& e) {
        oops::Log::warning() << "Could not read MetaData/dateTime: " << e.what() << std::endl;
        }

        // Assume the data is of type 'float'
        size_t nchans = (info.nchans > 0) ? info.nchans : 1;
        size_t totalSize = info.nobs * nchans;

        for (const std::string& varName : previewVars) {
            std::size_t slashPos = varName.find('/');
            if (slashPos == std::string::npos) {
               oops::Log::warning() << "Skipping variable with no group: " << varName << std::endl;
               continue; // Skip this variable
            }

            std::string group = varName.substr(0, slashPos);
            std::string variable = varName.substr(slashPos + 1);

            oops::Log::info() << "Reading data for group: '" << group
                              << "', variable: '" << variable << "'" << std::endl;

            // === 1. Create a vector to hold the data ===
            std::vector<float> obsData(totalSize);

            // === 2. Read the data into your vector === use ospace.get_db() inside a try/catch block
            try {
                // 1D data (no channels)
                // 2D data (with channels), get_db() will read all channels into the 1D vector
                ospace.get_db(group, variable, obsData);
            } catch (const ioda::Exception& e) {
              oops::Log::info() << "Skipping missing variable: " << group << "/" << variable << std::endl;
              continue; // Skip to the next variable
            }

            // === 3. Print values and store them for the report ===
            oops::Log::info() << "--- Displaying data for first " << previewCount
                              << "location of " << varName << " ---" << std::endl;

            const size_t locationsToPrint = previewCount;
            std::vector<std::string> dataLines;

            // Use std::min to avoid errors if you have fewer than 10 locations
            for (size_t i = 0; i < std::min(info.nobs, locationsToPrint); ++i) {
                if (nchans == 1) {                            // --- 1D Case (no channels) ---
                   std::stringstream line;                    // Build the string
                   line << "  Loc[" << i << "]: " << obsData[i];
                   oops::Log::info() << line.str() << std::endl;       // Log it
                   dataLines.push_back(line.str()); // Store it
                }  else {                                      // --- 2D Case (with channels) ---
                   std::stringstream line;
                   line << "  Loc[" << i << "]: ";
                   for (size_t j = 0; j < nchans; ++j) {
                       size_t flat_index = (i * nchans) + j;
                       line << obsData[flat_index] << " ";
                   }
                   oops::Log::info() << line.str() << std::endl; // Log it
                   dataLines.push_back(line.str()); // Store it
                }
            }
            info.previewData[varName] = dataLines; // Add to the map
        }
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

        // Serialize MetaDataVars vector
        size_t metaVarsCount = info.metaDataVars.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&metaVarsCount),
                     reinterpret_cast<const char*>(&metaVarsCount) + sizeof(size_t));
        for (const auto& var : info.metaDataVars) {
          size_t metaVarLen = var.size();
          buffer.insert(buffer.end(), reinterpret_cast<const char*>(&metaVarLen),
                       reinterpret_cast<const char*>(&metaVarLen) + sizeof(size_t));
          buffer.insert(buffer.end(), var.begin(), var.end());
        }

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

        // Serialize previewData map
        size_t mapSize = info.previewData.size();
        buffer.insert(buffer.end(), reinterpret_cast<const char*>(&mapSize),
                     reinterpret_cast<const char*>(&mapSize) + sizeof(size_t));
        for (const auto& pair : info.previewData) {
            // Serialize key (varName)
            size_t keyLen = pair.first.size();
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(&keyLen),
                         reinterpret_cast<const char*>(&keyLen) + sizeof(size_t));
            buffer.insert(buffer.end(), pair.first.begin(), pair.first.end());
            // Serialize value (vector of data lines)
            size_t vecSize = pair.second.size();
            buffer.insert(buffer.end(), reinterpret_cast<const char*>(&vecSize),
                         reinterpret_cast<const char*>(&vecSize) + sizeof(size_t));
            for (const auto& line : pair.second) {
                size_t lineLen = line.size();
                buffer.insert(buffer.end(), reinterpret_cast<const char*>(&lineLen),
                             reinterpret_cast<const char*>(&lineLen) + sizeof(size_t));
                buffer.insert(buffer.end(), line.begin(), line.end());
            }
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

        // Deserialize MetaDataVars vector
        size_t metaVarsCount;
        std::memcpy(&metaVarsCount, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        for (size_t j = 0; j < metaVarsCount; ++j) {
          size_t metaVarLen;
          std::memcpy(&metaVarLen, buffer.data() + pos, sizeof(size_t));
          pos += sizeof(size_t);
          std::string var(buffer.begin() + pos, buffer.begin() + pos + metaVarLen);
          info.metaDataVars.push_back(var);
          pos += metaVarLen;
        }

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

        // Deserialize previewData map
        size_t mapSize;
        std::memcpy(&mapSize, buffer.data() + pos, sizeof(size_t));
        pos += sizeof(size_t);
        for (size_t j = 0; j < mapSize; ++j) {
          // Deserialize key
          size_t keyLen;
          std::memcpy(&keyLen, buffer.data() + pos, sizeof(size_t));
          pos += sizeof(size_t);
          std::string key(buffer.begin() + pos, buffer.begin() + pos + keyLen);
          pos += keyLen;

          // Deserialize value (vector of strings)
          size_t vecSize;
          std::memcpy(&vecSize, buffer.data() + pos, sizeof(size_t));
          pos += sizeof(size_t);
          std::vector<std::string> dataLines;
          for (size_t k = 0; k < vecSize; ++k) {
            size_t lineLen;
            std::memcpy(&lineLen, buffer.data() + pos, sizeof(size_t));
            pos += sizeof(size_t);
            std::string line(buffer.begin() + pos, buffer.begin() + pos + lineLen);
            pos += lineLen;
            dataLines.push_back(line);
          }
          info.previewData[key] = dataLines;
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

    void writeOutputFile(const std::string& filename,
                         const std::vector<FileInfo>& fileInfos,
                         const std::vector<std::string>& previewVars,
                         size_t previewCount,
                         size_t indexChannel) const {
      std::ofstream outFile(filename);
      if (!outFile.is_open()) {
        throw eckit::Exception("Cannot open output file: " + filename);
      }

      // Write header
      outFile << "================================================================================\n";
      outFile << "                          IODA File Dump Report                             \n";
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

          auto satIt = info.previewData.find("MetaData/satelliteIdentifier_unique_sorted");
          if (satIt != info.previewData.end()) {
             outFile << "MetaData/satelliteIdentifier:\n";
             for (const std::string& line : satIt->second) {
                outFile << "    " << line << "\n";
             }
          }

          // print start and end date/Time
          if (info.previewData.find("MetaData/dateTime_unique_sorted") != info.previewData.end()) {
             const std::vector<std::string>& dateIDStrings = info.previewData.at("MetaData/dateTime_unique_sorted");

             if (!dateIDStrings.empty()) {
                outFile << "DateTime Range:\n";
                outFile << "  Start: " << dateIDStrings.front() << "\n";
                outFile << "  End:   " << dateIDStrings.back() << "\n";
             } else {
                outFile << "No valid dateTime entries found.\n";
             }
          }

          outFile << "Number of observations (nobs): " << info.nobs << "\n";
          outFile << "Number of records (nrecs): " << info.nrecs << "\n";
          outFile << "Number of channels (nchans): " << info.nchans << "\n";

          // print out metaData variables w/ size
          if (!info.metaDataVars.empty()) {
             outFile << "MetaData variables (" << info.metaDataVars.size() << "):\n";
             for (size_t i = 0; i < info.metaDataVars.size(); ++i) {
                outFile << "  " << (i + 1) << ". " << info.metaDataVars[i] << "\n";
             }
          }  else {
             outFile << "MetaData variables: None detected\n";
          }

          if (!info.obsValueVars.empty()) {
            outFile << "ObsValue variables (" << info.obsValueVars.size() << "):\n";
            for (size_t i = 0; i < info.obsValueVars.size(); ++i) {
              outFile << "  " << (i + 1) << ". " << info.obsValueVars[i] << "\n";
            }
          } else {
            outFile << "ObsValue variables: None detected\n";
          }

          // table view
          outFile << "Preview Data (" << previewCount << "):\n";
          const int colWidth = 24;

          if (!info.previewData.empty()) {
             // Extract variable names from previewData
             std::vector<std::string> keys;
             std::set<std::string> seenBaseNames;

             for (const auto& pair : info.previewData) {
                const std::string& fullName = pair.first;

                if (fullName == "MetaData/satelliteIdentifier_unique_sorted" ||
                    fullName ==  "MetaData/dateTime_unique_sorted") continue;

                // Strip channel suffix if present
                std::string baseName = fullName;
                size_t atPos = fullName.find("@channel_");
                if (atPos != std::string::npos) {
                    baseName = fullName.substr(0, atPos);
                    std::string indexChannel = fullName.substr(atPos + 9);  // after "@channel_"
                    if (indexChannel != "0") continue;  // skip non-zero channels
                }

                // Avoid duplicates
                if (seenBaseNames.count(baseName) == 0) {
                    keys.push_back(fullName);
                    seenBaseNames.insert(baseName);
                }
             }

             // Determine number of rows to print (use shortest vector)
             size_t nrows = std::numeric_limits<size_t>::max();
             for (const auto& key : keys) {
                nrows = std::min(nrows, info.previewData.at(key).size());
             }

             // Print header (strip group prefix like MetaData/)
             for (const auto& key : keys) {
                size_t slash = key.find('/');
                std::string label = (slash != std::string::npos) ? key.substr(slash + 1) : key;
                outFile << std::setw(colWidth) << std::left << label;
             }
             outFile << "\n";

             // Print aligned rows
             for (size_t i = 0; i < nrows; ++i) {
                for (const auto& key : keys) {
                    const std::vector<std::string>& values = info.previewData.at(key);

                    std::string val = (i < values.size()) ? values[i] : "";

                    // Trim any prefix like "float:" or "int:"
                    size_t colonPos = val.find(':');
                    std::string raw = (colonPos != std::string::npos) ? val.substr(colonPos + 1) : val;

                    // Tokenize the raw string into channel values
                    std::istringstream iss(raw);
                    std::vector<std::string> channels;
                    std::string token;
                    while (iss >> token) {
                           channels.push_back(token);
                    }


                    // Choose which channel to print (default: first channel)
                    if (indexChannel < channels.size()) {
                       outFile << std::setw(colWidth) << std::left << channels[indexChannel];
                    } else {
                       outFile << std::setw(colWidth) << std::left << "";  // empty if not available
                    }
                }
                outFile << "\n";
             }
          } else {
             outFile << "  No preview data stored.\n";
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

