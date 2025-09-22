#include "iodadump.h"
#include "oops/runs/Run.h"

// This application reads IODA files and dumps basic summary information
// to a formatted ASCII text file.
// It supports MPI execution and can process directories or lists of files.

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  dautils::IodaDump iodadump;
  return run.execute(iodadump);
}