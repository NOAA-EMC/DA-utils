#include "iodastats.h"
#include "oops/runs/Run.h"

// This application preprocesses observation space
// statistics and writes them to a new file for
// concatenation and future use/plotting

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  dautils::IodaStats iodastats;
  return run.execute(iodastats);
}
