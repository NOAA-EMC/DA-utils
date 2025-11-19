# ioda-dump.x — IODA file summary dumper

A small utility that reads IODA observation files (HDF5/NetCDF) and writes a formatted ASCII summary. It can scan a directory or process an explicit list of files and supports MPI execution.

Summary contents per file:
- Number of observations (nobs)
- Number of records (nrecs)
- Number of channels (nchans)
- List of variables in the `ObsValue/` group (by convention, variables with names starting with `ObsValue`)

The summary is written to a single text file you choose.

## Build

This program is part of `da-utils` and is built with Ecbuild/CMake.

- Target name: `ioda-dump.x`
- Language standard: C++17
- Linked libraries: `oops`, `ioda` (and their transitive deps, eckit, etc.)

It is automatically included by the top-level `CMakeLists.txt` via:
```
add_subdirectory( ioda-dump )
```

If you are building the whole bundle, the executable will appear under your build tree, typically `build/bin/ioda-dump.x` after a successful `make`.

## Usage

Basic usage is to pass a YAML configuration file:

```
./ioda-dump.x path/to/config.yaml
```

MPI is supported; for example with 4 processes:

```
mpirun -np 4 ./ioda-dump.x path/to/config.yaml
```

Files are distributed across ranks in a round-robin fashion, and all results are automatically gathered to rank 0 for output to the summary file.

## Configuration file (YAML)

The program reads a single YAML file. Required keys:

- `time window`: begin/end time used when opening the IODA ObsSpace
- Exactly one of:
  - `input directory`: directory to scan (non-recursive) for IODA files
  - `input files`: list of explicit file paths
- `output file`: path to the ASCII report to write (ensure the directory exists)

### Example: directory scan mode

```yaml
# Time window for ObsSpace (ISO-8601, UTC)
time window:
  begin: 2000-01-01T00:00:00Z
  end:   2030-12-31T23:59:59Z
  # optional; include lower bound
  bound to include: begin

# Scan a directory (non-recursive). Files with extensions .nc, .nc4, .h5, .hdf5 are included.
input directory: /path/to/ioda/files

# Output summary path (parent directory must exist)
output file: /path/to/output/ioda_summary.txt
```

### Example: explicit file list

```yaml
# Time window for ObsSpace
time window:
  begin: 2025-01-01T00:00:00Z
  end:   2025-12-31T23:59:59Z

# Provide files explicitly
input files:
  - /data/ioda/sst/file1.nc
  - /data/ioda/sst/file2.h5
  - /data/ioda/sst/file3.nc4

# Output summary path (parent directory must exist)
output file: /work/ioda_summary.txt
```

## Output format

A single ASCII file is produced containing a header and one section per input file, for example:

```
================================================================================
                          IODA File Summary Report                             
================================================================================
Generated on: 2025-09-30T12:34:56Z
Total files processed: 3
================================================================================

File: file1.nc
Full path: /path/to/ioda/file1.nc
Status: SUCCESS
Number of observations (nobs): 123456
Number of records (nrecs): 123456
Number of channels (nchans): 5
ObsValue variables (2):
  1. ObsValue/sst
  2. ObsValue/sea_ice_fraction
--------------------------------------------------------------------------------
... (repeated per file)

================================================================================
                               End of Report                                   
================================================================================
```

If a file fails to open or parse, the section shows `Status: FAILED` and the error message.

## Notes and limitations

- MPI distribution: input files are divided round-robin across ranks, and results are automatically gathered to rank 0 for output.
- Directory scanning is non-recursive and includes only regular files with extensions `.nc`, `.nc4`, `.h5`, `.hdf5`.
- Time window filtering is applied by `ioda::ObsSpace` if time metadata are present in the file.
- Ensure the parent directory of `output file` exists; otherwise opening the output will fail.

## Troubleshooting

- "Directory does not exist": check the `input directory` path.
- "Cannot open output file": create the parent directory or fix permissions.
- "Failed to process file ...": the file might not be an IODA v2 file, may be corrupt, or incompatible with the compiled IODA version.
- Build/link errors: verify your environment includes matching versions of `eckit`, `oops`, `ioda`, HDF5/NetCDF, and the required compilers/MPIs.

## Implementation details

- Backend: `ioda::ObsSpace` opened with the H5File engine.
- The program enumerates `ospace.listVariables()` and reports those with names beginning `ObsValue`.
- Counts reported are `ospace.nlocs()` (nobs), `ospace.nrecs()` (nrecs), and `ospace.nchans()` (nchans).
