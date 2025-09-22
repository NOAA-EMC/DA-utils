#!/bin/bash
#
# Example usage script for IODA dump utility
# This script demonstrates how to run the utility with different configurations
#

# Set paths (adjust for your environment)
IODA_DUMP_EXE="./bin/iodadump.x"
CONFIG_DIR="./test/testinput"
OUTPUT_DIR="./test/testrun"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "IODA Dump Utility - Example Usage"
echo "=================================="
echo

# Example 1: Single process, file list mode
echo "Example 1: Processing specific files (single process)"
echo "Command: $IODA_DUMP_EXE $CONFIG_DIR/iodadump.yaml"
echo
# $IODA_DUMP_EXE $CONFIG_DIR/iodadump.yaml

# Example 2: MPI mode, file list
echo "Example 2: Processing specific files (MPI, 4 processes)"
echo "Command: mpirun -np 4 $IODA_DUMP_EXE $CONFIG_DIR/iodadump.yaml"
echo
# mpirun -np 4 $IODA_DUMP_EXE $CONFIG_DIR/iodadump.yaml

# Example 3: Directory scanning mode
echo "Example 3: Scanning directory for IODA files"
echo "Command: $IODA_DUMP_EXE $CONFIG_DIR/iodadump_directory.yaml"
echo "Note: Update the directory path in iodadump_directory.yaml first"
echo
# $IODA_DUMP_EXE $CONFIG_DIR/iodadump_directory.yaml

echo "Configuration files:"
echo "  - $CONFIG_DIR/iodadump.yaml (file list mode)"
echo "  - $CONFIG_DIR/iodadump_directory.yaml (directory scan mode)"
echo
echo "Output files will be created in: $OUTPUT_DIR/"
echo
echo "For more information, see src/ioda-dump/README.md"