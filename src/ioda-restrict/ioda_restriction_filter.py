#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import argparse
import yaml
import os
import glob
import shutil
from datetime import datetime, timedelta

OBS_DIM = "Location"

# ----------------------------------------------------------------------
# Compress consecutive index ranges (used only in exprsrd mode)
# ----------------------------------------------------------------------
def compress_ranges(idx_list):
    if not idx_list:
        return []
    ranges = []
    start  = idx_list[0]
    prev   = idx_list[0]

    for x in idx_list[1:]:
        if x == prev + 1:
            prev = x
        else:
            if start == prev:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}-{prev}")
            start = x
            prev  = x
    if start == prev:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}-{prev}")
    return ranges

# ----------------------------------------------------------------------
# Recursively copy groups and variables
# ----------------------------------------------------------------------
def copy_group(in_group, out_group, mask):

    # Copy group-level attributes
    for attr in in_group.ncattrs():
        setattr(out_group, attr, getattr(in_group, attr))

    for var_name, var_in in in_group.variables.items():
        fill_value = getattr(var_in, "_FillValue", None)

        if fill_value is not None:
            var_out = out_group.createVariable(
                var_name, var_in.dtype, var_in.dimensions, fill_value=fill_value
            )
        else:
            var_out = out_group.createVariable(
                var_name, var_in.dtype, var_in.dimensions
            )

        data = var_in[:]

        if mask is None:
            var_out[:] = data
        elif OBS_DIM in var_in.dimensions:
            if data.ndim == 1:
                var_out[:] = data[mask]
            else:
                var_out[:] = data[mask, ...]
        else:
            var_out[:] = data

        for attr in var_in.ncattrs():
            if attr != "_FillValue":
                setattr(var_out, attr, getattr(var_in, attr))

    for grp_name, grp_in in in_group.groups.items():
        grp_out = out_group.createGroup(grp_name)
        copy_group(grp_in, grp_out, mask)

# ----------------------------------------------------------------------
# Copy entire file unchanged
# ----------------------------------------------------------------------
def copy_entire_file(infile, outfile):
    shutil.copy(infile, outfile)
    print(f"  Wrote (unchanged): {outfile}")

# ----------------------------------------------------------------------
# Extract date from path (exprsrd mode)
# ----------------------------------------------------------------------
def extract_date_from_path(input_dir):
    prefixes = ["gfs.", "gdas.", "gcdas."]
    for p in prefixes:
        idx = input_dir.find(p)
        if idx != -1:
            date_start = idx + len(p)
            date_str = input_dir[date_start:date_start + 8]
            return p, date_start, date_str
    raise ValueError(f"Could not find gfs., gdas., or gcdas. in path: {input_dir}")


def get_prev_48h_dir(input_dir):
    prefix, date_start, date_str = extract_date_from_path(input_dir)
    cur_date = datetime.strptime(date_str, "%Y%m%d")
    prev_date = cur_date - timedelta(hours=48)
    prev_date_str = prev_date.strftime("%Y%m%d")

    return (
        input_dir[:date_start]
        + prev_date_str
        + input_dir[date_start + 8:]
    )

# ----------------------------------------------------------------------
# Compute non-restricted mask (exprsrd mode)
# ----------------------------------------------------------------------
def compute_nonrestricted_mask(flag, exp):
    flag = np.ma.array(flag)
    exp  = np.ma.array(exp)

    n = len(flag)
    non_restricted = np.zeros(n, dtype=bool)

    rsrd_missing = np.ma.getmaskarray(flag)
    expr_missing = np.ma.getmaskarray(exp)
    rsrd_present = ~rsrd_missing
    expr_present = ~expr_missing

    missing_either = rsrd_missing | expr_missing
    exception_restricted = rsrd_present & expr_missing
    nr_missing = missing_either & ~exception_restricted
    non_restricted |= nr_missing

    expired = rsrd_present & expr_present & (exp == 48)
    non_restricted |= expired

    return non_restricted

# ----------------------------------------------------------------------
# Mode 1: RSRD filtering (atmos.nr)
# ----------------------------------------------------------------------
def process_rsrd_directory(input_dir, output_dir):
    print(f"[RSRD] Input directory:  {input_dir}")
    print(f"[RSRD] Output directory: {output_dir}")

    os.makedirs(output_dir, exist_ok=True)
    nc_files = sorted(glob.glob(os.path.join(input_dir, "*.nc")))

    if not nc_files:
        print("[RSRD] No .nc files found.")
        return

    for infile in nc_files:
        fname = os.path.basename(infile)
        outfile = os.path.join(output_dir, fname)

        print(f"\n[RSRD] Processing {fname}")

        with Dataset(infile, "r") as nc_in:
            loc_dim = nc_in.dimensions.get(OBS_DIM)

            if loc_dim is None or (loc_dim.isunlimited() and len(loc_dim) == 0):
                print("  No valid Location dimension — copying unchanged.")
                copy_entire_file(infile, outfile)
                continue
            
            md = nc_in.groups.get("MetaData", None)
            if md is None or "restrictionFlag" not in md.variables or "restrictionExpiration" not in md.variables:
                print("  Missing restriction variables — copying unchanged.")
                copy_entire_file(infile, outfile)
                continue

            else:
                flag = md["restrictionFlag"][:]
                exp  = md["restrictionExpiration"][:]

                if flag.size == 0 or exp.size == 0:
                    print("  Restriction arrays zero length — writing empty restricted file.")
                    nloc = len(nc_in.dimensions[OBS_DIM])
                    mask = np.zeros(nloc, dtype=bool)
                else:
                    flag_mask = np.ma.getmaskarray(flag)
                    exp_mask  = np.ma.getmaskarray(exp)
                    mask = flag_mask & exp_mask

            total = mask.size
            kept = int(mask.sum())
            dropped = total - kept

            print(f"  Total obs:   {total}")
            print(f"  Kept obs:    {kept}")
            print(f"  Dropped obs: {dropped}")

            with Dataset(outfile, "w") as nc_out:

                for attr in nc_in.ncattrs():
                    setattr(nc_out, attr, getattr(nc_in, attr))

                for dim_name, dim in nc_in.dimensions.items():
                    if dim_name == OBS_DIM:
                        # Preserve unlimited-ness
                        if dim.isunlimited():
                            nc_out.createDimension(dim_name, None)
                        else:
                            nc_out.createDimension(dim_name, kept)
                    else:
                        nc_out.createDimension(
                            dim_name,
                            None if dim.isunlimited() else len(dim)
                        )

                copy_group(nc_in, nc_out, mask)

            print(f"  Wrote: {outfile}")

# ----------------------------------------------------------------------
# Mode 2: EXPRSRD / Non-restricted filtering (atmos.us)
# ----------------------------------------------------------------------
def process_exprsrd_directory(prev_dir, output_dir):
    print(f"[NonRestrict] Previous directory: {prev_dir}")
    print(f"[NonRestrict] Output directory:   {output_dir}")

    os.makedirs(output_dir, exist_ok=True)
    nc_files = sorted(glob.glob(os.path.join(prev_dir, "*.nc")))

    if not nc_files:
        print("[NonRestrict] No .nc files found.")
        return

    for infile in nc_files:
        fname = os.path.basename(infile)
        outfile = os.path.join(output_dir, fname)

        print(f"\n[NonRestrict] Processing {fname}")

        with Dataset(infile, "r") as nc_in:
            loc_dim = nc_in.dimensions.get(OBS_DIM)
            if loc_dim is None or (loc_dim.isunlimited() and len(loc_dim) == 0):
                print("  No valid Location dimension — copying unchanged.")
                copy_entire_file(infile, outfile)
                continue
            
            md = nc_in.groups.get("MetaData", None)
            nloc = len(nc_in.dimensions[OBS_DIM])

            # Default: fail-closed empty mask
            non_restricted_mask = np.zeros(nloc, dtype=bool)
            kept = 0

            # Missing MetaData or missing restriction variables
            if (
                md is None
                or "restrictionFlag" not in md.variables
                or "restrictionExpiration" not in md.variables
            ):
                print("  Missing restriction variables — copying unchanged.")
                copy_entire_file(infile, outfile)
                continue

            else:

                flag = md["restrictionFlag"][:]
                exp  = md["restrictionExpiration"][:]

                # Fail closed: zero‑length arrays ---
                if flag.size == 0 or exp.size == 0:
                    print("  Restriction arrays zero length — writing empty filtered file.")

                else:
                    # Normal EXPRSRD filtering ---
                    non_restricted_mask = compute_nonrestricted_mask(flag, exp)
                    kept = int(np.sum(non_restricted_mask))
                    dropped = len(flag) - kept

                    print(f"  Total obs:          {len(flag)}")
                    print(f"  Non-restricted obs: {kept}")
                    print(f"  Restricted obs:     {dropped}")

                    print("  Unique RSRD / EXPRSRD patterns:")
                    unique_groups = {}

                    flag_mask = np.ma.getmaskarray(flag)
                    exp_mask  = np.ma.getmaskarray(exp)

                    for i in range(len(flag)):
                        fval = None if flag_mask[i] else int(flag[i])
                        eval = None if exp_mask[i] else int(exp[i])
                        key = (fval, eval)
                        unique_groups.setdefault(key, []).append(i)

                    for (fval, eval), idx_list in unique_groups.items():
                        idx_list_sorted = sorted(idx_list)
                        compressed = compress_ranges(idx_list_sorted)
                        count = len(idx_list_sorted)

                        ftxt = fval if fval is not None else "--"
                        etxt = eval if eval is not None else "--"

                        print(f"    RSRD = {ftxt}, EXPRSRD = {etxt}")
                        print(f"      idx ({count}) = {compressed}")
                        print()

            # If nothing kept, skip writing
            if kept == 0:
                print("  No non-restricted obs — skipping output.")
                continue

            with Dataset(outfile, "w") as nc_out:
                for attr in nc_in.ncattrs():
                    setattr(nc_out, attr, getattr(nc_in, attr))

                for dim_name, dim in nc_in.dimensions.items():
                    if dim_name == OBS_DIM:
                        nc_out.createDimension(dim_name, kept)
                    else:
                        nc_out.createDimension(
                            dim_name,
                            None if dim.isunlimited() else len(dim)
                        )

                copy_group(nc_in, nc_out, non_restricted_mask)

            print(f"  Wrote (non-restricted only): {outfile}")

# ----------------------------------------------------------------------
# Driver logic — RSRD always runs; EXPRSRD is skipped on the WCOSS2 production cluster.
# ----------------------------------------------------------------------
def run_rsrd_exprsrd(stats_yaml):
    with open(stats_yaml, "r") as f:
        stats = yaml.safe_load(f)

    input_dir = stats["input directory"]

    # --- 1. RSRD filter on current cycle ---
    output_nr = os.path.join(os.path.dirname(input_dir), "atmos.nr")
    print("\n=== Running RSRD filter (atmos.nr) ===")
    process_rsrd_directory(input_dir, output_nr)

    # --- Determine whether to run EXPRSRD filter ---
    dev_m = None
    this_m = None

    # Read designated backup machine
    try:
        with open("/lfs/h1/ops/prod/config/prodmachinefile") as f:
            for line in f:
                if "backup" in line:
                    parts = line.strip().split(":")
                    if len(parts) >= 2:
                        dev_m = parts[1]
                    break
    except (FileNotFoundError, OSError):
        print("  Cannot read prodmachinefile.")

    # Read current cluster name
    try:
        with open("/etc/cluster_name") as f:
            this_m = f.read().strip()
    except (FileNotFoundError, OSError):
        print("  Cannot read cluster_name.")

    # Decide whether EXPRSRD should run
    if dev_m is None or this_m is None:
        run_exprsrd = True
    else:
        run_exprsrd = (dev_m == this_m)

    print(f"\nCluster check: dev_m={dev_m}, this_m={this_m}, run_exprsrd={run_exprsrd}")


    # --- 2. EXPRSRD filter on previous 48h cycle ---
    if run_exprsrd:
        prev_dir = get_prev_48h_dir(input_dir)
        output_us = os.path.join(os.path.dirname(prev_dir), "atmos.us")
        print("\n=== Running EXPRSRD filter (atmos.us) ===")
        process_exprsrd_directory(prev_dir, output_us)
    else:
        print("\n=== Skipping EXPRSRD filter (cluster mismatch) ===")

# ----------------------------------------------------------------------
# CLI entry point
# ----------------------------------------------------------------------
def cli():
    parser = argparse.ArgumentParser(
        description="Run both RSRD and EXPRSRD filtering for IODA NetCDF files"
    )
    parser.add_argument(
        "-s", "--stats", required=True,
        help="stats.yaml file created by atmos_bufr_prepobs"
    )
    args = parser.parse_args()
    run_rsrd_exprsrd(args.stats)


# ----------------------------------------------------------------------
# Standard Python entry point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    cli()


