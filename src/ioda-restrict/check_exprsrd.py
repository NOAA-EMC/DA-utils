#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import argparse
import yaml
import os
import glob
from datetime import datetime, timedelta

OBS_DIM = "Location"


# ----------------------------------------------------------------------
# Compress consecutive index ranges
# ----------------------------------------------------------------------
def compress_ranges(idx_list):
    if not idx_list:
        return []
    ranges = []
    start = prev = idx_list[0]
    for x in idx_list[1:]:
        if x == prev + 1:
            prev = x
        else:
            if start == prev:
                ranges.append(f"{start}")
            else:
                ranges.append(f"{start}–{prev}")
            start = prev = x
    if start == prev:
        ranges.append(f"{start}")
    else:
        ranges.append(f"{start}–{prev}")
    return ranges


# ----------------------------------------------------------------------
# Copy groups with optional mask
# ----------------------------------------------------------------------
def copy_group(in_group, out_group, mask):
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
def copy_entire_file(nc_in, outfile):
    with Dataset(outfile, "w") as nc_out:
        for attr in nc_in.ncattrs():
            setattr(nc_out, attr, getattr(nc_in, attr))
        for dim_name, dim in nc_in.dimensions.items():
            nc_out.createDimension(
                dim_name,
                None if dim.isunlimited() else len(dim)
            )
        copy_group(nc_in, nc_out, mask=None)
    print(f"  Wrote (unchanged): {outfile}")


# ----------------------------------------------------------------------
# Extract date from path
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


# ----------------------------------------------------------------------
# Build previous (48h earlier) directory
# ----------------------------------------------------------------------
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
# Compute non-restricted mask
# ----------------------------------------------------------------------
def compute_nonrestricted_mask(flag, exp):
    flag = np.ma.array(flag)
    exp  = np.ma.array(exp)

    n = len(flag)
    non_restricted = np.zeros(n, dtype=bool)

    rsrd_missing = flag.mask
    expr_missing = exp.mask
    rsrd_present = ~flag.mask
    expr_present = ~exp.mask

    # Rule 1: missing either (except RSRD present + EXPRSRD missing)
    missing_either = rsrd_missing | expr_missing
    exception_restricted = rsrd_present & expr_missing
    nr_missing = missing_either & ~exception_restricted
    non_restricted |= nr_missing

    # Rule 2: expired (EXPRSRD == 48)
    expired = rsrd_present & expr_present & (exp == 48)
    non_restricted |= expired

    return non_restricted


# ----------------------------------------------------------------------
# Main processing
# ----------------------------------------------------------------------
def process_prev_directory(prev_dir, output_dir):
    print(f"[NonRestrict] Previous directory: {prev_dir}")
    print(f"[NonRestrict] Output directory:   {output_dir}")

    os.makedirs(output_dir, exist_ok=True)
    nc_files = sorted(glob.glob(os.path.join(prev_dir, "*.nc")))

    if not nc_files:
        print("[NonRestrict] No .nc files found.")
        return

    print(f"[NonRestrict] Found {len(nc_files)} NetCDF files.")

    for infile in nc_files:
        fname = os.path.basename(infile)
        outfile = os.path.join(output_dir, fname)

        print(f"\n[Processing] {fname}")

        with Dataset(infile, "r") as nc_in:

            loc_dim = nc_in.dimensions.get(OBS_DIM)
            if loc_dim is None or (loc_dim.isunlimited() and len(loc_dim) == 0):
                print("  No valid Location dimension — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            md = nc_in.groups.get("MetaData", None)
            if md is None:
                print("  No MetaData group — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            if "restrictionFlag" not in md.variables or "restrictionExpiration" not in md.variables:
                print("  Missing restriction variables — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            flag = md["restrictionFlag"][:]
            exp  = md["restrictionExpiration"][:]

#            print("  restrictionFlag sample:", flag[:10])
#            print("  restrictionExpiration sample:", exp[:10])

            non_restricted_mask = compute_nonrestricted_mask(flag, exp)
            restricted_mask = ~non_restricted_mask
            restricted_idx = np.where(restricted_mask)[0]

            total = len(flag)
            kept = int(np.sum(non_restricted_mask))
            dropped = total - kept

            print(f"  Total obs:          {total}")
            print(f"  Non-restricted obs: {kept}")
            print(f"  Restricted obs:     {dropped}")

            # ----------------------------------------------------------
            # Unique RSRD / EXPRSRD patterns (ONLY diagnostic kept)
            # ----------------------------------------------------------
            print("  Unique RSRD / EXPRSRD patterns:")

            unique_groups = {}

            for i in range(len(flag)):
                fval = None if flag.mask[i] else int(flag[i])
                eval = None if exp.mask[i] else int(exp[i])
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


            # ----------------------------------------------------------

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
# Main driver
# ----------------------------------------------------------------------
def main(stats_yaml):
    with open(stats_yaml, "r") as f:
        stats = yaml.safe_load(f)

    input_dir = stats["input directory"]
    prev_dir = get_prev_48h_dir(input_dir)

    output_dir = os.path.join(os.path.dirname(prev_dir), "atmos.us")
    os.makedirs(output_dir, exist_ok=True)

    print(f"[NonRestrict] stats.yaml:   {stats_yaml}")
    print(f"[NonRestrict] input_dir:    {input_dir}")
    print(f"[NonRestrict] prev_dir:     {prev_dir}")
    print(f"[NonRestrict] output_dir:   {output_dir}")

    process_prev_directory(prev_dir, output_dir)


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract non-restricted observations from previous (48h) IODA NetCDF files"
    )
    parser.add_argument(
        "-s", "--stats", required=True,
        help="stats.yaml file created by atmos_bufr_prepobs"
    )
    args = parser.parse_args()
    main(args.stats)

