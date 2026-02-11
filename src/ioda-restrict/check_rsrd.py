#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import argparse
import yaml
import os
import glob


OBS_DIM = "Location"   # confirmed by user


# ----------------------------------------------------------------------
# Helper: recursively copy groups and variables
# ----------------------------------------------------------------------
def copy_group(in_group, out_group, mask):
    """
    Recursively copy all variables and subgroups from in_group to out_group.
    If mask is None, copy data unchanged.
    """

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

        # No mask → copy unchanged
        if mask is None:
            var_out[:] = data

        # Masked copy
        elif OBS_DIM in var_in.dimensions:
            if data.ndim == 1:
                var_out[:] = data[mask]
            else:
                var_out[:] = data[mask, ...]
        else:
            var_out[:] = data

        # Copy attributes except _FillValue
        for attr in var_in.ncattrs():
            if attr == "_FillValue":
                continue
            setattr(var_out, attr, getattr(var_in, attr))

    # Recursively copy subgroups
    for grp_name, grp_in in in_group.groups.items():
        grp_out = out_group.createGroup(grp_name)
        copy_group(grp_in, grp_out, mask)


# ----------------------------------------------------------------------
# Helper: copy entire file unchanged
# ----------------------------------------------------------------------
def copy_entire_file(nc_in, outfile):
    with Dataset(outfile, "w") as nc_out:

        # Copy global attributes
        for attr in nc_in.ncattrs():
            setattr(nc_out, attr, getattr(nc_in, attr))

        # Copy dimensions
        for dim_name, dim in nc_in.dimensions.items():
            nc_out.createDimension(
                dim_name,
                None if dim.isunlimited() else len(dim)
            )

        # Copy all groups and variables unchanged
        copy_group(nc_in, nc_out, mask=None)

    print(f"  Wrote: {outfile}")


# ----------------------------------------------------------------------
# Main: process all files in a directory
# ----------------------------------------------------------------------
def process_directory(input_dir, output_dir):
    print(f"[ChkRestriction] Input directory:  {input_dir}")
    print(f"[ChkRestriction] Output directory: {output_dir}")

    os.makedirs(output_dir, exist_ok=True)

    nc_files = sorted(glob.glob(os.path.join(input_dir, "*.nc")))

    if not nc_files:
        print("[ChkRestriction] No .nc files found.")
        return

    print(f"[ChkRestriction] Found {len(nc_files)} NetCDF files.")

    for infile in nc_files:
        fname = os.path.basename(infile)
        outfile = os.path.join(output_dir, fname)

        print(f"\n[Processing] {fname}")

        with Dataset(infile, "r") as nc_in:

            # ----------------------------------------------------------
            # Case 1: Location dimension missing or zero-length
            # ----------------------------------------------------------
            loc_dim = nc_in.dimensions.get(OBS_DIM)

            if loc_dim is None:
                print("  No Location dimension — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            if loc_dim.isunlimited() and len(loc_dim) == 0:
                print("  Location dimension UNLIMITED with 0 length — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            # ----------------------------------------------------------
            # Case 2: Missing restriction variables
            # ----------------------------------------------------------
            md = nc_in.groups.get("MetaData", None)

            if (md is None or
                "restrictionFlag" not in md.variables or
                "restrictionExpiration" not in md.variables):

                print("  Missing restriction variables — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            # ----------------------------------------------------------
            # Case 3: restriction arrays exist but have zero length
            # ----------------------------------------------------------
            flag = md["restrictionFlag"][:]
            exp  = md["restrictionExpiration"][:]

            if flag.size == 0 or exp.size == 0:
                print("  Restriction arrays have zero length — copying unchanged.")
                copy_entire_file(nc_in, outfile)
                continue

            # ----------------------------------------------------------
            # Case 4: Normal filtering path
            # ----------------------------------------------------------
            print("  restrictionFlag unique:", np.unique(flag))
            print("  restrictionExpiration unique:", np.unique(exp))

            mask = flag.mask & exp.mask

            total = len(mask)
            kept = np.sum(mask)
            dropped = total - kept

            print(f"  Total obs:   {total}")
            print(f"  Kept obs:    {kept}")
            print(f"  Dropped obs: {dropped}")

            # Write filtered file
            with Dataset(outfile, "w") as nc_out:

                # Copy global attributes
                for attr in nc_in.ncattrs():
                    setattr(nc_out, attr, getattr(nc_in, attr))

                # Copy dimensions (Location shrinks)
                for dim_name, dim in nc_in.dimensions.items():
                    if dim_name == OBS_DIM:
                        nc_out.createDimension(dim_name, kept)
                    else:
                        nc_out.createDimension(
                            dim_name,
                            None if dim.isunlimited() else len(dim)
                        )

                # Copy groups + variables with mask
                copy_group(nc_in, nc_out, mask)

            print(f"  Wrote: {outfile}")

def main(stats_yaml):
    with open(stats_yaml, "r") as f:
        stats = yaml.safe_load(f)

    input_dir = stats["input directory"]
    output_dir = os.path.join(os.path.dirname(input_dir), "atmos.nr")

    os.makedirs(output_dir, exist_ok=True)

    process_directory(input_dir, output_dir)


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter IODA NetCDF files using stats.yaml")
    parser.add_argument("-s", "--stats", required=True, help="stats.yaml file created by atmos_bufr_prepobs")
    args = parser.parse_args()
    main(args.stats)
