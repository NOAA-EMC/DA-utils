#!/usr/bin/env python3
import argparse
import netCDF4 as nc
import numpy as np
import shutil
import os
import sys
import glob
import time
from pathlib import Path

import gsincdiag_to_ioda.gsi_ncdiag as gsid


def run_conv_obs(convfile, outdir, platforms, TotalBias=False):
    print("Processing:"+str(convfile))
    startt = time.time()
    Diag = gsid.Conv(convfile)
    Diag.read()
    Diag.toIODAobs(outdir, platforms=platforms, TotalBias=TotalBias)
    Diag.close()
    print("Time (OBS) %s[%s]: %.3g sec" % (convfile, ",".join(platforms), time.time() - startt))
    return 0


def run_radiances_obs(radfile, outdir, obsbias, TotalBias, qcvars, testrefs):
    print("Processing run_radiances_obs:%s" % radfile)
    startt = time.time()
    Diag = gsid.Radiances(radfile)
    Diag.read()
    Diag.toIODAobs(outdir, obsbias, qcvars, testrefs, TotalBias=TotalBias)
    Diag.close()
    print("Time (OBS) %s: %.3g sec" % (radfile, time.time() - startt))
    return 0


def run_oz_obs(ozfile, outdir, TotalBias=False):
    print("Processing:"+str(ozfile))
    startt = time.time()
    Diag = gsid.Ozone(ozfile)
    Diag.read()
    Diag.toIODAobs(outdir, TotalBias=TotalBias)
    print("Time (OBS) %s: %.3g sec" % (ozfile, time.time() - startt))
    return 0


def run_conv_geo(convfile, outdir):
    print("Processing:"+str(convfile))
    startt = time.time()
    Diag = gsid.Conv(convfile)
    Diag.read()
    Diag.toGeovals(outdir)
    print("Time (GEO) %s: %.3g sec" % (convfile, time.time() - startt))
    return 0


def run_radiances_geo(radfile, outdir):
    print("Processing run_radiances_geo:%s" % radfile)
    startt = time.time()
    Diag = gsid.Radiances(radfile)
    Diag.read()
    Diag.toGeovals(outdir)
    print("Time (GEO) %s: %.3g sec" % (radfile, time.time() - startt))
    return 0


def run_oz_geo(ozfile, outdir):
    print("Processing:"+str(ozfile))
    startt = time.time()
    Diag = gsid.Ozone(ozfile)
    Diag.read()
    Diag.toGeovals(outdir)
    print("Time (GEO) %s: %.3g sec" % (ozfile, time.time() - startt))
    return 0


def run_radiances_obsdiag(radfile, outdir):
    print("Processing run_radiances_obsdiag: %s" % radfile)
    startt = time.time()
    Diag = gsid.Radiances(radfile)
    Diag.read()
    Diag.toObsdiag(outdir)
    print("Time (DIAG) %s: %.3g sec" % (radfile, time.time() - startt))
    return 0

def proc_gsi_ncdiag(ObsDir=False, ObsBias=False, TotalBias=False, QCVars=False, TestRefs=False, DiagDir='./'):
    # process obs files
    if ObsDir:
        if not Path(ObsDir).is_dir():
            raise Exception("Obs dir: '%s' does not exist." % ObsDir)
        # conventional obs first
        # get list of conv diag files
        convfiles = glob.glob(DiagDir+'/*conv*')
        for convfile in convfiles:
            splitfname = convfile.split('/')[-1].split('_')
            if 'conv' in splitfname:
                i = splitfname.index('conv')
                c = "_".join(splitfname[i:i + 2])
            try:
                for p in gsid.conv_platforms[c]:
                    run_conv_obs(convfile, ObsDir, [p], TotalBias=TotalBias)
            except (KeyError, IndexError):
                pass
        # radiances next
        radfiles = glob.glob(DiagDir+'/diag*')
        for radfile in radfiles:
            process = False
            for p in gsid.rad_sensors:
                if p in radfile.split('/')[-1]:
                    process = True
            if process:
                run_radiances_obs(radfile, ObsDir, ObsBias, TotalBias, QCVars, TestRefs)
        # atmospheric composition observations
        # ozone
        for radfile in radfiles:
            process = False
            oz_sensors = gsid.oz_lay_sensors + gsid.oz_lev_sensors
            for p in oz_sensors:
                if p in radfile.split('/')[-1]:
                    process = True
            if process:
                run_oz_obs(radfile, ObsDir, TotalBias=TotalBias)

def combine_ges_anl_ioda(ges_ioda_file, anl_ioda_file, out_ioda_file):
    """Combine GSI IODA ges and anl files into one file."""
    if not os.path.exists(ges_ioda_file):
        raise FileNotFoundError(f"GSI IODA ges file not found: {ges_ioda_file}")
    if not os.path.exists(anl_ioda_file):
        raise FileNotFoundError(f"GSI IODA anl file not found: {anl_ioda_file}")
    
    # copy all groups, variables, and attributes from ges file to out file
    # then append groups and variables from anl file to out file
    # if a group or variable already exists in out file, skip it
    # if a variable exists in both files, the value from the ges file is kept
    with nc.Dataset(ges_ioda_file, 'r') as ges_ds, nc.Dataset(anl_ioda_file, 'r') as anl_ds:
        with nc.Dataset(out_ioda_file, 'w') as out_ds:
            # copy global attributes
            for attr_name in ges_ds.ncattrs():
                out_ds.setncattr(attr_name, ges_ds.getncattr(attr_name))
            # copy dimensions
            for dim_name, dim in ges_ds.dimensions.items():
                out_ds.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))
            for dim_name, dim in anl_ds.dimensions.items():
                if dim_name not in out_ds.dimensions:
                    out_ds.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))
            # copy variables from ges file
            for var_name, var in ges_ds.variables.items():
                out_var = out_ds.createVariable(var_name, var.datatype, var.dimensions)
                # copy variable attributes
                for attr_name in var.ncattrs():
                    out_var.setncattr(attr_name, var.getncattr(attr_name))
                # copy variable data
                out_var[:] = var[:]
            # copy groups from ges file
            for group_name in ges_ds.groups:
                ges_group = ges_ds.groups[group_name]
                out_group = out_ds.createGroup(group_name)
                # copy group attributes
                for attr_name in ges_group.ncattrs():
                    out_group.setncattr(attr_name, ges_group.getncattr(attr_name))
                # copy dimensions
                for dim_name, dim in ges_group.dimensions.items():
                    if dim_name not in out_group.dimensions:
                        out_group.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))
                # copy variables
                for var_name, var in ges_group.variables.items():
                    out_var = out_group.createVariable(var_name, var.datatype, var.dimensions)
                    # copy variable attributes
                    for attr_name in var.ncattrs():
                        if attr_name != '_FillValue':  # avoid issues with _FillValue
                            out_var.setncattr(attr_name, var.getncattr(attr_name))
                    # copy variable data
                    out_var[:] = var[:]
            # copy variables from anl file
            for var_name, var in anl_ds.variables.items():
                if var_name not in out_ds.variables:
                    out_var = out_ds.createVariable(var_name, var.datatype, var.dimensions)
                    # copy variable attributes
                    for attr_name in var.ncattrs():
                        out_var.setncattr(attr_name, var.getncattr(attr_name))
                    # copy variable data
                    out_var[:] = var[:]
            # copy groups from anl file
            for group_name in anl_ds.groups:
                if group_name not in out_ds.groups:
                    anl_group = anl_ds.groups[group_name]
                    out_group = out_ds.createGroup(group_name)
                    # copy group attributes
                    for attr_name in anl_group.ncattrs():
                        out_group.setncattr(attr_name, anl_group.getncattr(attr_name))
                    # copy dimensions
                    for dim_name, dim in anl_group.dimensions.items():
                        if dim_name not in out_group.dimensions:
                            out_group.createDimension(dim_name, (len(dim) if not dim.isunlimited() else None))
                    # copy variables
                    for var_name, var in anl_group.variables.items():
                        if var_name not in out_group.variables:
                            out_var = out_group.createVariable(var_name, var.datatype, var.dimensions)
                            # copy variable attributes
                            for attr_name in var.ncattrs():
                                if attr_name != '_FillValue':  # avoid issues with _FillValue
                                    out_var.setncattr(attr_name, var.getncattr(attr_name))
                            # copy variable data
                            out_var[:] = var[:]
            # create ombg and oman variables in out file
            ges_group_hofx = ges_ds.groups.get('GsiHofXBcGes')
            ges_group_obs = ges_ds.groups.get('ObsValue')
            anl_group_hofx = anl_ds.groups.get('GsiHofXBcAnl')
            for var_name, var in ges_group_hofx.variables.items():
                # get observed value
                obvalue = ges_group_obs.variables[var_name][:]
                # compute and save ombg
                out_group_bg = out_ds.createGroup("ombg")
                out_var_bg = out_group_bg.createVariable(var_name, var.datatype, var.dimensions)
                # copy variable attributes
                for attr_name in var.ncattrs():
                    if attr_name != '_FillValue':  # avoid issues with _FillValue
                        out_var_bg.setncattr(attr_name, var.getncattr(attr_name))
                out_var_bg[:] = obvalue - var[:]
                # compute and save oman
                out_group_an = out_ds.createGroup("oman")
                out_var_an = out_group_an.createVariable(var_name, var.datatype, var.dimensions)
                # copy variable attributes
                for attr_name in var.ncattrs():
                    if attr_name != '_FillValue':  # avoid issues with _FillValue
                        out_var_an.setncattr(attr_name, var.getncattr(attr_name))
                    # copy variable data
                out_var_an[:] = obvalue - anl_group_hofx.variables[var_name][:]

        print(f"Combined GSI IODA file created: {out_ioda_file}")

if __name__ == "__main__":
    ScriptName = os.path.basename(sys.argv[0])
    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("input_dir", help="Path to concatenated GSI diag files")
    ap.add_argument("-o", "--obs_dir",
                    help="Path to directory to output observations")
    ap.add_argument("-g", "--geovals_dir",
                    help="Path to directory to output observations")
    ap.add_argument("-d", "--obsdiag_dir",
                    help="Path to directory to output observations")
    ap.add_argument("-b", "--add_obsbias", default=False, action='store_true',
                    help="Add ObsBias group to output observations")
    ap.add_argument("--add_total_bias", default=False, action='store_true',
                    help="Add bias generated from adjusted minus unadjusted to output observations")
    ap.add_argument("-q", "--add_qcvars", default=False, action='store_true',
                    help="Add QC variables to output observations")
    ap.add_argument("-r", "--add_testrefs", default=False, action='store_true',
                    help="Add TestReference group to output observations")

    MyArgs = ap.parse_args()

    DiagDir = MyArgs.input_dir

    proc_gsi_ncdiag(ObsDir=MyArgs.obs_dir, ObsBias=MyArgs.add_obsbias,
                    TotalBias=MyArgs.add_total_bias, QCVars=MyArgs.add_qcvars,
                    TestRefs=MyArgs.add_testrefs, DiagDir=DiagDir)
