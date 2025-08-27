"""fetch_data.py: Module is dedicated of data fetching using pydarn module and fittofitacf (binary) modules"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import datetime as dt
import argparse
from dateutil import parser as prs
import bz2
import pydarn
import gc

from get_map_grid_data import FetchMap, to_xarray
from get_fit_data import FetchData
import fit_utils as futils
import boxcar_filter as BF
import pandas as pd

import json

def check_folders(args):
    args.tmp_folder_store = args.tmp_folder_store.format(user=args.param_file_name)
    if not os.path.exists(args.tmp_folder_store): os.system("mkdir -p " + args.tmp_folder_store)
    args.tmp_folder_check = args.tmp_folder_check.format(user=args.param_file_name)
    if not os.path.exists(args.tmp_folder_check): os.system("mkdir -p " + args.tmp_folder_check)
    return args

def fetch_fit_level_data(args):
    args = check_folders(args)
    fname = args.tmp_folder_store + args.file_name_format.format(rad=args.rad, start=args.start_date.strftime("%Y%m%d%H%M"),\
            end=args.end_date.strftime("%Y%m%d%H%M"))
    if not os.path.exists(fname):
        fdata = FetchData( args.rad, [args.start_date, args.end_date], args.file_type )
        if "fit" == args.file_type: fdata.convert2fitacf("raw/")
        beams, scans = fdata.fetch_data(by="scan")
        print(" Length(scans) - %d"%(len(scans)))
        os.system("rm -rf raw/*.fit*")
        if len(scans) > 0:
            if args.boxcar:
                fscans = BF.doFilter(scans, thresh=args.med_filt["thresh"], cores=args.cores) 
                if args.save_type=="pandas":
                    fbeams = []
                    for s in fscans:
                        fbeams.extend(s.beams)
                    fo = fdata.convert_to_pandas(fbeams, s_params=args.scalers, v_params=args.vectors)
                    fo = futils.update_geo_location(args.rad, fo)
                    fo.to_csv(fname.replace(args.rad, args.rad+"_mfil"), header=True, index=False,
                              float_format="%g")
                elif args.save_type=="netCDE4":
                    ds = fdata.to_xarray(fscans)
            if args.save_type=="pandas":
                o = fdata.convert_to_pandas(beams, s_params=args.scalers, v_params=args.vectors)
                o = futils.update_geo_location(args.rad, o)
                o.to_csv(fname, header=True, index=False, float_format="%g")
            elif args.save_type=="netCDF4":
                ds = fdata.to_xarray(scans)
                ds.to_netcdf(fname)
        else: print(" Error - Data not found")
    return

def fetch_map_level_data(args):
    args = check_folders(args)
    dates = args.dates
    if len(dates) == 0:
        dn = args.start_date
        while dn <= args.end_date:
            dates.append(dn)
            dn += dt.timedelta(1)
    if len(args.scalers)+len(args.vectors) > 0: 
        print(f" Fetching additional parameters: {'-'.join(args.scalers+args.vectors)}")
    for hemi in args.hemis:
        for d in dates:
            try:
                fm = FetchMap([d], hemi, args.file_type, args.filestr)
                if args.data_type in ["map", "map2", "mapex"]: fm.fetch_map_files()
                elif args.data_type in ["cnvmap"]: fm.fetch_cnvmap_files()
                else: print(f" Wrong type found {args.file_type}!")
                start = d if args.start_date is None else args.start_date 
                end = d+dt.timedelta(1) if args.end_date is None else args.end_date 
                if args.start_date is None: 
                    fname_store = args.tmp_folder_store + args.file_name_format.format(
                        hemi=hemi, dn=start.strftime("%Y%m%d")
                    )
                    fname_check = args.tmp_folder_check + args.file_name_format.format(
                        hemi=hemi, dn=start.strftime("%Y%m%d")
                    )
                else: 
                    fname_store = args.tmp_folder_store + args.file_name_format.format(
                        hemi=hemi, dn=start.strftime("%Y%m%d.%H%M-")+end.strftime("%H%M")
                    )
                    fname_check = args.tmp_folder_check + args.file_name_format.format(
                        hemi=hemi, dn=start.strftime("%Y%m%d.%H%M-")+end.strftime("%H%M")
                    )
                if (not os.path.exists(fname_store)) and (not os.path.exists(fname_check)):
                    obj = dict()
                    if len(args.scalers)+len(args.vectors) > 0: obj["sv_o"] = fm.get_maps(start, end, 
                                                                                          args.scalers,
                                                                                          args.vectors)
                    if ("summary" in args.grid_params.keys()) or ("records" in args.grid_params.keys()):
                        obj["summ_o"], obj["reco_o"] = fm.get_grids(start, end, args.grid_params["summary"],
                                                                    args.grid_params["records"])
                    if len(args.pev_params) > 0: obj["pev_o"] = fm.calcFitCnvs(start, end, args.pot_lat_min,
                                                                               args.cores, args.pev_params, 
                                                                               args.plots)
                    ds = to_xarray(obj, args.pev_params, args.scalers, args.vectors, args.grid_params)
                    ds.to_netcdf(fname_store)
                    if ".csv" in fname_store:
                        to_csv(fname_store, obj, args.pev_params)
                    os.system("rm -rf raw/*")
            except:
                import traceback
                traceback.print_exc()
                print(f" Exception occured - {d}")
            gc.collect()
    return

def to_csv(fname, obj, pev_params):
    """
    Convert obj files to csv files
    """
    if "sv_o" in obj:
        obj["sv_o"].to_csv(fname.replace(".csv", "-sv.csv"), index=False, header=True)
    if "summ_o" in obj:
        obj["summ_o"].to_csv(fname.replace(".csv", "-sum.csv"), index=False, header=True)
        obj["reco_o"].to_csv(fname.replace(".csv", "-rec.csv"), index=False, header=True)
    if "pev_o" in obj:
        o = pd.DataFrame()
        for i in obj["pev_o"]:
            x = pd.DataFrame()
            L = len(i["pot"]["pot_arr"].ravel())
            x["start_time"], x["end_time"], x["hemi_str"], x["coords"] = (
                [i["stime"]]*L,
                [i["etime"]]*L,
                [i["hemi_str"]]*L,
                [i["coords"]]*L,
            )
            x["pot(V)"], x["lat_cntr"], x["lon_cntr"], x["mlt_cntr"] = (
                i["pot"]["pot_arr"].ravel()*1e3,
                i["pot"]["lat_cntr"].ravel(),
                i["pot"]["lon_cntr"].ravel(),
                i["pot"]["mlt_cntr"].ravel()
            )
            o = pd.concat([o, x])
        o.to_csv(fname.replace(".csv", "-pev.csv"), index=False, header=True)
    return