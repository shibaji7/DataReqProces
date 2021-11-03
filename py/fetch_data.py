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

from get_map_grid_data import FetchMap, to_xarray
from get_fit_data import FetchData
import fit_utils as futils
import boxcar_filter as BF

import json

def fetch_fit_level_data(args):
    args.tmp_folder = args.tmp_folder + args.param_file_name + "/"
    if not os.path.exists(args.tmp_folder): os.system("mkdir -p " + args.tmp_folder)
    fname = args.tmp_folder + args.file_name_format.format(rad=args.rad, start=args.start_date.strftime("%Y%m%d%H%M"),\
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
                    fo.to_csv(fname.replace(args.rad, args.rad+"_mfil"), header=True, index=False, float_format="%g")
                elif args.save_type=="netCDE4":
                    ds = fdata.to_xarray(fscans)
            if args.save_type=="pandas":
                o = fdata.convert_to_pandas(beams, s_params=args.scalers, v_params=args.vectors)
                o = futils.update_geo_location(args.rad, o)
                o.to_csv(fname, header=True, index=False, float_format="%g")
            elif args.save_type=="netCDE4":
                ds = fdata.to_xarray(scans)
        else: print(" Error - Data not found")
    return

def fetch_map_level_data(args):
    args.tmp_folder = args.tmp_folder + args.param_file_name + "/"
    if not os.path.exists(args.tmp_folder): os.system("mkdir -p " + args.tmp_folder)
    dates = args.dates
    if len(dates) == 0:
        dn = args.start_date
        while dn <= args.end_date:
            dates.append(dn)
            dn += dt.timedelta(1)
    if len(args.scalers)+len(args.vectors) > 0: print(f" Fetching additional parameters: {'-'.join(args.scalers+args.vectors)}")
    for hemi in args.hemis:
        fm = FetchMap(dates, hemi, args.file_type, args.filestr)
        if args.data_type in ["map", "map2", "mapex"]: fm.fetch_map_files()
        elif args.data_type in ["cnvmap"]: fm.fetch_cnvmap_files()
        else: print(f" Wrong type found {args.file_type}!")
        for d in dates:
            start = d if args.start_date is None else args.start_date 
            end = d+dt.timedelta(1) if args.end_date is None else args.end_date 
            if args.start_date is None: fname = args.tmp_folder + args.file_name_format.format(hemi=hemi,
                                                                                               dn=start.strftime("%Y%m%d"))
            else: fname = args.tmp_folder + args.file_name_format.format(hemi=hemi, 
                                                                         dn=start.strftime("%Y%m%d.%H%M-")+\
                                                                         end.strftime("%H%M"))
            print(f" File stores - \n\t{fname}")
            if not os.path.exists(fname):
                obj = dict()
                if len(args.scalers)+len(args.vectors) > 0: obj["sv_o"] = fm.get_maps(start, end, args.scalers, args.vectors)
                if ("summary" in args.grid_params.keys()) or ("records" in args.grid_params.keys()):
                    obj["summ_o"], obj["reco_o"] = fm.get_grids(start, end, args.grid_params["summary"], args.grid_params["records"])
                if len(args.pev_params) > 0: obj["pev_o"] = fm.calcFitCnvs(start, end, args.pot_lat_min, args.cores, args.pev_params)
                ds = to_xarray(obj, args.pev_params, args.scalers, args.vectors, args.grid_params)
                ds.to_netcdf(fname)
                os.system("rm -rf raw/*")
    return