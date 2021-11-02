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
        if args.file_type in ["map", "map2", "mapex"]: fm.fetch_map_files()
        elif args.file_type in ["cnvmap"]: fm.fetch_cnvmap_files()
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
#             if not os.path.exists(fmap):
#                 if len(j["scalers"])+len(j["vectors"]) > 0: 
#                     d = fm.get_maps(args.start, args.end, j["scalers"], j["vectors"])
#                     if len(d) == 0: print(f" No data for additional parameters!")
#                     else: d.to_csv(fmap, index=False, header=True, float_format="%g")
#                 else: print(" No additional parameters!")
#             else: print(f" File exists {fmap}!")

                if len(args.pev_params) > 0:
                    o = fm.calcFitCnvs(start, end, args.pot_lat_min, args.cores, args.pev_params)
                    if len(o) == 0: print(f" No data for PEV params!")
                    else:
                        ds = to_xarray(o, args.pev_params)
                        ds.to_netcdf(fname)
                else: print(" No PEV parameters!")
            
#             if not os.path.exists(fgrid_summ) or not os.path.exists(fgrid_reco):
#                 if ("summary" in j["grid_params"].keys() or "records" in j["grid_params"].keys()):
#                     summ, reco = fm.get_grids(args.start, args.end, j["grid_params"]["summary"], j["grid_params"]["records"])
#                     if len(summ) == 0: print(f" No data for additional grid parameters!")
#                     else:
#                         col_map = dict(zip(j["grid_params"]["records"], j["grid_params"]["records_map"]))
#                         summ = summ.rename(columns=col_map)
#                         summ.to_csv(fgrid_summ, index=False, header=True, float_format="%g")
#                     if len(reco) == 0: print(f" No data for additional grid parameters!")
#                     else: reco.to_csv(fgrid_reco, index=False, header=True, float_format="%g")
#                 else: print(" No additional grid parameters!")
#             else: print(f" File exists {fgrid_reco} and(or) {fgrid_summ}!")
    return