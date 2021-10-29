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
import rad_fov
import boxcar_filter as BF

import json

def update_geo_location(rad, o):
    hdw = pydarn.read_hdw_file(rad)
    sgate = 0
    egate = hdw.gates
    sbeams, ebeam = 0, hdw.beams
    rfov = rad_fov.CalcFov(hdw=hdw, ngates=egate)
    glat, glon = [], []
    for i, rec in o.iterrows():
        glat, glon = rfov.latFull[rec["bmnum"], rec["slist"]], rfov.lonFull[rec["bmnum"], rec["slist"]]
    o["glat"], o["glon"] = glat, glon
    return o

def fetch_fit_level_data(args):
    args.tmp_folder = args.tmp_folder + args.param_file_name + "/"
    with open("params/%s.json"%args.param_file_name, "r") as f: j = json.loads("\n".join(f.readlines()))
    o, fo = [], []
    if not os.path.exists(args.tmp_folder): os.system("mkdir -p " + args.tmp_folder)
    fname = args.tmp_folder + args.file_name_format.format(rad=args.rad, start=args.start.strftime("%Y%m%d%H%M"),\
            end=args.end.strftime("%Y%m%d%H%M"))
    if not os.path.exists(fname):
        if "fit" == args.file_type:
            fdata = FetchData( args.rad, [args.start, args.end], args.file_type )
            fdata.convert2fitacf(args.tmp_folder)
            beams, scans = fdata.fetch_data(by="scan")
            o = fdata.convert_to_pandas(beams)
            os.system("rm -rf " + args.tmp_folder + "*.fit*")
        elif "fitacf" in args.file_type:
            fdata = FetchData( args.rad, [args.start, args.end], args.file_type )
            beams, scans = fdata.fetch_data(by="scan")
            o = fdata.convert_to_pandas(beams)
        else: print(" Error - Data not found")
        print(" Length(scans) - %d"%(len(scans)))
        if len(j["med_filt"]) and len(scans) > 3:
            print(" Into median filtering ...")
            fscans = BF.doFilter(scans, thresh=j["med_filt"]["thresh"], cores=j["med_filt"]["cores"]) 
            fbeams = []
            for s in fscans:
                fbeams.extend(s.beams)
            fo = fdata.convert_to_pandas(fbeams)
        if len(o) > 0:
            print(" Saving to file -", fname)
            o = update_geo_location(args.rad, o)
            if os.path.exists(fname): os.system("rm " + fname)
            if j["print_desc"]:
                txt = "=======================================\n"
                txt += " Parameter Description:\n"
                txt += "".join([" %s - %s\n"%(k, j["desc"][k]) for k in j["keys"]])
                txt += "=======================================\n"
                with open(fname, "w") as f: f.write(txt)
            o[j["keys"]].to_csv(fname, mode="a", header=True, index=False, float_format="%g")
            with open(fname, "rb") as src, open(fname+".bz2", "wb") as dest: dest.write(bz2.compress(src.read()))
            os.system("rm " + fname)
        if len(fo) > 0: 
            fname = fname.replace(args.rad , args.rad+"_mf")
            print(" Saving to file -", fname)
            fo = update_geo_location(args.rad, fo)
            if os.path.exists(fname): os.system("rm " + fname)
            if j["print_desc"]:
                txt = "=======================================\n"
                txt += " Parameter Description:\n"
                txt += "".join([" %s - %s\n"%(k, j["desc"][k]) for k in j["keys"]])
                txt += " mFilter thresh - %.2f \n"%j["med_filt"]["thresh"]
                txt += "=======================================\n"
                with open(fname, "w") as f: f.write(txt)
            fo[j["keys"]].to_csv(fname, mode="a", header=True, index=False, float_format="%g")
            with open(fname, "rb") as src, open(fname+".bz2", "wb") as dest: dest.write(bz2.compress(src.read()))
            os.system("rm " + fname)
    return

def fetch_map_level_data(args):
    args.tmp_folder = args.tmp_folder + args.param_file_name + "/"
    dn = args.start
    dates = []
    while dn <= args.end:
        dates.append(dn)
        dn += dt.timedelta(1)
    if not os.path.exists(args.tmp_folder): os.system("mkdir -p " + args.tmp_folder)
    with open("params/%s.json"%args.param_file_name, "r") as f: j = json.loads("\n".join(f.readlines()))
    fparam = args.tmp_folder + args.file_name_format.format(pev_params="_".join(j["pev_params"]), hemi=j["hemi"],
                                                            start=args.start.strftime("%Y%m%d%H%M"),
                                                            end=args.end.strftime("%Y%m%d%H%M"))
    fmap = args.tmp_folder + args.file_name_format.format(pev_params="map", hemi=j["hemi"],
                                                            start=args.start.strftime("%Y%m%d%H%M"),
                                                            end=args.end.strftime("%Y%m%d%H%M")).replace(".nc", ".csv")
    fgrid_summ = args.tmp_folder + args.file_name_format.format(pev_params="grid_summ", hemi=j["hemi"],
                                                            start=args.start.strftime("%Y%m%d%H%M"),
                                                            end=args.end.strftime("%Y%m%d%H%M")).replace(".nc", ".csv")
    fgrid_reco = args.tmp_folder + args.file_name_format.format(pev_params="grid_reco", hemi=j["hemi"],
                                                            start=args.start.strftime("%Y%m%d%H%M"),
                                                            end=args.end.strftime("%Y%m%d%H%M")).replace(".nc", ".csv")
    
    if len(j["scalers"])+len(j["vectors"]) > 0: print(f" Fetching additional parameters: {'-'.join(j['scalers']+j['vectors'])}")
    print(f" File stores - \n\t{fparam}\n\t{fmap}\n\t{fgrid_summ}\n\t{fgrid_reco}")
    if not os.path.exists(fparam) or not os.path.exists(fmap)\
            or not os.path.exists(fgrid_summ) or not os.path.exists(fgrid_reco): 
        fm = FetchMap(dates, j["hemi"], args.file_type, j["filestr"])
        if args.file_type in ["map", "map2", "mapex"]: fm.fetch_map_files()
        elif args.file_type == "cnvmap": fm.fetch_cnvmap_files()
        else: print(f" Wrong type found {args.file_type}!")
            
        if not os.path.exists(fmap):
            if len(j["scalers"])+len(j["vectors"]) > 0: 
                d = fm.get_maps(args.start, args.end, j["scalers"], j["vectors"])
                if len(d) == 0: print(f" No data for additional parameters!")
                else: d.to_csv(fmap, index=False, header=True, float_format="%g")
            else: print(" No additional parameters!")
        else: print(f" File exists {fmap}!")
            
        if not os.path.exists(fparam):
            if len(j["pev_params"]) > 0:
                o = fm.calcFitCnvs(args.start, args.end, j["pot_lat_min"], j["cores"], j["pev_params"])
                if len(o) == 0: print(f" No data for PEV params!")
                else:
                    ds = to_xarray(o, j["pev_params"])
                    ds.to_netcdf(fparam)
            else: print(" No PEV parameters!")
        else: print(f" File exists {fparam}!")
            
        if not os.path.exists(fgrid_summ) or not os.path.exists(fgrid_reco):
            if ("summary" in j["grid_params"].keys() or "records" in j["grid_params"].keys()):
                summ, reco = fm.get_grids(args.start, args.end, j["grid_params"]["summary"], j["grid_params"]["records"])
                if len(summ) == 0: print(f" No data for additional grid parameters!")
                else:
                    col_map = dict(zip(j["grid_params"]["records"], j["grid_params"]["records_map"]))
                    summ = summ.rename(columns=col_map)
                    summ.to_csv(fgrid_summ, index=False, header=True, float_format="%g")
                if len(reco) == 0: print(f" No data for additional grid parameters!")
                else: reco.to_csv(fgrid_reco, index=False, header=True, float_format="%g")
            else: print(" No additional grid parameters!")
        else: print(f" File exists {fgrid_reco} and(or) {fgrid_summ}!")
    return