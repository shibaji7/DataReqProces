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

from get_sd_data import FetchData
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

def fetch_data(args):
    with open("params/%s.json"%args.param_file_name, "r") as f: j = json.loads("\n".join(f.readlines()))
    o, fo = [], []
    if not os.path.exists(args.tmp_folder): os.system("mkdir " + args.tmp_folder)
    fname = args.tmp_folder + args.file_name_format.format(rad=args.rad, start=args.start.strftime("%Y%m%d%H%M"),\
            end=args.end.strftime("%Y%m%d%H%M")) + args.store_type
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
    print(" Medfilter %d, Scans - %d"%(args.med_filter, len(scans)))
    if args.med_filter and len(scans) > 3:
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
        o[j["keys"]].to_csv(fname, mode="a", header=True, index=False, float_format="%.1f")
        if j["zip"]:
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
        fo[j["keys"]].to_csv(fname, mode="a", header=True, index=False, float_format="%.1f")
        if j["zip"]:
            with open(fname, "rb") as src, open(fname+".bz2", "wb") as dest: dest.write(bz2.compress(src.read()))
            os.system("rm " + fname)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", default=dt.datetime(2010,12,12,12), help="Start date", type=prs.parse)
    parser.add_argument("-e", "--end", default=dt.datetime(2010,12,13), help="End date", type=prs.parse)
    parser.add_argument("-f", "--file_type", default="fitacf", help="File type", type=str)
    parser.add_argument("-r", "--rad", default="gbr", help="Radar code", type=str)
    parser.add_argument("-st", "--store_type", default=".csv", help="Save type", type=str)
    parser.add_argument("-fn", "--file_name_format", default="{rad}_{start}_{end}", help="Save type", type=str)
    parser.add_argument("-pm", "--param_file_name", default="Prikryl", help="Parameter file name", type=str)
    parser.add_argument("-tmp", "--tmp_folder", default="tmp/", help="Local folder to save files", type=str)
    parser.add_argument("-mfil", "--med_filter", default=True, help="Run codes for median filter aswell", type=bool)
    args = parser.parse_args()
    for k in vars(args).keys():
        print("     " + k + "->" + str(vars(args)[k]))
    fetch_data(args)
    os.system("rm -rf .empty __pycache__")
