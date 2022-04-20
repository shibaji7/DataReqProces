"""fetch.py: Module is dedicated of data fetching using pydarn module"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys
sys.path.extend(["py/"])
import os
import datetime as dt
import argparse
from dateutil import parser as prs

import json

from fetch_data import fetch_fit_level_data, fetch_map_level_data

def load_param_json(param_file_name):
    print(f" Read param file params/{param_file_name}.json")
    with open("params/%s.json"%param_file_name, "r") as f: o = json.loads("\n".join(f.readlines()))
    # Set datetime objects from json
    if o["start_date"] is not None: o["start_date"] = prs.parse(o["start_date"])
    if o["end_date"] is not None: o["end_date"] = prs.parse(o["end_date"])
    if (o["dates"] is not None) and (len(o["dates"]) > 0): o["dates"] = [prs.parse(d) for d in o["dates"]]
    # Set save file type
    if ".nc" in o["file_name_format"]: o["save_type"] = "netCDF4"
    elif ".csv" in o["file_name_format"]: o["save_type"] = "pandas"
    # Set median filter
    o["boxcar"] = True if ("med_filt" in o.keys()) and ("thresh" in o["med_filt"].keys())\
            and (o["med_filt"]["thresh"] is not None) else False
    if "plots" not in o.keys(): o["plots"] = {}
    return o

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", default=None, help="Start date", type=prs.parse)
    parser.add_argument("-e", "--end", default=None, help="End date", type=prs.parse)
    parser.add_argument("-r", "--rad", default="gbr", help="Radar code", type=str)
    parser.add_argument("-f", "--file_type", default=None, help="File type other than fitacf or map2", type=str)
    parser.add_argument("-pm", "--param_file_name", default="Lei", help="Parameter file name", type=str)
    parser.add_argument("-tmp", "--tmp_folder", default="tmp/", help="Local folder to save files", type=str)
    args = parser.parse_args()
    o = load_param_json(args.param_file_name)
    args.file_type = o["data_type"] if args.file_type is None else args.file_type
    for k in list(o.keys()):
        setattr(args, k, o[k])
    for k in vars(args).keys():
        print("     " + k + "->" + str(vars(args)[k]))
    if args.data_type in ["fitacf"]: fetch_fit_level_data(args)
    elif args.data_type in ["map", "grid", "cnvmap", "mapex", "map2"]: fetch_map_level_data(args)
    os.system("rm -rf .empty __pycache__")
