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
    with open("params/%s.json"%param_file_name, "r") as f: j = json.loads("\n".join(f.readlines()))
    return j

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", default=dt.datetime(2011,10,11,12), help="Start date", type=prs.parse)
    parser.add_argument("-e", "--end", default=dt.datetime(2011,10,11,15), help="End date", type=prs.parse)
    parser.add_argument("-r", "--rad", default="gbr", help="Radar code", type=str)
    parser.add_argument("-pm", "--param_file_name", default="Lulu", help="Parameter file name", type=str)
    parser.add_argument("-tmp", "--tmp_folder", default="tmp/", help="Local folder to save files", type=str)
    args = parser.parse_args()
    j = load_param_json(args.param_file_name)
    args.file_type = j["file_type"]
    args.file_name_format = j["file_name_format"]
    for k in vars(args).keys():
        print("     " + k + "->" + str(vars(args)[k]))
    if args.file_type in ["fit", "fitacf"]: fetch_fit_level_data(args)
    if args.file_type in ["map2", "mapex", "cnvmap"]: fetch_map_level_data(args)
    os.system("rm -rf .empty __pycache__")