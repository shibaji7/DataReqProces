#!/usr/bin/env python

"""validate_website.py: Validate based on website data."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import pandas as pd
import datetime as dt
import netCDF4 as nc
from netCDF4 import num2date
import numpy as np

def get_gridded_parameters(
    q, xparam="mlat", yparam="mlon", zparam="Potential"
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(xparam, yparam)
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z

def read_csv_file(fname):
    """
    Read CSV file
    """
    o = []
    with open(fname, "r") as f:
        lines = f.readlines()
        for l in lines[2:]:
            l = list(filter(None, l.replace("\n", "").split(" ")))
            o.append({
                "mlat": float(l[2]),
                "mlon": float(l[3]),
                "EField_north": float(l[4]), 
                "EField_east": float(l[5]), 
                "Fitted_Vel_North": float(l[6]), 
                "Fitted_Vel_East": float(l[7]), 
                "Potential": float(l[8])
            })
    o = pd.DataFrame.from_records(o)
    return o

def read_nc_files(fname, time):
    """
    Read NC file from local
    """
    d = nc.Dataset(fname)
    stime = num2date(
        d.variables["fparam.stime"][:].squeeze(),
        d.variables["fparam.stime"].units,
        calendar=d.variables["fparam.stime"].calendar,
        only_use_cftime_datetimes=False
    )
    stime = [
        dt.datetime.strptime(i.isoformat().split(".")[0],"%Y-%m-%dT%H:%M:%S") 
        for i in stime
    ]
    date_index = stime.index(time)
    pot_arr = d.variables["fparam.pot_arr"][date_index,:,:]
    print(np.min(pot_arr), np.max(pot_arr), np.max(pot_arr)-np.min(pot_arr))
    mlat = d.variables["fparam.lat_pot"][:]
    mlon = d.variables["fparam.lon_pot"][:]
    o = pd.DataFrame()
    o["mlat"], o["mlon"], o["Potential"] = (
        mlat.ravel(), mlon.ravel(), pot_arr.ravel()*1e3
    )
    o = o[
        (o.mlat>=50.) &
        (o.mlon<=358.)
    ]
    return o

def validate(csv_fname, nc_fname, time):
    """
    Validate the errors
    """
    o_csv = read_csv_file(csv_fname)
    o_nc = read_nc_files(nc_fname, time)
    csv_mlon, csv_mlat, csv_pot = get_gridded_parameters(o_csv)
    nc_mlon, nc_mlat, nc_pot = get_gridded_parameters(o_nc)
    diff = np.abs(nc_pot - csv_pot)
    print(f"Date/{time}: Max({'%.3f'%np.max(diff)}), Min({'%.3f'%np.min(diff)}), Mean({'%.3f'%np.mean(diff)})")
    return

validate("tmp/Sai/2202.csv", "tmp/Sai/20170907.north.nc", dt.datetime(2017,9,7,22,2))
validate("tmp/Sai/2342.csv", "tmp/Sai/20170907.north.nc", dt.datetime(2017,9,7,23,42))