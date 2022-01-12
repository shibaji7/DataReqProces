#!/usr/bin/env python

"""plotlib.py: utility module for Costom Carto py geoaxes to plot data on aacgmv2 coordinates."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use(["science", "ieee"])

import sys
sys.path.extend(["py/"])

import xarray
import os
import datetime as dt
import argparse
from dateutil import parser as prs
import numpy as np
import pandas as pd
import aacgmv2

from get_fit_data import txt2csv

import cartopy.crs as ccrs
import cartopy
import sdcarto

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def get_gridded_parameters(q, xparam, yparam, zparam, r=0, rounding=False):
    """
    Method convert dataframe to 2D array.
    """
    plotParamDF = q[ [xparam, yparam, zparam] ]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    else:
        plotParamDF[xparam] = plotParamDF[xparam].tolist()
        plotParamDF[yparam] = plotParamDF[yparam].tolist()
    plotParamDF = plotParamDF.groupby( [xparam, yparam] ).mean().reset_index()
    plotParamDF = plotParamDF[ [xparam, yparam, zparam] ].pivot( xparam, yparam )
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y  = np.meshgrid( x, y )
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
            np.isnan(plotParamDF[zparam].values),
            plotParamDF[zparam].values)
    return X,Y,Z.T

def _add_colorbar(fig, ax, bounds, colormap, label=""):
    """
    Add a colorbar to the right of an axis.
    """
    import matplotlib as mpl
    pos = ax.get_position()
    cpos = [pos.x1 + 0.06, pos.y0 + 0.0125,
            0.015, pos.height * 0.5]                # this list defines (left, bottom, width, height
    cax = fig.add_axes(cpos)
    norm = matplotlib.colors.BoundaryNorm(bounds, colormap.N)
    cb2 = matplotlib.colorbar.ColorbarBase(cax, cmap=colormap,
                                           norm=norm,
                                           ticks=bounds,
                                           spacing="uniform",
                                           orientation="vertical")
    cb2.set_label(label)
    return

def add_axes(fig, num, proj, coords, date, vlim, nlim, ty):
    ax = fig.add_subplot(num, projection="sdcarto",\
                         map_projection = proj,\
                         coords=coords, plot_date=date)
    ax.overaly_coast_lakes()
    ax.set_extent([-180, 180, 50, 90], crs=cartopy.crs.PlateCarree())
    
    bounds = np.linspace(vlim[0], vlim[1], nlim)
    cmap = matplotlib.cm.get_cmap("jet")
    cmap.set_bad("w", alpha=0.0)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    
    plt_lons = np.arange( 0, 361, 30 )
    mark_lons = np.arange( 0, 360, 30 )
    plt_lats = np.arange(30,90,10)
    gl = ax.gridlines(crs=cartopy.crs.Geodetic(), linewidth=0.5)
    gl.xlocator = mticker.FixedLocator(plt_lons)
    gl.ylocator = mticker.FixedLocator(plt_lats)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.n_steps = 90
    ax.mark_latitudes(plt_lats, fontsize=10, color="r")
    ax.mark_longitudes(plt_lons, fontsize=10, color="darkblue")
    
    ax.text(0.01, 1.1, "Date: "+date.strftime("%Y-%m-%d")+"\n"+\
             "Time: %s UT"%(args.sdate.strftime("%H:%M")+"-"+(date+dt.timedelta(minutes=2)).strftime("%H:%M")),
             ha="left", va="center", transform=ax.transAxes)
    ax.text(0.99, 1.1, "Coords: $%s$"%coords.replace("_", "-"), ha="right", va="center", transform=ax.transAxes)
    ax.text(1.05, 0.99, ty, ha="center", va="top", transform=ax.transAxes)
    return ax, bounds, norm, cmap

def plot_map_grid_level_data(args):
    print(f" Plot operations- {args.plot_type}")
    ds = xarray.open_dataset(args.nc_file)
    stime = [dt.datetime.utcfromtimestamp( (d - np.datetime64("1970-01-01T00:00:00Z","s")) / np.timedelta64(1, "s"))
            for d in ds.coords["fparam.stime"].values]
    etime = [dt.datetime.utcfromtimestamp( (d - np.datetime64("1970-01-01T00:00:00Z","s")) / np.timedelta64(1, "s"))
            for d in ds.coords["fparam.etime"].values]
    t_id = stime.index(args.sdate)
    geodetic = ccrs.Geodetic()
    orthographic = ccrs.NorthPolarStereo()
    for p in args.plot_types:    
        if p == "pot":
            vlim, sep = [-40,40], 21
            if args.comp is not None: 
                fig = plt.figure(dpi=150, figsize=(8,4))
                ax0, _, _, _ = add_axes(fig, 121, orthographic, args.coords, args.sdate, vlim, sep, "py")
                ax1, bounds, norm, cmap = add_axes(fig, 122, orthographic, args.coords, args.sdate, vlim, sep, "web")
            else: 
                fig = plt.figure(dpi=150, figsize=(4,4))
                ax0, bounds, norm, cmap = add_axes(fig, 111, orthographic, args.coords, args.sdate, vlim, sep, "py")
            
            lat, lon = ds.coords["fparam.lat_pot"].values, ds.coords["fparam.lon_pot"].values
            
            if args.coords == "geo":
                lon = np.mod( (lon + 180), 360 ) - 180
                glat, glon = np.zeros_like(lat), np.zeros_like(lon)
                for i in range(lat.shape[0]):
                    glat[i,:], glon[i,:], _ = aacgmv2.convert_latlon_arr(lat[i,:], lon[i,:], 300, args.sdate ,"A2G")
                lat, lon = glat, glon
            pot = ds.data_vars["fparam.pot_arr"].values[t_id,:,:]
            XYZ = orthographic.transform_points(geodetic, lon, lat)
            ax0.contourf(XYZ[:,:,0], XYZ[:,:,1], pot, cmap=cmap, vmax=vlim[1], 
                        vmin=vlim[0], transform=orthographic, alpha=0.8)
            if args.comp is None: _add_colorbar(fig, ax0, bounds, cmap, r"Potential ($\Phi_{pc}$), kV")
            else:
                o = txt2csv(args.comp)
                mlat, mlon, mpot = get_gridded_parameters(o, "mlat", "mlon", "Potential")
                mpot /= 1000.
                XYZ = orthographic.transform_points(geodetic, mlon, mlat)
                ax1.contourf(XYZ[:,:,0], XYZ[:,:,1], mpot, cmap=cmap, vmax=vlim[1], 
                        vmin=vlim[0], transform=orthographic, alpha=0.8)
                _add_colorbar(fig, ax1, bounds, cmap, r"Potential ($\Phi_{pc}$), kV")
            
        fname = "tmp/%s.%s.png"%(args.sdate.strftime("%Y%m%dT%H%M"), p)
        fig.savefig(fname)
        plt.close()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-nc", "--nc_file", default="tmp/Lulu/20101011.1200-1500.north.nc", help="netCDF File name")
    parser.add_argument("-s", "--sdate", default=dt.datetime(2010,10,11,13,34), help="Start date to plot", type=prs.parse)
    parser.add_argument("-e", "--edate", default=dt.datetime(2010,10,11,15), help="End date to plot", type=prs.parse)
    parser.add_argument("-p", "--plot_type", default="pot", help="Plot types", type=str)
    parser.add_argument("-c", "--coords", default="aacgmv2", help="Coordinate types [geo, aacgmv2, aacgmv2_mlt]", type=str)
    parser.add_argument("-cp", "--comp", default=None, help="Comparison file", type=str)
    args = parser.parse_args()
    args.plot_types = args.plot_type.split(",")
    for k in vars(args).keys():
        print("     " + k + "->" + str(vars(args)[k]))
    if "pot" in args.plot_types or "efield" in args.plot_types\
            or "vel" in args.plot_types or "grid" in args.plot_types:
        plot_map_grid_level_data(args)
    else: print(f" Invalid plot operation {args.plot_type}!")
    os.system("rm -rf .empty __pycache__")