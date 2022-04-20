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

import cartopy.crs as ccrs
import cartopy
import sdcarto

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def convert_to_map_lat_lon(xs, ys, _from, _to):
    lat, lon = [], []
    for x, y in zip(xs, ys):
        _lon, _lat = _to.transform_point(x, y, _from)
        lat.append(_lat)
        lon.append(_lon)
    return lat, lon

class MapPlot(object):
    """
    Plot data from map(ex) files
    """

    def __init__(self, rec, ax, hemi="north", maxVelScale=1000.0, min_vel=0.0, map_type="map"):
        self.rec = rec
        self.radEarth = 6371.0
        self.lenFactor = 500.0
        self.radEarthMtrs = self.radEarth * 1000.0
        self.maxVelPlot = maxVelScale
        self.min_vel = min_vel
        self.hemi = hemi
        self.ini_figure(ax)
        self.overlayCnvCntrs()
        return

    def ini_figure(self, ax):
        """
        Instatitate figure and axes labels
        """
        proj = cartopy.crs.NorthPolarStereo() if self.hemi == "north" else cartopy.crs.SouthPolarStereo()
        self.ax = ax
        #self.fig = plt.figure(dpi=300, figsize=(3.5,3.5))
        #self.ax = self.fig.add_subplot(111, projection="sdcarto", map_projection = proj,
        #        coords=self.rec["coords"], plot_date=self.rec["stime"])
        self.ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        if self.hemi == "north": self.ax.set_extent([-180, 180, 50, 90], crs=cartopy.crs.PlateCarree())
        else: self.ax.set_extent([-180, 180, -90, -50], crs=cartopy.crs.PlateCarree())
        plt_lons = np.arange( 0, 361, 15 )
        mark_lons = np.arange( 0, 360, 15 )
        plt_lats = np.arange(40,90,10) if self.hemi == "north" else np.arange(-90,-40,10)
        gl = self.ax.gridlines(crs=cartopy.crs.Geodetic(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        self.ax.mark_latitudes(plt_lats, fontsize="small", color="darkblue")
        self.ax.mark_longitudes(plt_lons, fontsize="small", color="darkblue")
        self.ax.text(0.5, 0.95, self.date_string(), ha="center", va="center",
                transform=self.ax.transAxes, fontsize="medium")
        self.proj = proj
        self.geo = ccrs.Geodetic()
        self.ax.text(-0.02, 0.99, "Coord: MLT", ha="center", va="top",
                transform=self.ax.transAxes, fontsize="x-small", rotation=90)
        model_details = ""
        if "pot.drop" in self.rec.keys(): model_details += r"$\Phi_{pc}=%d$ kV"%(self.rec["rec"]["pot.drop"]/1e3) + "\n"
        if "latmin" in self.rec.keys(): model_details += r"$\Lambda_{HM}=%d^{\circ}$"%np.min(self.rec["rec"]["latmin"]) + "\n"
        if "vector.mlon" in self.rec.keys(): model_details += r"$N_{vc}=%d$"%len(self.rec["rec"]["vector.mlon"]) + "\n"
        if "stid" in self.rec.keys(): model_details += r"$N_{rads}=%d$"%len(self.rec["rec"]["stid"])
        self.ax.text(0.05, 0.05, model_details, ha="left", va="bottom",
                transform=self.ax.transAxes, fontsize="small")
        return

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d/%b/%Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime, etime = self.rec["stime"], self.rec["etime"]
        date_str = "{:{dd} {tt}} -- ".format(stime, dd=dfmt, tt=tfmt)
        if etime.date() == stime.date(): date_str = "{:s}{:{tt}} UT".format(date_str, etime, tt=tfmt)
        else: date_str = "{:s}{:{dd} {tt}} UT".format(date_str, etime, dd=dfmt, tt=tfmt)
        return date_str

    def overlayHMB(self, hmbCol="Gray"):
        xs, ys = 15*aacgmv2.convert_mlt(self.rec["boundary.mlon"], self.rec["stime"], m2a=False),\
                self.rec["boundary.mlat"]
        lat, lon = convert_to_map_lat_lon(xs, ys, self.geo, self.proj)
        self.ax.plot(lon, lat, linewidth=0.8, linestyle="-", color="darkgreen", transform=self.proj)
        self.ax.plot(lon, lat, linewidth=0.4, linestyle="--", color="k", zorder=4.0, transform=self.proj)
        return

    def overlayGridVel(self):
        return
    
    def overlayMapModelVel(self, pltColBar=True, colorBarLabelSize="small", colMap=cm.jet):
        return

    def overlayMapFitVel(self, pltColBar=True, pltModelBar=True, colorBarLabelSize="xx-small",
            colMap=cm.hot ):
        return

    def overlayCnvCntrs(self, zorder=2, line_color="k", line_width=0.6, font_size="small",
            plot_label=True):
        lat_cntr = self.rec["pot"]["lat_cntr"]
        xs = self.rec["pot"]["lon_cntr"]
        pot_arr = self.rec["pot"]["pot_arr"]
        lon_cntr = np.zeros_like(xs)
        for i in range(xs.shape[0]):
            lon_cntr[i,:] = 15*aacgmv2.convert_mlt(xs[i,:], self.rec["stime"], m2a=False)
        XYZ = self.proj.transform_points(self.geo, lon_cntr, lat_cntr)
        p_lat, p_lon = lat_cntr.ravel()[np.argmax(pot_arr)], lon_cntr.ravel()[np.argmax(pot_arr)]
        p_lat, p_lon = convert_to_map_lat_lon([p_lon], [p_lat], self.geo, self.proj)
        self.ax.text(p_lon[0], p_lat[0], "+", fontsize=9)
        n_lat, n_lon = lat_cntr.ravel()[np.argmin(pot_arr)], lon_cntr.ravel()[np.argmin(pot_arr)]
        n_lat, n_lon = convert_to_map_lat_lon([n_lon], [n_lat], self.geo, self.proj)
        self.ax.text(n_lon[0], n_lat[0], "--", fontsize=9)
        cp = self.ax.contour(XYZ[:,:,0], XYZ[:,:,1], pot_arr, colors=line_color, linewidths=line_width,
                locator=LinearLocator(9), transform=self.proj)
        self.ax.clabel(cp, inline=1, fontsize=4, fmt="%d", colors="darkblue")
        return


def txt2csv(fname="tmp/pot.txt", linestart=13):
    """
    Convert potential data to csv for plotting.
    """
    o = []
    with open(fname, "r") as f: lines = f.readlines()[linestart:]
    for l in lines:
        l = list(filter(None, l.split(" ")))
        d = l[9].replace("\n", "")
        x = dict(
                    mlat = float(l[2]),
                    mlon = float(l[3]),
                    EField_north = float(l[4]),
                    EField_east = float(l[5]),
                    Fitted_Vel_North = float(l[6]),
                    Fitted_Vel_East = float(l[7]),
                    Potential = float(l[8]),
                    date = dt.datetime.strptime(d, "%Y-%m-%d/%H:%M:%S")
                )
        o.append(x)
    o = pd.DataFrame.from_records(o)
    return o


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


def parse_NETCDF_file(fname, sdate):
    print(f" Parse NETCDF files {fname}")
    ds = xarray.open_dataset(fname)
    stime = [dt.datetime.utcfromtimestamp( (d - np.datetime64("1970-01-01T00:00:00Z","s")) / np.timedelta64(1, "s"))
            for d in ds.coords["fparam.stime"].values]
    etime = [dt.datetime.utcfromtimestamp( (d - np.datetime64("1970-01-01T00:00:00Z","s")) / np.timedelta64(1, "s"))
            for d in ds.coords["fparam.etime"].values]
    t_id = stime.index(sdate)
    lat, lon = ds.coords["fparam.lat_pot"].values, ds.coords["fparam.lon_pot"].values
    pot = ds.data_vars["fparam.pot_arr"].values[t_id,:,:]
    return lat, lon, pot

def plot_map_grid_level_data(args):
    print(f" Plot operations- {args.plot_type}")
    dat = {}
    if args.nc_file and os.path.exists(args.nc_file): dat["lat"], dat["lon"], dat["pot"] = parse_NETCDF_file(args.nc_file, args.sdate)
    if args.ascii_file and os.path.exists(args.ascii_file): dat["o"] = txt2csv(args.ascii_file)
    geodetic = ccrs.Geodetic()
    orthographic = ccrs.NorthPolarStereo()
    for p in args.plot_types:    
        if p == "pot":
            comp = True if (("o" in dat.keys()) and ("pot" in dat.keys())) else False
            vlim, sep = [-40,40], 21
            rec = {"pot": {}, "stime":args.sdate, "etime": args.edate, "coords":args.coords}
            if comp:
                fig = plt.figure(dpi=300, figsize=(8,4))
                ax0 = fig.add_subplot(121, projection="sdcarto", map_projection = orthographic,
                        coords=args.coords, plot_date=args.sdate)
                ax1 = fig.add_subplot(122, projection="sdcarto", map_projection = orthographic,
                        coords=args.coords, plot_date=args.sdate)
            else: 
                fig = plt.figure(dpi=300, figsize=(4,4))
                ax0 = fig.add_subplot(111, projection="sdcarto", map_projection = orthographic,
                        coords=args.coords, plot_date=args.sdate)
            if comp:
                rec["pot"]["lat_cntr"], rec["pot"]["lon_cntr"], rec["pot"]["pot_arr"] = dat["lat"], dat["lon"], dat["pot"]
                mp = MapPlot(rec, ax0)
                mlat, mlon, mpot = get_gridded_parameters(dat["o"], "mlat", "mlon", "Potential")
                mpot /= 1000.
                rec["pot"]["lat_cntr"], rec["pot"]["lon_cntr"], rec["pot"]["pot_arr"] = mlat, mlon, mpot
                mp = MapPlot(rec, ax1)
            elif ("pot" in dat.keys()):
                rec["pot"]["lat_cntr"], rec["pot"]["lon_cntr"], rec["pot"]["pot_arr"] = dat["lat"], dat["lon"], dat["pot"]
                mp = MapPlot(rec, ax0)
            elif ("o" in dat.keys()):
                mlat, mlon, mpot = get_gridded_parameters(dat["o"], "mlat", "mlon", "Potential")
                mpot /= 1000.
                rec["pot"]["lat_cntr"], rec["pot"]["lon_cntr"], rec["pot"]["pot_arr"] = mlat, mlon, mpot
                mp = MapPlot(rec, ax0)
        fname = "tmp/%s.%s.png"%(args.sdate.strftime("%Y%m%dT%H%M"), p)
        fig.savefig(fname)
        plt.close()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-nc", "--nc_file", default=None, help="netCDF File name", type=str)
    parser.add_argument("-s", "--sdate", default=dt.datetime(2010,10,11,13,34), help="Start date to plot", type=prs.parse)
    parser.add_argument("-e", "--edate", default=dt.datetime(2010,10,11,15), help="End date to plot", type=prs.parse)
    parser.add_argument("-p", "--plot_type", default="pot", help="Plot types", type=str)
    parser.add_argument("-c", "--coords", default="aacgmv2_mlt", help="Coordinate types [aacgmv2, aacgmv2_mlt]", type=str)
    parser.add_argument("-af", "--ascii_file", default="tmp/pot.txt", help="ASCII File", type=str)
    args = parser.parse_args()
    args.plot_types = args.plot_type.split(",")
    for k in vars(args).keys():
        print("     " + k + "->" + str(vars(args)[k]))
    if "pot" in args.plot_types or "efield" in args.plot_types\
            or "vel" in args.plot_types or "grid" in args.plot_types:
        plot_map_grid_level_data(args)
    else: print(f" Invalid plot operation {args.plot_type}!")
    os.system("rm -rf .empty __pycache__")
