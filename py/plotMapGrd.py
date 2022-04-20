#!/usr/bin/env python

"""plotMapGrd.py: module to plot map plots."""

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
import matplotlib as mpl
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys
sys.path.extend(["py/"])

import xarray
import os
import datetime as dt
import argparse
from dateutil import parser as prs
import numpy as np
import scipy as sp
import pandas as pd
import aacgmv2

import cartopy.crs as ccrs
import cartopy
import sdcarto

import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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
    
    def __init__(self, rec, hemi="north", maxVelScale=1000.0, min_vel=0.0, map_type="map"):
        self.rec = rec
        self.radEarth = 6371.0
        self.lenFactor = 500.0
        self.radEarthMtrs = self.radEarth * 1000.0
        self.maxVelPlot = maxVelScale
        self.min_vel = min_vel
        self.hemi = hemi
        self.ini_figure()
        return
    
    def ini_figure(self):
        """
        Instatitate figure and axes labels
        """
        proj = cartopy.crs.NorthPolarStereo() if self.hemi == "north" else cartopy.crs.SouthPolarStereo()
        
        self.fig = plt.figure(dpi=300, figsize=(3.5,3.5))
        self.ax = self.fig.add_subplot(111, projection="sdcarto", map_projection = proj,
                                  coords=self.rec["coords"], plot_date=self.rec["stime"])
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
        if "pot.drop" in self.rec["rec"].keys(): 
            model_details += r"$\Phi_{pc}=%d$ kV"%(self.rec["rec"]["pot.drop"]/1e3) + "\n"
        if "latmin" in self.rec["rec"].keys():
            model_details += r"$\Lambda_{HM}=%d^{\circ}$"%np.min(self.rec["rec"]["latmin"]) + "\n"
        if "vector.mlon" in self.rec["rec"].keys():
            model_details += r"$N_{vc}=%d$"%len(self.rec["rec"]["vector.mlon"]) + "\n"
        if "stid" in self.rec["rec"].keys():
            model_details += r"$N_{rads}=%d$"%len(self.rec["rec"]["stid"])
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
        xs, ys = 15*aacgmv2.convert_mlt(self.rec["rec"]["boundary.mlon"], self.rec["stime"], m2a=False),\
                    self.rec["rec"]["boundary.mlat"]
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
        norm = mpl.colors.Normalize(self.min_vel, self.maxVelPlot)
        mlats_plot = self.rec["vel_efield"]["mlats"]
        mlons_plot = self.rec["vel_efield"]["mlons"]
        vel_mag = self.rec["vel_efield"]["vel_mag"]
        vel_azm = self.rec["vel_efield"]["vel_azm"]
        self.mapFitPltStrt = []
        self.mapFitPltVec = []
        for nn, nn_mlats in enumerate(mlats_plot):
            start_lat, start_lon = mlats_plot[nn], mlons_plot[nn]
            start_lon = 15*aacgmv2.convert_mlt(start_lon, self.rec["stime"], m2a=False)
            vec_len = vel_mag[nn] * self.lenFactor / self.radEarth / 1000.0
            end_lat = np.arcsin(np.sin(np.deg2rad(nn_mlats)) * np.cos(vec_len) +
                               np.cos(np.deg2rad(nn_mlats)) * np.sin(vec_len) *
                               np.cos(np.deg2rad(vel_azm[nn])))
            end_lat = np.degrees(end_lat)
            del_lon = np.arctan2(np.sin(np.deg2rad(vel_azm[nn])) *
                                 np.sin(vec_len) * np.cos(np.deg2rad(nn_mlats)),
                                 np.cos(vec_len) - np.sin(np.deg2rad(nn_mlats))
                                 * np.sin(np.deg2rad(end_lat)))
            end_lon = mlons_plot[nn] + np.degrees(del_lon)
            end_lon = 15*aacgmv2.convert_mlt(end_lon, self.rec["stime"], m2a=False)
            ylat_start, xlon_start = convert_to_map_lat_lon([start_lon], [start_lat], self.geo, self.proj)
            ylat_end, xlon_end = convert_to_map_lat_lon([end_lon], [end_lat], self.geo, self.proj)
            self.mapFitPltStrt.append(self.ax.scatter(xlon_start, ylat_start, 
                                                        c=[vel_mag[nn]], s=2.0,
                                                        vmin=self.min_vel,
                                                        vmax=self.maxVelPlot, 
                                                        alpha=0.7, cmap=colMap,
                                                        zorder=5.0,
                                                        edgecolor="none"))
            map_color = colMap(norm(vel_mag[nn]))
            self.mapFitPltVec.append(self.ax.plot([xlon_start, xlon_end],
                                                    [ylat_start, ylat_end], 
                                                    color=map_color, lw=0.4))
        if pltColBar:
            ax = inset_axes(self.ax, width="3%", height="15%", loc="upper left")
            cbar = mpl.pyplot.colorbar(self.mapFitPltStrt[0], cax=ax, orientation="vertical")
            cbar.ax.tick_params(labelsize=colorBarLabelSize) 
            vlabel = "Vel [$m s^{-1}$]"
            cbar.set_label(vlabel, size=colorBarLabelSize)
            v = 1000. * self.lenFactor / self.radEarth / 1000.0
            self.ax.plot([0.025,0.025+v], [0.75,0.75], transform=self.ax.transAxes, 
                         ls="-", color="darkred", lw=0.8)
            self.ax.scatter([0.025],[0.75], transform=self.ax.transAxes, s=2, color="darkred", alpha=0.7,
                           zorder=5.0)
            self.ax.text(0.025, 0.76, r"1000 $ms^{-1}$", ha="left", va="bottom", 
                         transform=self.ax.transAxes, fontsize="xx-small")
        if pltModelBar:
            By, Bz = self.rec["rec"]["IMF.By"], self.rec["rec"]["IMF.Bz"]
            #lim = 15*(2+int(np.max(np.abs([By, Bz]))/15))
            lim = 20
            ax = inset_axes(self.ax, width="15%", height="20%", loc="upper right")
            ax.axhline(0, ls="-", lw=0.4, color="k")
            ax.axvline(0, ls="-", lw=0.4, color="k")
            ax.set_ylim(-lim,lim)
            ax.set_xlim(-lim,lim)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.patch.set_alpha(0.0)
            if lim < 1000: 
                w = 0.05
                ax.arrow(x=0, y=0, dx=By, dy=Bz, width=w, head_width=w*7, 
                         head_length=w*5, color="darkred")
            else: lim = 10.
            model_txt = "OMNI IMF\n"+ "Stat Mod: RG96\n"\
                        + self.rec["rec"]["model.angle"] + ", " + r"$%s$"%self.rec["rec"]["model.level"]
            self.ax.text(0.9, 0.75, model_txt, ha="center", va="top", 
                     transform=self.ax.transAxes, fontsize=4, color="darkred")
            ax.text(0.9, 0.4, "+Y", ha="center", va="center", 
                    transform=ax.transAxes, fontsize=3.5, color="darkred")
            ax.text(0.6, 0.95, "+Z (%d nT)"%lim, ha="left", va="center", 
                    transform=ax.transAxes, fontsize=3.5, color="darkred")
            for key, spine in ax.spines.items():
                spine.set_visible(False)
            ax.set_xticks([])
            ax.set_yticks([])
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
    
    def set_radars(self):
        import pydarn
        date = self.rec["stime"]
        for stid in self.rec["rec"]["stid"]:
            self.ax.overlay_radar(stid)
            #self.ax.overlay_fov(stid)
        return
    
    def save(self, fname):
        self.fig.savefig(fname, bbox_inches="tight")
        return

if __name__ == "__main__":
    rec = dict(stime=dt.datetime(2016,1,1), etime=dt.datetime(2016,1,1,0,2), coords="aacgmv2_mlt")
    MapPlot(rec).save("out.png")
