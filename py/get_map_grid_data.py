#!/usr/bin/env python

"""get_map_grid_data.py: module is dedicated to fetch map2, mapex, grid2, grd, gridex data from files."""

__author__ = "Chakraborty, S."
__copyright__ = "Copyright 2020, SuperDARN@VT"
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
import numpy as np
import scipy
import pandas as pd
import datetime as dt
import glob
import bz2
import gzip
import pydarn
import pydarnio
import configparser
import shutil
import xarray

from plotMapGrd import MapPlot

import multiprocessing as mp
from functools import partial

class FetchMap(object):
    """
    Fetch map level data [map, mapex, cnvmap]
    """

    def __init__(self, dates, hemi, file_type="map2", 
                 _filestr="/sd-data/{year}/{file_type}/{hemi}/{date}.{hemi}.{file_type}.bz2",
                radEarth=6371.0, lenFactor=500.):
        self.dates = dates
        self.hemi = hemi
        self.file_type = file_type
        self._filestr = _filestr
        
        # set up some initial parameters
        self.radEarth = radEarth
        # This is used to change the length of the vector on the plot
        self.lenFactor = lenFactor
        self.radEarthMtrs = self.radEarth * 1000.0
        self.records = None
        return
    
    def fetch_map_files(self):
        """
        Read mapex and map2 files
        """
        self.files = []
        for d in self.dates:
            f = self._filestr.format(year=d.year, hemi=self.hemi,file_type=self.file_type,
                                     date=d.strftime("%Y%m%d"))
            fs = glob.glob(f)
            if len(fs) > 0: self.files.append(fs[0])
            else: print(f" File not exists, {f}!")
        return
    
    def fetch_cnvmap_files(self):
        """
        Read and copy cnvmaps
        """
        if not os.path.exists("raw/"): os.system("mkdir raw/")
        self.files = []
        for d in self.dates:
            f = self._filestr.format(year=d.year, hemi=self.hemi,file_type=self.file_type, 
                                     date=d.strftime("%Y%m%d"))
            fs = glob.glob(f)
            if len(fs) > 0:
                f = fs[0]
                shutil.copy(f, "raw/")
                dest = "raw/" + f.split("/")[-1]
                self.files.append(dest.replace(".bz2", ""))
                os.system("bzip2 -d " + dest)
            else: print(f" File not exists, {f}!")
        return
    
    def fetch_records(self):
        if self.records == None:
            self.records = []
            for f in self.files:
                if ("cnvmap" in f) or ("mapex" in f):
                    reader = pydarn.SuperDARNRead()
                    recs = reader.read_dmap(f)
                else:
                    with bz2.open(f) as fp: ds = fp.read()
                    reader = pydarnio.SDarnRead(ds, True)
                    recs = reader.read_map()
                self.records.extend(recs)
            if self.file_type == "cnvmap": os.system("rm -rf raw/*")
        return self.records
    
    def get_grids(self, start, end, summary=[], records=[]):
        """
        Fetch gridex, grid2 content
        """
        print(" Fetch grid records.")
        self.summ, self.reco = pd.DataFrame(), pd.DataFrame()
        grids = self.fetch_records()
        for r in grids:
            stime = dt.datetime(r["start.year"], r["start.month"], r["start.day"], r["start.hour"],
                                r["start.minute"], int(r["start.second"]))
            etime = dt.datetime(r["end.year"], r["end.month"], r["end.day"], r["end.hour"],
                                r["end.minute"], int(r["end.second"]))
            o = pd.DataFrame(r, columns=summary)
            o["stime"], o["etime"] = stime, etime
            self.summ = pd.concat([self.summ, o])
            if "vector.mlat" in r: 
                o = pd.DataFrame(r, columns=records)
                o["stime"], o["etime"] = stime, etime
                self.reco = pd.concat([self.reco, o])
        self.summ = self.summ.reset_index().drop(columns=["index"])
        self.reco = self.reco.reset_index().drop(columns=["index"])
        self.summ = self.summ[(self.summ.stime>=start) & (self.summ.stime<=end)]
        self.reco = self.reco[(self.reco.stime>=start) & (self.reco.stime<=end)]
        return self.summ, self.reco
    
    def get_maps(self, start, end, scalers=["pot.drop"], vectors=[]):
        """
        Fetch mapex, map2 file content
        """
        print(" Fetch map records.")
        self.reco = pd.DataFrame()
        records = self.fetch_records()
        for r in records:
            stime = dt.datetime(r["start.year"], r["start.month"], r["start.day"], r["start.hour"],
                                r["start.minute"], int(r["start.second"]))
            etime = dt.datetime(r["end.year"], r["end.month"], r["end.day"], r["end.hour"],
                                r["end.minute"], int(r["end.second"]))
            if len(vectors)>0: o = pd.DataFrame(r, columns=vectors)
            else: o = pd.DataFrame()
            L = 1 if len(o) == 0 else len(o)
            o["stime"], o["etime"] = [stime]*L, [etime]*L
            for p in scalers:
                o[p] = [r[p]]*L
            self.reco = pd.concat([self.reco, o])
        self.reco = self.reco.sort_values(by="stime")
        self.reco = self.reco[(self.reco.stime>=start) & (self.reco.stime<=end)].reset_index().\
                    drop(columns=["index"])
        return self.reco
    
    def calcFitCnvVel(self, rec):
        """
        Calculate fitted convection velocity magnitude and azimuth from
        map data (basically coefficients of the fit)
        """
        stime, etime, hemi, r = rec["stime"], rec["etime"], rec["hemi"], rec["rec"]
        
        if "vector.mlat" in r:
            hemi_str = "north" if hemi==1 else "south"

            # get the standard location/LoS(grid) Vel parameters.
            mlats, mlons = r["vector.mlat"], r["vector.mlon"] 
            vels, azms = r["vector.vel.median"], r["vector.kvect"]
            
            # Some important parameters from fitting.
            coeff_fit = np.array(r["N+2"])
            order_fit = r["fit.order"]
            lat_shft_fit = r["lat.shft"]
            lon_shft_fit = r["lon.shft"]
            lat_min_fit = r["latmin"]

            # Set up more parameters for getting the fitted vectors
            # the absolute part is for the southern hemisphere
            theta = np.deg2rad(90.0 - np.absolute(mlats))
            theta_max = np.deg2rad(90.0 - np.absolute(lat_min_fit))

            # Now we need the adjusted/normalized values of the theta such that
            # full range of theta runs from 0 to pi.  At this point if you are
            # wondering why we are doing this, It would be good to refer Mike's
            # paper
            alpha = np.pi / theta_max
            theta_prime = alpha * theta
            x = np.cos(theta_prime)

            # Here we evaluate the associated legendre polynomials..from order 0
            # to order_fit we use scipy.special.lpmn() function to get the assciated
            # legendre polynomials...but it doesnt accept an array...so do loop
            # calculate the leg.pol for each value of x and append these arrays to
            # a new array
            for j,xj in enumerate(x):
                plm_temp = scipy.special.lpmn(order_fit, order_fit, xj)
                if j == 0: plm_fit = np.append([plm_temp[0]], [plm_temp[0]], axis=0)
                else: plm_fit = np.append(plm_fit, [plm_temp[0]], axis=0)

            # we need to remove the first part/subarray/element (whatever you want
            # to call it) of this array its there twice....look at j==0 part.
            plm_fit = np.delete(plm_fit, 0, 0)
            phi = np.deg2rad(mlons)

            # now do the index legender part,
            # We are doing Associated Legendre Polynomials but for each polynomial
            # we have two coefficients one for cos(phi) and the other for sin(phi),
            # so we do spherical harmonics for a real valued function using
            # sin(phi) and cos(phi) rather than exp(i*phi).

            # we use a lambda function for the index legender part, since we use
            # it in other places as well.  A good thing about python is this lambda
            # functions..u dont have to define another function for this.
            indexLgndr = lambda l, m : (m == 0 and l**2) or \
            ((l != 0) and (m != 0) and l**2 + 2 * m - 1) or 0
            kmax = indexLgndr(order_fit, order_fit)

            # set up arrays and small stuff for the eFld coeffs calculation
            theta_ecoeffs = np.zeros((kmax + 2, len(theta)))
            phi_ecoeffs = np.zeros((kmax + 2, len(theta)))

            qprime = np.array(np.where(theta_prime != 0.0))
            qprime = qprime[0]
            q = np.array(np.where(theta != 0.0))
            q = q[0]

            # finally get to converting coefficients for the potential into
            # coefficients for elec. Field
            coeff_fit_flat = coeff_fit.flatten()
            for m in range(order_fit + 1):
                for l in range(m, order_fit + 1):
                    k3 = indexLgndr(l, m)
                    k4 = indexLgndr(l, m)

                    if k3 >= 0:
                        theta_ecoeffs[k4, qprime] = theta_ecoeffs[k4, qprime] - \
                            coeff_fit_flat[k3] * alpha * l * \
                            np.cos(theta_prime[qprime]) \
                            / np.sin(theta_prime[qprime]) / self.radEarthMtrs
                        phi_ecoeffs[k4, q] = phi_ecoeffs[k4, q] - \
                            coeff_fit_flat[k3 + 1] * m / np.sin(theta[q]) / \
                            self.radEarthMtrs
                        phi_ecoeffs[k4 + 1, q] = phi_ecoeffs[k4 + 1, q] + \
                            coeff_fit_flat[k3] * m / np.sin(theta[q]) / \
                            self.radEarthMtrs

                    if l < order_fit:
                        k1 = indexLgndr(l+1, m)
                    else:
                        k1 = -1

                    k2 = indexLgndr(l, m)

                    if k1 >= 0:
                        theta_ecoeffs[k2, qprime] = theta_ecoeffs[k2, qprime] + \
                            coeff_fit_flat[k1] * alpha * (l + 1 + m) / \
                            np.sin(theta_prime[qprime]) / self.radEarthMtrs

                    if m > 0:
                        if k3 >= 0:
                            k3 = k3 + 1
                        k4 = k4 + 1

                        if k1 >= 0:
                            k1 = k1 + 1
                        k2 = k2 + 1

                        if k3 >= 0:
                            theta_ecoeffs[k4, qprime] = theta_ecoeffs[k4, qprime] \
                                - coeff_fit_flat[k3] * alpha * l * \
                                np.cos(theta_prime[qprime]) / \
                                np.sin(theta_prime[qprime]) / self.radEarthMtrs

                        if k1 >= 0:
                            theta_ecoeffs[k2, qprime] = theta_ecoeffs[k2, qprime] \
                                + coeff_fit_flat[k1] * alpha * (l + 1 + m) / \
                                np.sin(theta_prime[qprime]) / self.radEarthMtrs
            # Calculate the Elec. fld positions where
            theta_ecomp = np.zeros(theta.shape)
            phi_ecomp = np.zeros(theta.shape)

            for m in range(order_fit + 1):
                for l in range(m, order_fit + 1):
                    k = indexLgndr(l, m)
                    # Now in the IDL code we use plm_fit[:,l,m] instead of
                    # plm_fit[:,m,l] like here, this is because we have a different
                    # organization of plm_fit due to the way scipy.special.lpmn
                    # stores values in arrays...
                    if m == 0:
                        theta_ecomp = theta_ecomp + theta_ecoeffs[k,:] * \
                                     plm_fit[:,m,l]
                        phi_ecomp = phi_ecomp + phi_ecoeffs[k,:] * plm_fit[:,m,l]
                    else:
                        theta_ecomp = theta_ecomp + theta_ecoeffs[k,:] * \
                            plm_fit[:,m,l] * np.cos(m * phi) + \
                            theta_ecoeffs[k+1,:] * plm_fit[:,m,l] * np.sin(m * phi)
                        phi_ecomp = phi_ecomp + phi_ecoeffs[k,:] * \
                            plm_fit[:,m,l] * np.cos(m * phi) + \
                            phi_ecoeffs[k+1,:] * plm_fit[:,m,l] * np.sin(m * phi)

            # Store the two components of EFld into a single array
            efield_fit = np.append([theta_ecomp], [phi_ecomp], axis=0)

            # We'll calculate Bfld magnitude now, need to initialize some more
            # stuff
            alti = 300.0 * 1000.0
            b_fld_polar = -0.62e-4
            b_fld_mag = b_fld_polar * (1.0 - 3.0 * alti / self.radEarthMtrs) \
                * np.sqrt(3.0 * np.square(np.cos(theta)) + 1.0) / 2

            # get the velocity components from E-field
            vel_fit_vecs = np.zeros(efield_fit.shape)
            vel_fit_vecs[0,:] = efield_fit[1,:] / b_fld_mag
            vel_fit_vecs[1,:] = -efield_fit[0,:] / b_fld_mag

            vel_mag = np.sqrt(np.square(vel_fit_vecs[0,:]) +
                              np.square(vel_fit_vecs[1,:]))
            vel_chk_zero_inds = np.where(vel_mag != 0.0)
            vel_chk_zero_inds = vel_chk_zero_inds[0]

            vel_azm = np.zeros(vel_mag.shape)

            if len(vel_chk_zero_inds) == 0:
                vel_mag = np.array([0.0])
                vel_azm = np.array([0.0])
            else:
                if hemi == -1: vel_azm[vel_chk_zero_inds] =\
                    np.rad2deg(np.arctan2(vel_fit_vecs[1,vel_chk_zero_inds],                                                                                         vel_fit_vecs[0,vel_chk_zero_inds]))
                else: vel_azm[vel_chk_zero_inds] = np.rad2deg(np.arctan2(vel_fit_vecs[1,vel_chk_zero_inds],
                                                                         -vel_fit_vecs[0,vel_chk_zero_inds]))
        else: mlats, mlons, vel_mag, vel_azm, efield_fit = np.zeros((1))*np.nan, np.zeros((1))*np.nan,\
                                            np.zeros((1))*np.nan, np.zeros((1))*np.nan, np.zeros((2,1))*np.nan
        return mlats, mlons, vel_mag, vel_azm, efield_fit
    
    def calcCnvPots(self, rec, pot_lat_min=30.):
        """
        Calculate equipotential contour values from map data (basically
        coefficients of the fit)
        """
        stime, etime, hemi, r = rec["stime"], rec["etime"], rec["hemi"], rec["rec"]
        lat_step, lon_step = 1., 2.
        num_lats = int((90.0 - pot_lat_min) / lat_step)
        num_longs = int(360.0 / lon_step) + 1
        if "vector.mlat" in r:
            hemi_str = "north" if hemi==1 else "south"

            # get the standard location parameters.
            mlats, mlons = r["vector.mlat"], r["vector.mlon"]
            
            # Some important parameters from fitting.
            coeff_fit = np.array(r["N+2"])
            order_fit = r["fit.order"]
            lat_shft_fit = r["lat.shft"]
            lon_shft_fit = r["lon.shft"]
            lat_min_fit = r["latmin"]

            # Set up more parameters for getting the fitted vectors
            theta_max = np.deg2rad(90.0 - np.absolute(lat_min_fit))

            # we set up a grid to evaluate potential on...    
            zat_arr = np.array(range(num_lats)) * lat_step + pot_lat_min
            zat_arr = zat_arr * hemi
            zon_arr = np.array(range(num_longs))* lon_step

            # Right now create a grid kinda stuff with lats and lons
            grid_arr = np.zeros((2, num_lats * num_longs))
            grid_arr[0, :] = np.array([zat_arr.tolist()]*num_longs).ravel()
            grid_arr[1, :] = np.array([[x]*num_lats for x in zon_arr]).ravel()
            #counter1 = 0
            #for lo in zon_arr :
            #    for la in zat_arr :
            #        grid_arr[1, counter1] = lo
            #        counter1 = counter1 + 1
            #print(grid_arr[1,:180].tolist())
            #print(.tolist()[:180])
            
            # Now we need to convert a few things to spherical coordinates
            theta = np.deg2rad(90.0 - np.abs(grid_arr[0,:]))
            phi = np.deg2rad(grid_arr[1,:])

            # Now we need the adjusted/normalized values of the theta such that
            # full range of theta runs from 0 to pi.  At this point if you are
            # wondering why we are doing this, refer Mike's paper (REF NEEDED)
            alpha = np.pi / theta_max
            x = np.cos(alpha * theta)

            # Here we evaluate the associated legendre polynomials..from order 0 to
            # order_fit.  We use scipy.special.lpmn() function to get the assciated
            # legendre polynomials...but it doesn't accept an array...so do loop
            # calculate the leg.pol for each value of x and append these arrays to
            # a new array
            for j,xj in enumerate(x):
                plm_temp = scipy.special.lpmn(order_fit, order_fit, xj)

                if j == 0:
                    plm_fit = np.append([plm_temp[0]], [plm_temp[0]], axis=0)
                else:
                    plm_fit = np.append(plm_fit, [plm_temp[0]], axis=0)

            # we need to remove the first part/subarray/element (whatever you want
            # to call it) of this array. It's there twice, look at j==0 part.
            plm_fit = np.delete(plm_fit, 0, 0)

            # Get to evaluating the potential
            lmax = plm_fit.shape
            lmax = lmax[1]
            v = np.zeros(phi.shape)

            # we use a lambda function for the index legender part, since we use it
            # in other places as well.
            indexLgndr = lambda l,m : (m == 0 and l**2) or \
                ((l != 0) and (m != 0) and l**2 + 2*m - 1) or 0

            coeff_fit_flat = coeff_fit.flatten()
            for m in range(lmax):
                for l in range(m, lmax):
                    k = indexLgndr(l, m)
                    if m == 0:
                        v = v + coeff_fit_flat[k] * plm_fit[:,0,l]
                    else:
                        v = v + \
                            coeff_fit_flat[k]*np.cos(m * phi) * plm_fit[:,m,l] + \
                            coeff_fit_flat[k+1]*np.sin(m * phi) * plm_fit[:,m,l]

            pot_arr = np.zeros((num_longs, num_lats))
            pot_arr = np.reshape(v, pot_arr.shape) / 1000.0

            # lat_shft_fit and lon_shft_fit are almost always zero
            # but in case they are not... we print out a message...
            # you need an extra bit of code to account for the lat shift
            if lat_shft_fit == 0.0:
                q = np.array(np.where(np.abs(zat_arr) <= np.abs(lat_min_fit)))
                q = q[0]

                if len(q) != 0:
                    pot_arr[:,q] = 0
            else:
                estr = "LatShift is not zero, need to rewrite code for that, \
                        {:s}currently continuing assuming it is zero".format(estr)
                print(estr)

            grid_arr[1,:] = (grid_arr[1,:] + lon_shft_fit)

            lat_cntr = grid_arr[0,:].reshape((num_longs, num_lats))
            lon_cntr = grid_arr[1,:].reshape((num_longs, num_lats))
        else: lat_cntr, lon_cntr, pot_arr = np.zeros((num_longs, num_lats))*np.nan,\
            np.zeros((num_longs, num_lats))*np.nan, np.zeros((num_longs, num_lats))*np.nan
        return lat_cntr, lon_cntr, pot_arr
    
    def proc(self, rec, pot_lat_min=30., pev_params=["pot", "efield", "vel"], plots={}):
        """
        Compute E-field and Pot
        """
        print(f" Processing vel, eField, and pot for [{rec['hemi_str']}]: {rec['stime']}-to-{rec['etime']}")
        rec["N_rads"], rec["N_vecs"] = len(rec["rec"]["stid"]), np.sum(rec["rec"]["nvec"])
        if "efield" in pev_params or "vel" in pev_params: 
            mlats, mlons, vel_mag, vel_azm, efield_fit = self.calcFitCnvVel(rec)
            rec["vel_efield"] = {}
            rec["vel_efield"]["mlats"], rec["vel_efield"]["mlons"], rec["vel_efield"]["vel_mag"],\
                    rec["vel_efield"]["vel_azm"], rec["vel_efield"]["efield_fit"] = mlats, mlons, vel_mag,\
                                                                    vel_azm, efield_fit
        if "pot" in pev_params: 
            rec["pot"] = {}
            lat_cntr, lon_cntr, pot_arr = self.calcCnvPots(rec, pot_lat_min)
            rec["pot"]["lat_cntr"], rec["pot"]["lon_cntr"], rec["pot"]["pot_arr"] = lat_cntr, lon_cntr, pot_arr
        rec["coords"] = "aacgmv2_mlt"
        if (len(plots) > 0) and ("map" in plots.keys()) and plots["map"]: self.map_plot(rec, plots["map"]) 
        del rec["rec"]
        return rec
    
    def map_plot(self, rec, ftag):
        fname = ftag.format(date=rec["stime"].strftime("%Y%m%d-%H%M"), hemi=rec["hemi_str"][0].upper())
        mp = MapPlot(rec, hemi=rec["hemi_str"])
        mp.overlayHMB()
        mp.overlayCnvCntrs()
        mp.overlayMapFitVel()
        mp.set_radars()
        mp.save(fname)
        return
    
    def calcFitCnvs(self, start=None, end=None, pot_lat_min=30., cores=24, 
                    pev_params=["pot", "efield", "vel"], plots={}):
        record_list = []
        records = self.fetch_records()
        hemi = 1 if self.hemi=="north" else -1
        hemi_str = self.hemi
        for r in records:
            stime = dt.datetime(r["start.year"], r["start.month"], r["start.day"], r["start.hour"],
                                r["start.minute"], int(r["start.second"]))
            etime = dt.datetime(r["end.year"], r["end.month"], r["end.day"], r["end.hour"],
                                r["end.minute"], int(r["end.second"]))
            if (start is None) and (end is None):
                record_list.append({"stime":stime, "etime": etime, "rec":r, "hemi": hemi, "hemi_str": hemi_str})
            else:
                if (stime >= start) and (etime <= end):
                    record_list.append({"stime":stime, "etime": etime, "rec":r, "hemi": hemi, 
                                        "hemi_str": hemi_str})                
        o = []
        p0 = mp.Pool(cores)
        partial_filter = partial(self.proc, pot_lat_min=pot_lat_min, pev_params=pev_params,plots=plots)
        for rec in p0.map(partial_filter, record_list):
            o.append(rec)
        return o
    
def to_xarray(obj, pev_params, scalers, vectors, grid_params):
    """
    Convert to X-array
    """
    var = dict()
    crds = dict()
    atrs = dict()
    
    if len(pev_params) > 0: var, crds, atrs = to_xarray_pev(obj["pev_o"], pev_params, var, crds, atrs)
    if len(scalers) + len(vectors) > 0: var, crds, atrs = to_xarray_map(obj["sv_o"], scalers, 
                                                                        vectors, var, crds, atrs)
    if len(grid_params.keys()) > 0: var, crds, atrs = to_xarray_grd(obj["summ_o"], obj["reco_o"], 
                                                                    grid_params, var, crds, atrs)
    ds = xarray.Dataset(
            coords=crds,            
            data_vars=var,
            attrs=atrs,
    )
    print(ds)
    return ds

def to_xarray_map(mo, scalers, vectors, var, crds, atrs):
    crds["map.stime"] = ("map.time", mo.stime)
    crds["map.etime"] = ("map.time", mo.etime)
    for p in scalers+vectors:
        var["map."+p] = ("map.time", mo[p])
    atrs["map.desciption"] = "Processed Map data from VT SuperDARN (2021)\
            ------------------------------------------------------------\
            Parameter extension: [map]\
            ------------------------------------------------------------\
            @Powered by pydarn"
    return var, crds, atrs

def to_xarray_grd(so, ro, grid_params, var, crds, atrs):
    if ("summary" in grid_params.keys()) and (len(grid_params["summary"]) > 0):
        crds["grd.summary.stime"] = ("grd.summary.time", so.stime)
        crds["grd.summary.etime"] = ("grd.summary.time", so.etime)
        for p in grid_params["summary"]:
            var["grd.summary."+p] = ("grd.summary.time", so[p])
    if ("records" in grid_params.keys()) and (len(grid_params["records"]) > 0):
        crds["grd.records.stime"] = ("grd.records.time", ro.stime)
        crds["grd.records.etime"] = ("grd.records.time", ro.etime)
        for p in grid_params["records"]:
            var["grd.records."+p.replace("vector.","")] = ("grd.records.time", ro[p])
    atrs["grd.desciption"] = "Processed Grid data from VT SuperDARN (2021)@Powered by pydarn"            
    atrs["param.ext"] = "Parameter extension: [grd]"
    atrs["grd.summary"] = "Holds information about the data processing"
    atrs["grd.records"] = "Holds grid data records"
    return var, crds, atrs

def to_xarray_pev(o, pev_params, var, crds, atrs):
    stime, etime, hemi, N_vecs, N_rads = [], [], [], [], []
    max_ev_len = 0
    for j, i in enumerate(o):
        N_vecs.append(i["N_vecs"])
        N_rads.append(i["N_rads"])
        stime.append(i["stime"])
        etime.append(i["etime"])
        hemi.append(i["hemi_str"])
        if ("pot" in pev_params) and (j==0):  pot_arr_shape = i["pot"]["pot_arr"].shape
        if ("efield" in pev_params) or ("vel" in pev_params): 
            max_ev_len = max_ev_len if max_ev_len >= len(i["vel_efield"]["mlons"])\
                    else len(i["vel_efield"]["mlons"])
    
    stime, etime, hemi = list(set(stime)), list(set(etime)), list(set(hemi))
    stime.sort()
    etime.sort()
    
    if "pot" in pev_params:
        pot_arr, lat_cntr, lon_cntr = np.zeros((len(stime), pot_arr_shape[0], pot_arr_shape[1])), None, None
    if "vel" in pev_params or "efield" in pev_params:
        mlons, mlats = np.zeros((len(stime), max_ev_len))*np.nan, np.zeros((len(stime), max_ev_len))*np.nan
        if "vel" in pev_params: vel_mag, vel_azm = np.zeros((len(stime), max_ev_len))*np.nan,\
                    np.zeros((len(stime), max_ev_len))*np.nan
        if "efield" in pev_params: efield_fit = np.zeros((len(stime), 2, max_ev_len))
    for j, i in enumerate(o):
        if "pot" in pev_params:
            if j == 0: lat_cntr, lon_cntr = i["pot"]["lat_cntr"], i["pot"]["lon_cntr"]
            pot_arr[stime.index(i["stime"]), :, :] = i["pot"]["pot_arr"]
        if "vel" in pev_params or "efield" in pev_params:
            L = len(i["vel_efield"]["mlats"])
            mlats[stime.index(i["stime"]), :L] = i["vel_efield"]["mlats"]
            mlons[stime.index(i["stime"]), :L] = i["vel_efield"]["mlons"]
            if "vel" in pev_params: 
                vel_mag[stime.index(i["stime"]), :L] = i["vel_efield"]["vel_mag"]
                vel_azm[stime.index(i["stime"]), :L] = i["vel_efield"]["vel_azm"]
            if "efield" in pev_params:
                efield_fit[stime.index(i["stime"]), :, :L] = i["vel_efield"]["efield_fit"]
    
    crds["fparam.hemisphere"] = ("fparam.hemi", hemi)
    crds["fparam.stime"] = ("fparam.time", stime)
    crds["fparam.etime"] = ("fparam.time", etime)
    crds["fparam.n_rads"] = ("fparam.time", N_rads)
    crds["fparam.n_vecs"] = ("fparam.time", N_vecs)
    
    var["fparam.N_rads"] = (["fparam.time"], N_rads)
    var["fparam.N_vecs"] = (["fparam.time"], N_vecs)
    
    if "pot" in pev_params: 
        crds["fparam.lat_pot"] = (["fparam.pot_x","fparam.pot_y"], lat_cntr.astype(int))
        crds["fparam.lon_pot"] = (["fparam.pot_x","fparam.pot_y"], lon_cntr.astype(int))
        var["fparam.pot_arr"] = (["fparam.time", "fparam.pot_x","fparam.pot_y"], pot_arr)
    if "vel" in pev_params or "efield" in pev_params:
        crds["max_efield_vel_len"] = ("fparam.max_ev_len", range(max_ev_len))
        var["fparam.mlats"] = (["fparam.time", "fparam.max_ev_len"], mlats)
        var["fparam.mlons"] = (["fparam.time", "fparam.max_ev_len"], mlons)
        if "vel" in pev_params: 
            var["fparam.vel_mag"] = (["fparam.time", "fparam.max_ev_len"], vel_mag)
            var["fparam.vel_azm"] = (["fparam.time", "fparam.max_ev_len"], vel_azm)
        if "efield" in pev_params:
            var["fparam.efield_fit_theta"] = (["fparam.time", "fparam.max_ev_len"], efield_fit[:,0,:])
            var["fparam.efield_fit_phi"] = (["fparam.time", "fparam.max_ev_len"], efield_fit[:,1,:])
    atrs["fparam.desciption"] = "Processed %s data from VT SuperDARN (2021)@Powered by pydarn"%("-".join(pev_params))
    atrs["param.ext"] = "Parameter extension: [fparam]"
    atrs["fparam.efield_fit_theta"] = "efield north [V/m]"
    atrs["fparam.efield_fit_phi"] = "efield east [V/m]"
    atrs["fparam.vel_mag"] = "velocity magnitude [m/s]"
    atrs["fparam.vel_azm"] = "velocity azimuth [degree]"
    atrs["fparam.mlats"] = "magnetic latitudes [degrees; for fitted efields and gridded velocities]"
    atrs["fparam.mlons"] = "magnetic longitudes [degrees; for fitted efields and gridded velocities]"
    atrs["fparam.stime"] = "start time [datetime]"
    atrs["fparam.etime"] = "end time [datetime]"
    atrs["fparam.lat_pot"] = "magnetic latitudes [degrees; for fitted potentials]"
    atrs["fparam.lon_pot"] = "magnetic longitudes [degrees; for fitted potentials]"
    atrs["fparam.pot_arr"] = "fitted potential [kV]"            
    atrs["fparam.n_rads"] = "number of radar station"            
    atrs["fparam.n_vecs"] = "number of associated data vectors"            
    return var, crds, atrs
