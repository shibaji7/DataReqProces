import pandas as pd
import os
import datetime as dt

import glob
import traceback
import pydarn
import bz2

fetcher = "data/Sur,Dibyendu"
hemis = ["north", "south"]
if not os.path.exists(fetcher): os.system("mkdir -p " + fetcher)

def read_map_files(dates, hemi):
    for d in dates:
        year, month, day = d.year, d.month, d.day
        fname = ("/sd-data/{:d}/map2/%s/{:d}{:02d}{:02d}.%s.map2.bz2".format(year, year, month, day))%(hemi, hemi)
        if os.path.exists(fname):
            print("Fname - ", fname)
            try:
                with bz2.open(fname) as fp: fs = fp.read()
                reader = pydarn.SuperDARNRead(fs, True)
                _md = reader.read_map()
                cpcps = [i["pot.drop"] for i in _md]
                times = [d + dt.timedelta(hours=int(i["start.hour"])) + dt.timedelta(minutes=int(i["start.minute"])) + 
                        dt.timedelta(seconds=int(i["start.second"])) for i in _md]
                x = pd.DataFrame()
                x["cpcps"], x["times"] = cpcps, times
                x.to_csv(fetcher + ("/{:d}{:02d}{:02d}_%s.csv".format(year, month, day))%hemi, header=True, index=False)
            except: traceback.print_exc()
        else: print("File not found for - ", d, hemi)
    return

def create_hourly_averages(dates=[dt.datetime(2000,1,1), dt.datetime(2020,12,31)]):
    start, end = dates[0], dates[1]
    for hemi in hemis:
        dn = start
        while dn <= end:
            year, month, day = dn.year, dn.month, dn.day
            fname = fetcher + ("/{:d}{:02d}{:02d}_%s.csv".format(year, month, day))%hemi
            if not os.path.exists(fname): 
                read_map_files([dn], hemi)
                if os.path.exists(fname):
                    df = pd.read_csv(fname, parse_dates=["times"])
                    hour = pd.to_timedelta(df["times"].dt.hour, unit="H") 
                    df = df.groupby(hour).mean().reset_index()
                    df["times"] = [dn + u for u in df.times]
                    df.to_csv(fname, header=True, index=False)
            dn += dt.timedelta(days=1)
        files = glob.glob(fetcher + "/*_%s.csv"%hemi)
        X = pd.DataFrame()
        files.sort()
        print("Number of files - ", len(files))
        for f in files:
            X = pd.concat([X, pd.read_csv(f)])
        X.to_csv(fetcher + "/%s.csv"%hemi, index=False, header=True)
    return

if __name__ == "__main__":
    indv = False
    if indv:
        dates = [dt.datetime(2004,7,22), dt.datetime(2012,3,9), dt.datetime(2013,3,17), dt.datetime(2015,3,17),
                dt.datetime(2017,9,7), dt.datetime(2017,9,8), dt.datetime(2017,9,9), dt.datetime(2015,6,22), dt.datetime(2015,6,23)]
        for hemi in hemis:
            read_map_files(dates, hemi)
        os.system("zip -r %s.zip %s/*"%(fetcher, fetcher))
    else: create_hourly_averages()
    os.system("rm -rf *.log")
    pass
