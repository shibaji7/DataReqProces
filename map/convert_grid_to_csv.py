import pandas as pd
import os
import datetime as dt

import traceback
import pydarn
import bz2

fetcher = "data/Starr,Gregory"
hemis = ["north", "south"]
if not os.path.exists(fetcher): os.system("mkdir -p " + fetcher)
comments = """# Key-Value description:
              # ----------------------
              # vector_mlat: Magnetic Latitude
              # vector_mlon: Magnetic Longitude
              # vector_kvect: Magnetic Azimuth 
              # vector_stid: Station identifier 
              # vector_channel: Channel number 
              # vector_index: Grid cell index
              # vector_vel_median: Weighted mean velocity magnitude
              # vector_vel_sd: Velocity standard deviation
              # vector_pwr_median: Weighted mean power
              # vector_pwr_sd: Power standard deviation
              # vector_wdt_median: Weighted mean spectral width
              # vector_wdt_sd: Standard deviation of spectral width\n"""

def read_map_files(dates, hemi):
    for d in dates:
        year, month, day = d.year, d.month, d.day
        fname = ("/sd-data/{:d}/map2/%s/{:d}{:02d}{:02d}.%s.map2.bz2".format(year, year, month, day))%(hemi, hemi)
        tofname = fetcher + ("/{:d}{:02d}{:02d}_%s.csv".format(year, month, day))%hemi
        if os.path.exists(fname) and not os.path.exists(tofname):
            print("Fname - ", fname)
            try:
                with bz2.open(fname) as fp: fs = fp.read()
                reader = pydarn.SuperDARNRead(fs, True)
                _md = reader.read_map()
                objs = []
                keys = ["vector.mlat", "vector.mlon", "vector.kvect", "vector.stid", "vector.channel", "vector.index",
                        "vector.vel.median", "vector.vel.sd", "vector.pwr.median", "vector.pwr.sd", "vector.wdt.median", 
                        "vector.wdt.sd"]
                for i in _md:
                    if "vector.mlat" in i.keys():
                        num = len(i["vector.mlat"])
                        tm = d + dt.timedelta(hours=int(i["start.hour"])) + dt.timedelta(minutes=int(i["start.minute"])) +\
                            dt.timedelta(seconds=int(i["start.second"]))
                        for _n in range(num):
                            obj = {"time": tm}
                            for k in keys:
                                if k in i.keys(): obj[k.replace(".", "_")] = i[k][_n]
                                else: obj[k.replace(".", "_")] = np.nan
                            objs.append(obj)
                x = pd.DataFrame.from_records(objs)
                #fname = fetcher + ("/{:d}{:02d}{:02d}_%s.csv".format(year, month, day))%hemi
                with open(tofname, "w") as f:
                    f.write(comments.replace(" ", ""))
                    x.to_csv(f, header=True, index=False, float_format="%.2f")
            except: traceback.print_exc()
        else: print("File not found for - ", fname)
    return

if __name__ == "__main__":
    ok_zip = False
    Is = (dt.datetime(2020,12,31)-dt.datetime(2014,12,31)).days + 1
    dates = [dt.datetime(2014,12,31) + dt.timedelta(days=i) for i in range(Is)]
    for hemi in hemis:
        read_map_files(dates, hemi)
    os.system("rm -rf *.log")
    if ok_zip: os.system("zip -r %s.zip %s"%(fetcher, fetcher))
    pass
