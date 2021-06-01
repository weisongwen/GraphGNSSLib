from std_msgs.msg import String
from numpy import *
import math
import llh2ecef #ecef coordinate to llh coordinate
import ecef2llh #ecef coordinate to llh coordinate
a=6378137.0
b=6356752.314
# llh= [22.30419607114303, 114.1812441949352, 121.42342057159652]
# line_index= 156 itera_xyz= [-2418380.163540058, 5385854.5559571795, 2405656.8610494756]

def enu2ecef(originllh, enu):
    e = float(enu[0])
    n = float(enu[1])
    u = float(enu[2])
    lon = float(originllh[1]) * 3.1415926 / 180.0;
    lat = float(originllh[0]) * 3.1415926 / 180.0;

    oxyz = []
    oxyz = llh2ecef.llh2xyz(originllh)
    ox = float(oxyz[0])
    oy = float(oxyz[1])
    oz = float(oxyz[2])

    # a1 = ox - sin(lon) * e - cos(lon) * sin(lat) * n + cos(lon) * cos(lat) * u
    # a2 = oy + cos(lon) * e - sin(lon) * sin(lat) * n + cos(lat) * sin(lon) * u
    # a3 = oz + cos(lat) * n + sin(lat) * u

    oxyz[0] = ox - sin(lon) * e - cos(lon) * sin(lat) * n + cos(lon) * cos(lat) * u
    oxyz[1] = oy + cos(lon) * e - sin(lon) * sin(lat) * n + cos(lat) * sin(lon) * u
    oxyz[2] = oz + cos(lat) * n + sin(lat) * u
    return oxyz;