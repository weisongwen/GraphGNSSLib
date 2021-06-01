from std_msgs.msg import String
from numpy import *
import math
a=6378137.0
b=6356752.314
# llh= [22.30419607114303, 114.1812441949352, 121.42342057159652]
# line_index= 156 itera_xyz= [-2418380.163540058, 5385854.5559571795, 2405656.8610494756]
def llh2xyz(llh):
    xyz=[]
    lat = float(llh[0]) * 3.1415926 / 180.0
    lon = float(llh[1]) * 3.1415926 / 180.0
    hr = float(llh[2])
    n = a * a / sqrt(a * a * cos(lat) * cos(lat) + b * b * sin(lat) * sin(lat))
    Rx = (n + hr) * cos(lat) * cos(lon)
    Ry = (n + hr) * cos(lat) * sin(lon)
    Rz = (b * b / (a * a) * n + hr) * sin(lat)
    xyz.append(float(Rx))
    xyz.append(float(Ry))
    xyz.append(float(Rz))
    return xyz