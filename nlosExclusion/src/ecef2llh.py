from std_msgs.msg import String
from numpy import *
import math
def xyz2llh(xyz):
    # %XYZ2LLH	Convert from ECEF cartesian coordinates to
    # %               latitude, longitude and height.  WGS-84
    # %
    # %	llh = XYZ2LLH(xyz)
    # %
    # %    INPUTS
    # %	xyz(1) = ECEF x-coordinate in meters
    # %	xyz(2) = ECEF y-coordinate in meters
    # %	xyz(3) = ECEF z-coordinate in meters
    # %
    # %    OUTPUTS
    # %	llh(1) = latitude in radians
    # %	llh(2) = longitude in radians
    # %	llh(3) = height above ellipsoid in meters
    #
    # %	Reference: Understanding GPS: Principles and Applications,
    # %	           Elliott D. Kaplan, Editor, Artech House Publishers,
    # %	           Boston, 1996.
    # %
    # %	M. & S. Braasch 10-96
    # %	Copyright (c) 1996 by GPSoft
    # %	All Rights Reserved.
    # %
    x = float(xyz[0])
    y = float(xyz[1])
    z = float(xyz[2])
    x2 = x**2
    y2 = y**2
    z2 = z**2

    a = 6378137.0000 # earth radius in meters
    b = 6356752.3142 # earth semiminor in meters
    e = sqrt (1-(b/a)**2)
    b2 = b*b
    e2 = e**2
    ep = e*(a/b)
    r = sqrt(x2+y2)
    r2 = r*r
    E2 = a**2 - b**2
    F = 54*b2*z2
    G = r2 + (1-e2)*z2 - e2*E2
    c = (e2*e2*F*r2)/(G*G*G)
    s = ( 1 + c + sqrt(c*c + 2*c) )**(1/3)
    P = F / (3 * ((s+1/s+1)**2) * G*G)
    Q = sqrt(1+2*e2*e2*P)
    ro = -(P*e2*r)/(1+Q) + sqrt((a*a/2)*(1+1/Q) - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2)
    tmp = (r - e2*ro)**2
    U = sqrt( tmp + z2 )
    V = sqrt( tmp + (1-e2)*z2 )
    zo = (b2*z)/(a*V)

    height = U*( 1 - b2/(a*V) )

    lat = math.atan((z + ep*ep*zo)/r)

    temp = math.atan(y/x)
    if x >=0:
        long = temp
    elif (x < 0) and (y >= 0):
        long = pi + temp
    else:
        long = temp - pi
    llh=[]
    llh.append(float(lat)*(180/pi))
    llh.append(float(long)*(180/pi))
    llh.append(float(height))
    # print 'lngth of llh',len(llh)
    return llh