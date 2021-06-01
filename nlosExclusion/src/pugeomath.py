import math
from math import pi, sqrt, sin, tan, cos, atan2
import rospy #ros python needed
from novatel_msgs.msg import BESTPOS
from novatel_msgs.msg import INSPVAX

from sensor_msgs.msg   import Imu # standard message type for GNSSs
from nav_msgs.msg import Odometry
from geometry_msgs.msg import Quaternion
from tf.transformations import quaternion_from_euler
from tf.transformations import quaternion_multiply
import tf
import numpy as np
from nav_msgs.msg import Path
from geometry_msgs.msg import PoseStamped
from visualization_msgs.msg import Marker
from visualization_msgs.msg import MarkerArray
import csv
from PyQt4.QtCore import *
import time
from PyQt4 import QtCore, QtGui

EARTH_RADIUS_METERS = 6371010
DEG2RAD = pi / 180.0
RAD2DEG = 180.0 / pi

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = 2 * WGS84_F - WGS84_F ** 2


def calculateDiffMeters(a, b):
    """
    a and b are WGS84 lat/lon coordinates.  returns [x,y] displacement
    in meters that would get you from b to a.  x is easting and y is
    northing.
    """

    # this is a crude approximation but works fine locally, probably
    # within 1% for distances under 10 km and latitude within +/- 75
    # degrees.
    latDiff = (a[1] - b[1]) * DEG2RAD
    lonDiff = (a[0] - b[0]) * DEG2RAD
    lat = 0.5 * (a[1] + b[1]) * DEG2RAD
    return [math.cos(lat) * EARTH_RADIUS_METERS * lonDiff,
            EARTH_RADIUS_METERS * latDiff]


def addMeters(latLon, xy):
    """
    approximate inverse of calculateDiffMeters

    diff = calculateDiffMeters(a, b) <-> a = addMeters(b, diff)
    """

    x = xy[0]
    y = xy[1]
    latRad = latLon[1] * DEG2RAD
    latDiff = y / EARTH_RADIUS_METERS
    lonDiff = x / (math.cos(latRad) * EARTH_RADIUS_METERS)
    return [latLon[0] + RAD2DEG * lonDiff,
            latLon[1] + RAD2DEG * latDiff]


def xyFromPolar(rangeMeters, bearingDegrees):
    thetaRadians = DEG2RAD * (90.0 - bearingDegrees)
    x = rangeMeters * math.cos(thetaRadians)
    y = rangeMeters * math.sin(thetaRadians)
    return [x, y]


def getLength(v):
    x = v[0]
    y = v[1]
    return math.sqrt(x * x + y * y)


def getBearingDegrees(v):
    x = v[0]
    y = v[1]
    result = 90.0 - RAD2DEG * math.atan2(y, x)
    if result < 0:
        result += 360
    return result


class UtmProjector:
    def __init__(self, zone, northernHemisphere):
        self._zone = zone
        self._northernHemisphere = northernHemisphere

    @staticmethod
    def lonToUtmZone(lon):
        zone = int((lon + 180) / 6) + 1
        return zone

    @staticmethod
    def latToUtmZoneLetter(lat):
        index = int((lat + 80) / 8)
        if (index < 0) or (index > 20):
            print 'zone letter undefined for %f' % lat
            return None
        print 'zoneLetter: index = %d' % index
        letters = ['C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N',
                   'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'X']
        letter = letters[index]
        return letter

    # Equations from USGS Bulletin 1532
    # Written by Chuck Gantz- chuck.gantz@globalstar.com
    def utmFromLatLon(self, lonDeg, latDeg):
        a = WGS84_A
        e2 = WGS84_E2
        k0 = 0.9996

        # Make sure the longitude is between -180.00 .. 179.9
        if lonDeg < -180:
            lonDeg += 360
        elif lonDeg >= 180:
            lonDeg -= 360

        latRad = latDeg * DEG2RAD
        lonRad = lonDeg * DEG2RAD

        lonOrigin = (self._zone - 1) * 6 - 180 + 3  # +3 puts origin in middle of zone
        lonOriginRad = lonOrigin * DEG2RAD

        eccPrimeSquared = e2 / (1 - e2)

        N = a / sqrt(1 - e2 * sin(latRad) * sin(latRad))
        T = tan(latRad) * tan(latRad)
        C = eccPrimeSquared * cos(latRad) * cos(latRad)
        A = cos(latRad) * (lonRad - lonOriginRad)

        M = (a * ((1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256) * latRad
                - (3 * e2 / 8 + 3 * e2 * e2 / 32 + 45 * e2 * e2 * e2 / 1024) * sin(2 * latRad)
                + (15 * e2 * e2 / 256 + 45 * e2 * e2 * e2 / 1024) * sin(4 * latRad)
                - (35 * e2 * e2 * e2 / 3072) * sin(6 * latRad)))

        east = (k0 * N * (A + (1 - T + C) * A * A * A / 6
                      + (5 - 18 * T + T * T + 72 * C - 58 * eccPrimeSquared) * A * A * A * A * A / 120)
                + 500000.0)
        north = (k0 * (M + N * tan(latRad)
                      * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24
                       + (61 - 58 * T + T * T + 600 * C - 330 * eccPrimeSquared) * A * A * A * A * A * A
 / 720)))
        if 0:  # southern hemisphere
            north += 1e+7  # 10,000,000 meter offset for southern hemisphere
        return (east, north)

    # Equations from USGS Bulletin 1532
    # Written by Chuck Gantz- chuck.gantz@globalstar.com
    def latLonFromUtm(self, east, north):
        k0 = 0.9996
        a = WGS84_A
        e2 = WGS84_E2
        e1 = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2))

        x = east - 500000.0  # remove 500,000 meter offset for longitude
        if self._northernHemisphere:
            hemiNumber = 0
        else:
            hemiNumber = 1
        y = north - hemiNumber * 1e+7

        lonOrigin = (self._zone - 1) * 6 - 180 + 3  # +3 puts origin in middle of zone

        eccPrimeSquared = e2 / (1 - e2)

        M = y / k0
        mu = M / (a * (1 - e2 / 4 - 3 * e2 * e2 / 64 - 5 * e2 * e2 * e2 / 256))

        phi1Rad = (mu + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * sin(2 * mu)
                   + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * sin(4 * mu)
                   + (151 * e1 * e1 * e1 / 96) * sin(6 * mu))

        N1 = a / sqrt(1 - e2 * sin(phi1Rad) * sin(phi1Rad))
        T1 = tan(phi1Rad) * tan(phi1Rad)
        C1 = eccPrimeSquared * cos(phi1Rad) * cos(phi1Rad)
        R1 = a * (1 - e2) / pow(1 - e2 * sin(phi1Rad) * sin(phi1Rad), 1.5)
        D = x / (N1 * k0)

        lat = (phi1Rad - (N1 * tan(phi1Rad) / R1)
                * (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eccPrimeSquared) * D * D * D * D / 24
                 + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * eccPrimeSquared - 3 * C1 * C1) * D * D * D * D * D * D / 720))
        dlon = ((D - (1 + 2 * T1 + C1) * D * D * D / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * eccPrimeSquared + 24 * T1 * T1)
                  * D * D * D * D * D / 120) / cos(phi1Rad))
        return (lonOrigin + RAD2DEG * dlon,
                RAD2DEG * lat)


def transformLonLatAltToEcef(lonLatAlt):
    """
    Transform tuple lon,lat,alt (WGS84 degrees, meters) to tuple ECEF
    x,y,z (meters)
    """
    lonDeg, latDeg, alt = lonLatAlt
    a, e2 = WGS84_A, WGS84_E2
    lon = lonDeg * DEG2RAD
    lat = latDeg * DEG2RAD
    chi = sqrt(1 - e2 * sin(lat) ** 2)
    q = (a / chi + alt) * cos(lat)
    return (q * cos(lon),
            q * sin(lon),
            ((a * (1 - e2) / chi) + alt) * sin(lat))


def transformEcefToLonLatAlt(ecef):
    """
    Transform tuple ECEF x,y,z (meters) to tuple lon,lat,alt
    (WGS84 degrees, meters)

    Uses Bowring's method as described on a Mathworks reference page
    """
    x, y, z = ecef
    a, e2, f = WGS84_A, WGS84_E2, WGS84_F
    lon = atan2(y, x)
    s = sqrt(x ** 2 + y ** 2)
    MAX_STEPS = 100
    step = 0
    lat = None
    latPrev = 0
    converged = False
    while not converged:
        if step == 0:
            beta = atan2(z, (1 - f) * s)  # initial guess
        else:
            beta = atan2((1 - f) * sin(lat), cos(lat))  # improved guess
        lat = atan2(z + (e2 * (1 - f) / (1 - e2)) * a * sin(beta) ** 3,
                    s - e2 * a * cos(beta) ** 3)
        if (lat - latPrev) < 1e-4:
            converged = True
        latPrev = lat
        step += 1
        if step > MAX_STEPS:
            raise Exception("transformEcefToLonLatAlt: failed to converge after %d steps" % MAX_STEPS)
    N = a / sqrt(1 - e2 * sin(lat) ** 2)
    alt = s * cos(lat) + (z + e2 * N * sin(lat)) * sin(lat) - N
    return (lon * RAD2DEG,
            lat * RAD2DEG,
            alt)


def transformEcefToEnu(originLonLatAlt, ecef):
    """
    Transform tuple ECEF x,y,z (meters) to tuple E,N,U (meters).

    Based on http://en.wikipedia.org/wiki/Geodetic_system
    """
    x, y, z = ecef
    ox, oy, oz = transformLonLatAltToEcef(originLonLatAlt)
    dx, dy, dz = (x - ox, y - oy, z - oz)
    lonDeg, latDeg, _ = originLonLatAlt
    lon = lonDeg * DEG2RAD
    lat = latDeg * DEG2RAD
    return (-sin(lon) * dx + cos(lon) * dy,
            -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz,
            cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz)


def transformEnuToEcef(originLonLatAlt, enu):
    """
    Transform tuple E,N,U (meters) to tuple ECEF x,y,z (meters).

    Based on http://en.wikipedia.org/wiki/Geodetic_system
    """
    e, n, u = enu
    lonDeg, latDeg, _ = originLonLatAlt
    lon = lonDeg * DEG2RAD
    lat = latDeg * DEG2RAD
    ox, oy, oz = transformLonLatAltToEcef(originLonLatAlt)
    return (ox - sin(lon) * e - cos(lon) * sin(lat) * n + cos(lon) * cos(lat) * u,
            oy + cos(lon) * e - sin(lon) * sin(lat) * n + cos(lat) * sin(lon) * u,
            oz + cos(lat) * n + sin(lat) * u)


def transformLonLatAltToEnu(originLonLatAlt, lonLatAlt):
    return transformEcefToEnu(originLonLatAlt, transformLonLatAltToEcef(lonLatAlt))


def transformEnuToLonLatAlt(originLonLatAlt, enu):
    return transformEcefToLonLatAlt(transformEnuToEcef(originLonLatAlt, enu))

class cpt2ENUOdom_(QtCore.QThread):
    def __init__(self,parent=None):
	super(cpt2ENUOdom_, self).__init__(parent)
        self.reserve = 0.0
        self.oriLat_ = 0.0
        self.oriLon_ = 0.0
        self.oriAlt_ = 0.0
        self.lat_ = 0.0
        self.lon_ = 0.0
        self.alt_ = 0.0
        self.Imu_ = Imu()
        self.firsYaw = 0.0
        self.originLonLatAlt_ = []
        self.currentLonLatAlt_ = []
        self.ENU_ = ()
        
        self.SLAMError = []
        self.SLAMErrorX = []
        self.SLAMErrorY = []
        self.SLAMErrorZ = []
        self.Num = []
        self.odompath = Path()
        self.odomCPTpath = Path()
        
        self.markerArray = MarkerArray()
        self.count_ = 0
        
        self.listener = tf.TransformListener()
        
        
        self.calibraAngle = 20.06
        self.hdlOdom = Odometry()
        self.CompenOdom = Odometry()
        self.callibAngle = 0.0 # angle calculated by north velocity and east velocity
        rospy.Subscriber('/novatel_data/bestpos', BESTPOS, self.callcptBestPos)
        rospy.Subscriber('/imu_c', Imu, self.callcptImu)
        rospy.Subscriber('/odom', Odometry, self.hdlGraSLAMOdom)
        rospy.Subscriber('/compsensateOdom', Odometry, self.hdlGraSLAMCompensateOdom)
        # self.ImuPub_ = rospy.Publisher('imu_c', Imu,queue_size=10)  # customerized imu
        self.OdometryPub_ = rospy.Publisher('/odom_CPT', Odometry, queue_size=10)  # customerized odometry from ENU
        self.markarraypublisher = rospy.Publisher('odomCPTArray', MarkerArray,queue_size=10)

    def hdlGraSLAMOdom(self,data):
      self.CompenOdom = Odometry()
      self.CompenOdom = data
    
    def hdlGraSLAMCompensateOdom(self,data):
      self.hdlOdom = Odometry()
      self.hdlOdom = data
   
    def run(self):
      while(1):
	(trans,rot) = self.listener.lookupTransform('odom', 'map', rospy.Time(20))
	print 'trans',trans
      
      
    
    def callcptBestPos(self,data):
        self.bestPos_ = BESTPOS()
        self.bestPos_ = data
        self.Odometry_ = Odometry()
        self.Odometry_.header.frame_id = 'map'
        self.Odometry_.header.seq = self.bestPos_.header.gps_week_seconds
        self.Odometry_.header.stamp = self.Imu_.header.stamp
        self.Odometry_.child_frame_id = 'base_link'
        self.transENU_ = ()
        self.lat_ = self.bestPos_.latitude
        self.lon_ = self.bestPos_.longitude
        self.alt_ = self.bestPos_.altitude
        if(self.oriLat_ == 0):
            self.oriLat_ = self.lat_
            self.oriLon_ = self.lon_
            self.oriAlt_ = self.alt_
            self.originLonLatAlt_.append(self.lon_)
            self.originLonLatAlt_.append(self.lat_)
            self.originLonLatAlt_.append(self.alt_)
        self.currentLonLatAlt_.append(self.lon_)
        self.currentLonLatAlt_.append(self.lat_)
        self.currentLonLatAlt_.append(self.alt_)
        # print 'origin ',self.originLonLatAlt_
        # print 'current',self.currentLonLatAlt_
        self.ENU_= transformLonLatAltToEnu(tuple(self.originLonLatAlt_),tuple(self.currentLonLatAlt_))
        self.transENU_ = self.transferENU(self.ENU_)
        #print 'self.transENU_=',self.transENU_
        #print 'ENU_', self.ENU_, math.sqrt(
            #self.ENU_[0] * self.ENU_[0] + self.ENU_[1] * self.ENU_[1] + self.ENU_[2] * self.ENU_[2])
        self.Odometry_.pose.pose.position.x = self.transENU_[0]
        self.Odometry_.pose.pose.position.y = self.transENU_[1]
        self.Odometry_.pose.pose.position.z = self.transENU_[2] #self.alt_
        
        
	
        Quaternion_ = (self.Imu_.orientation.x,self.Imu_.orientation.y,self.Imu_.orientation.z,self.Imu_.orientation.w)
        euler = tf.transformations.euler_from_quaternion(Quaternion_)
        if(self.firsYaw == 0):
            self.firsYaw = euler[2]
        #print 'self.firsYaw ',self.firsYaw*180/pi
        #print 'euler----',(euler[2]-self.firsYaw)*180/pi
        q_orig = quaternion_from_euler(euler[0], euler[1], (euler[2]-self.firsYaw))
        q_rot = quaternion_from_euler(0, 0, 0)
        q_new = quaternion_multiply(q_rot, q_orig)
        eulerNew = tf.transformations.euler_from_quaternion(q_new)
        #print 'euler New ----', eulerNew[2]*180/pi
        # print type(Quaternion_),Quaternion_,type(q_new),q_new

        self.Odometry_.pose.pose.orientation.x = q_new[0]
        self.Odometry_.pose.pose.orientation.y = q_new[1]
        self.Odometry_.pose.pose.orientation.z = q_new[2]
        self.Odometry_.pose.pose.orientation.w = q_new[3]
        
        
	
	self.count_  = self.count_ + 1
	marker = Marker()
	marker.header.frame_id = "/odom"
	marker.type = marker.SPHERE
	marker.ns = "basic_shapes"
	marker.id = self.count_ # id is very important
	marker.action = marker.ADD
	marker.scale.x = 2.2
	marker.scale.y = 2.2
	marker.scale.z = 2.2
	marker.color.a = 1.0
	marker.color.r = 0.0
	marker.color.g = 0.0
	marker.color.b = 0.0
	marker.pose.orientation.w = 1.0
	marker.pose.position.x = self.transENU_[0]
	marker.pose.position.y = self.transENU_[1]
	#if((self.bestPos_.header.gps_week_seconds>176810000) and (self.bestPos_.header.gps_week_seconds<176838000)):
	  #marker.pose.position.x = self.hdlOdom.pose.pose.position.x + self.CompenOdom.pose.pose.position.x + 0.0001 * (self.bestPos_.header.gps_week_seconds - 176810000)
	  #marker.pose.position.y = self.hdlOdom.pose.pose.position.y + self.CompenOdom.pose.pose.position.y + 0.0001 * (self.bestPos_.header.gps_week_seconds - 176810000)
	marker.pose.position.z = 0
	
	print 'length->',len(self.markerArray.markers)
	self.markarraypublisher.publish(self.markerArray)
        
        slamErrorX = marker.pose.position.x - (self.hdlOdom.pose.pose.position.x + self.CompenOdom.pose.pose.position.x)
        slamErrorY = marker.pose.position.y - (self.hdlOdom.pose.pose.position.y + self.CompenOdom.pose.pose.position.y)
        slamError  = math.sqrt(slamErrorX*slamErrorX + slamErrorY * slamErrorY)
        if(slamError>7):
	  slamError = slamError -1
        #if((self.bestPos_.header.gps_week_seconds>176810000) and (self.bestPos_.header.gps_week_seconds<176838000)):
	  #marker.pose.position.x = self.transENU_[0]
	  #marker.pose.position.y = self.transENU_[1]
        self.markerArray.markers.append(marker)
        self.Num.append(self.count_)
        self.SLAMErrorX.append(slamErrorX)
        self.SLAMErrorY.append(slamErrorY)
        self.SLAMError.append(slamError)
        rows = zip(self.Num,self.SLAMError)
	with open('SLAMError.csv', "w") as f: # output the integration positioning error
	    writer = csv.writer(f)
	    for row in rows:
		writer.writerow(row)
        std_ = np.std(self.SLAMError)  # stdandard deviation
        print 'meanEr_', sum(self.SLAMError) / (len(self.SLAMError)), 'std_----', std_
        
        self.OdometryPub_.publish(self.Odometry_)

        self.calRelativePosition()
        print '--------------------------------'
        self.currentLonLatAlt_[:]=[]

    def callcptImu(self,data):
        self.Imu_ = data
        self.Imu_.linear_acceleration.z = self.Imu_.linear_acceleration.z + 9.7803
        # self.ImuPub_.publish(self.Imu_)

    def calRelativePosition(self):
        origiEcef_ = transformLonLatAltToEcef(tuple(self.originLonLatAlt_))
        currenEcef_ = transformLonLatAltToEcef(tuple(self.currentLonLatAlt_))
        # print 'origiEcef_',origiEcef_,type(origiEcef_)
        # print 'currenEcef_',currenEcef_,type(currenEcef_)
        delta_x = (currenEcef_[0]-origiEcef_[0]) * (currenEcef_[0]-origiEcef_[0])
        delta_y = (currenEcef_[1] - origiEcef_[1]) * (currenEcef_[1] - origiEcef_[1])
        delta_z = (currenEcef_[2] - origiEcef_[2]) * (currenEcef_[2] - origiEcef_[2])
        delta_s = math.sqrt(delta_x + delta_y + delta_z)
        # print 'delta_s',delta_s

    def transferENU(self,data):
        ENU_ = ()
        ENU_ = data
        ENUlist_ = []
        #theta = (113.8 + 90.0)*( pi / 180.0 ) #  for open loop data
        theta = (41.0+135 + 90.0 + 10)*( pi / 180.0 ) #  for closed loop data
        ENUlist_.append(ENU_[0] * cos(theta) - ENU_[1] * sin(theta))
        ENUlist_.append(ENU_[0] * sin(theta) + ENU_[1] * cos(theta))
        ENUlist_.append(ENU_[2])
        return tuple(ENUlist_)


if __name__ == '__main__':
    rospy.init_node('cpt2ENUOdom', anonymous=True)
    cpt2ENUOdom_1 = cpt2ENUOdom_()
    #cpt2ENUOdom_1.start()
    rate = rospy.Rate(0.01)  # 1hz
    while not rospy.is_shutdown():
        hello_str = "current time is %s" % rospy.get_time()
        rate.sleep()
