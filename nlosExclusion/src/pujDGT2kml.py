#! /usr/bin/env python
# -*- coding=utf-8 -*-
# finished by Dr. WEN
"""
    Function: read lat lon height from jingdong standard dataset
    Welson Wen, Ph.D.
    https://sites.google.com/view/weisongwen/about-me

    subcribe : '/ublox_gps_node/fix' latitude longitude

"""
import rospy
import rospy
from pykml.factory import KML_ElementMaker as KML 
import csv # csv reading needed library
from lxml import etree 

from sensor_msgs.msg   import NavSatFix # standard message type for GNSSs\
import llh2ecef #ecef coordinate to llh coordinate
import ecef2llh #ecef coordinate to llh coordinate
import enu2ecef

import math
from numpy import *
import math


class pullh2kml_eval():
    def __init__(self):
        rospy.Subscriber('/fix2', NavSatFix, self.callublox_llh)
        self.oriLLH = []
        self.oriLLH.append(22.3111737354)
        self.oriLLH.append(114.16931259)
        self.oriLLH.append(10)
        self.ENU = []
        self.lat_ = [] # used to save latitude
        self.lon_ = [] # used to save longitude
        self.GPS_Week_Second = 0.0
        self.writeToKML = 0.0
        self.standardGT2kml()
        

    def callublox_llh(self,data):
        self.navsatfix_ = NavSatFix()
        self.navsatfix_ = data
        self.lat_.append(float(self.navsatfix_.latitude))
        print 'len(self.lat_)',len(self.lat_)
        self.lon_.append(float(self.navsatfix_.longitude))
        self.GPS_Week_Second = self.navsatfix_.header.stamp
        print 'GPS_Week_Second',self.GPS_Week_Second

    def standardGT2kml(self):
        print 'begin read standard data'
        self.Fcsv_GNSS = csv.reader(open('/home/wws/map_trajectory.csv','r'))  # read csv context to csv_reader variable
        for rowCsv in self.Fcsv_GNSS:
            # self.lat_.append(float(rowCsv[1]))
            # self.lon_.append(float(rowCsv[2]))
            east = rowCsv[1]
            north = rowCsv[2]
            up = rowCsv[3]
            prex_ = float(east)
            prey_ = float(north)

            origin_azimuth = 348.747632943 + 180 - 1.35 # 1.5
            theta = -1 * (origin_azimuth - 90)*( 3.141592 / 180.0 )
            # east = (prex_ * cos(theta) - prey_ * sin(theta)) - 1.3
            # north = (prex_ * sin(theta) + prey_ * cos(theta)) + 1.3
            east = (prex_ * cos(theta) - prey_ * sin(theta)) + 1.6 # best suit for data in moko east 
            north = (prex_ * sin(theta) + prey_ * cos(theta)) + 1.3 # best suit moko east


            self.ENU = []
            self.ENU.append(east)
            self.ENU.append(north)
            self.ENU.append(up)
            self.ecef = enu2ecef.enu2ecef(self.oriLLH, self.ENU)
            print 'self.ENU-> ',self.ENU,'\n'
            print 'self.ecef-> ',self.ecef,'\n'
            self.llh_cal = []
            self.llh_cal = ecef2llh.xyz2llh(self.ecef)
            self.lat_.append(float(self.llh_cal[0]))
            self.lon_.append(float(self.llh_cal[1]))
            print 'self.llh_cal-> ',self.llh_cal,'\n'
        print 'finish tranversal standard data',len(self.lat_)

if __name__ == '__main__':
    rospy.init_node('pullh2kml_evaluublox', anonymous=True)
    pullh2kml_eval_ =pullh2kml_eval()
    rate = rospy.Rate(0.002)#
    preTim = 0.0
    while not rospy.is_shutdown():
	if( len(pullh2kml_eval_.lon_)>5):
            pullh2kml_eval_.writeToKML = 1
            # 使用第一个点创建Folder
            fold = KML.Folder(KML.Placemark(
                KML.Point(KML.coordinates(str(pullh2kml_eval_.lon_[0]) + ',' + str(pullh2kml_eval_.lat_[0]) + ',0'))
            )
            )
            # 将剩余的点追加到Folder中
            for i in range(1, len(pullh2kml_eval_.lon_)):
                fold.append(KML.Placemark(
                    KML.Point(KML.coordinates(str(pullh2kml_eval_.lon_[i]) + ',' + str(pullh2kml_eval_.lat_[i]) + ',0')))
                )
            # 使用etree将KML节点输出为字符串数据
            content = etree.tostring(etree.ElementTree(fold), pretty_print=True)
            # 保存到文件，然后就可以在Google地球中打开了
            with open('/home/wws/map_trajectory.kml', 'w') as fp:
                fp.write(content)

        preTim = pullh2kml_eval_.GPS_Week_Second
