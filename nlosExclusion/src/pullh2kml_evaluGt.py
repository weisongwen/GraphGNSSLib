#! /usr/bin/env python
# -*- coding=utf-8 -*-
# finished by Dr. WEN
"""
    Function: subscribe GNSS related rostopic and save it into a KML file
    Welson Wen, Ph.D.
    https://sites.google.com/view/weisongwen/about-me

    subcribe : '/novatel_data/bestpos' latitude longitude

"""
from lxml import etree  #将KML节点输出为字符串
# import xlrd             #操作Excel
import rospy
from pykml.factory import KML_ElementMaker as KML #使用factory模块
import csv # csv reading needed library
# import pandas as pd
from novatel_msgs.msg import BESTPOS

class pullh2kml_eval():
    def __init__(self):
        rospy.Subscriber('/novatel_data/bestpos', BESTPOS, self.callcptBestPos_llh)
        self.lat_ = [] # used to save latitude
        self.lon_ = [] # used to save longitude
        self.GPS_Week_Second = 0.0
        self.writeToKML = 0.0

    def callcptBestPos_llh(self,data):
        self.bestPos_ = BESTPOS()
        self.bestPos_ = data
        self.lat_.append(float(self.bestPos_.latitude))
        print 'len(self.lat_)',len(self.lat_)
        self.lon_.append(float(self.bestPos_.longitude))
        self.GPS_Week_Second = self.bestPos_.header.gps_week_seconds
        print 'GPS_Week_Second',self.GPS_Week_Second


if __name__ == '__main__':
    rospy.init_node('pullh2kml_evaluGt', anonymous=True)
    pullh2kml_eval_ =pullh2kml_eval()
    rate = rospy.Rate(0.002)#
    preTim = 0.0
    while not rospy.is_shutdown():
        #rate.sleep()
        print 'GPS Time ',preTim,pullh2kml_eval_.GPS_Week_Second
        #if( (preTim > 0) and (preTim == pullh2kml_eval_.GPS_Week_Second) and (pullh2kml_eval_.writeToKML ==0)):
	if( len(pullh2kml_eval_.lon_)>5):
            print 'write llh to kml'
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
            with open('/home/wws/CV_GNSS/Gt.kml', 'w') as fp:
                fp.write(content)

        preTim = pullh2kml_eval_.GPS_Week_Second