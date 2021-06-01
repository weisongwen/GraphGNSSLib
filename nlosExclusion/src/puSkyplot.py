#!/usr/bin/env python
# license removed for brevity
"""
    Nlos exclusion show
    Welson Wen, Ph.D.
    https://sites.google.com/view/weisongwen/about-me
"""
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as  NavigationToolbar
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import rospy
from sensor_msgs.msg import LaserScan
from sensor_msgs.msg import PointCloud2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd # pandas to pd
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import sys
import  math
from matplotlib.patches import Circle
import csv # csv reading needed library
from nlosExclusion.msg import GNSS_Raw,GNSS_Raw_Array,exclusionSatNum # ros msg
from geometry_msgs.msg import Quaternion, Point, Pose, Twist,PoseArray # ros message needed
from PyQt4 import QtCore, QtGui
import puGNSSPosCal
import time
from novatel_msgs.msg import INSPVAX
from novatel_msgs.msg import BESTPOS
import llh2ecef # llh to ecef
import ecef2llh #ecef coordinate to llh coordinate

class puSkyplot(QMainWindow): # plot skyplot
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('puSkyplot')
        self.create_main_frame()
        self.fontSize = 20
        self.fontColor = 'g'
        # self.drawBase()
        self.on_draw()

    def save_plot(self):
        pass

    def on_about(self):
        pass

    def on_pick(self, event):
        pass

    def on_draw(self):
        # self.axes_1.clear()
        # self.axes_1.grid(True)
        self.canvas.draw()

    def drawBase(self): # draw base for skyplot
        self.axes_2.set_aspect('equal')
        self.axes_2.add_artist(Circle((0, 0), 50, color='lightgrey'))  # self.circle = Circle((0, 0), 50)
        self.axes_2.add_artist(Circle((0, 0), 40, color='silver'))
        self.axes_2.add_artist(Circle((0, 0), 30, color='darkgray'))
        self.axes_2.add_artist(Circle((0, 0), 20, color='gray'))
        self.axes_2.add_artist(Circle((0, 0), 10, color='dimgrey'))
        self.axes_2.plot([0, 50], [0, 0], linewidth='1', color='gray')  # draw a line from (0,0) to (50,0)
        self.axes_2.text(51, 0, 'E', fontdict={'size': self.fontSize, 'color': self.fontColor})  # draw the 'E'
        self.axes_2.plot([0, 0], [0, 50], linewidth='1', color='gray')  # draw a line from (0,0) to (0,50)
        self.axes_2.text(0, 51, 'N', fontdict={'size': self.fontSize, 'color': self.fontColor})  # draw the 'N'
        self.axes_2.plot([0, -50], [0, 0], linewidth='1', color='gray')  # draw a line from (0,0) to (-50,0)
        self.axes_2.text(-54, 0, 'W', fontdict={'size': self.fontSize, 'color': self.fontColor})  # draw the 'W'
        self.axes_2.plot([0, 0], [0, -50], linewidth='1', color='gray')  # draw a line from (0,0) to (0,-50)
        self.axes_2.text(0, -54, 'S', fontdict={'size': self.fontSize, 'color': self.fontColor})  # draw the 'S'
        # tilt and azimuth
        self.axes_2.plot([0, 43.3], [0, 25], linewidth='1', color='gray')
        self.axes_2.text(45, 26, '60', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, 25], [0, 43.3], linewidth='1', color='gray')
        self.axes_2.text(26, 45, '30', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, -43.3], [0, 25], linewidth='1', color='gray')
        self.axes_2.text(-51, 28, '300', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, -25], [0, 43.3], linewidth='1', color='gray')
        self.axes_2.text(-27, 48, '330', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, -43.3], [0, -25], linewidth='1', color='gray')
        self.axes_2.text(-51, -30, '240', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, -25], [0, -43.3], linewidth='1', color='gray')
        self.axes_2.text(-29, -48, '210', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, 43.3], [0, -25], linewidth='1', color='gray')
        self.axes_2.text(48, -28, '120', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.plot([0, 25], [0, -43.3], linewidth='1', color='gray')
        self.axes_2.text(27, -47, '150', fontdict={'size': self.fontSize, 'color': self.fontColor})
        # draw elevation indicators
        self.axes_2.text(2, 13, '72', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.text(2, 23, '54', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.text(2, 33, '36', fontdict={'size': self.fontSize, 'color': self.fontColor})
        self.axes_2.text(2, 43, '18', fontdict={'size': self.fontSize, 'color': self.fontColor})

    def create_main_frame(self):
        self.main_frame = QWidget()
        self.dpi = 100
        self.fig = Figure((10.0, 10.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.axes_2 = self.fig.add_subplot(111)
        self.axes_2.set_xlim([-60, 60])
        self.axes_2.set_ylim([-60, 60])
        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

class puDouDeckBusS(puSkyplot): # callback double-decker, yaw and GNSS data for creating Skyplot
    def __init__(self):
        puSkyplot.__init__(self)
        self.window_size=20
        self.calAngNS = 0.0  # define calibrationa angle of car yaw and earth north -105.57 for static test
        rospy.Subscriber('double_decker_parameters', PoseArray,self.callDouDeckBusS)  # initialize subscriber to "double_decker_parameters", which contains azimuth and elevetion of double-decker bus boundary
        rospy.Subscriber('GNSS_', GNSS_Raw_Array, self.CallGNSS_1) # subscribe the data provided by roabag file
        rospy.Subscriber('novatel_data/inspvax', INSPVAX, self.callINSPVAX)
        rospy.Subscriber('exclusionSatNum', exclusionSatNum, self.callExcludedSatNum)
        rospy.Subscriber('velodyne_points_0', PointCloud2, self.PointCloudCall)
        self.GNSS_update = 0.0
        self.azim_ = []
        self.elev_ = []
        self.satIdx = []
        self.bouSide1X = []
        self.bouSide1Y = []
        self.bouSide2X = []
        self.bouSide2Y = []
        self.bouSide1X_ = []
        self.bouSide1Y_ = []
        self.bouSide2X_ = []
        self.bouSide2Y_ = []
        self.posArr = []  # create a list to save double-decker bus boundary information
        self.excludedSatList = []
        self.satColor='b'
        self.satColorLos  = 'g'
        self.satColorNlos = 'r'
        self.maxCircRadi = 50

    def callINSPVAX(self,data):
        inspvax_ = INSPVAX()
        inspvax_ = data
        self.calAngNS = (float(360.0 - inspvax_.azimuth))

    def callExcludedSatNum(self,data):
        self.excludedSatList = []
        self.excludedSatList = data.satlist

    def callDouDeckBusS(self, data):
        self.posArr = data.poses  # save double-decker bus boundary information
        self.bouSide1X[:] = []
        self.bouSide1Y[:] = []
        self.bouSide2X[:] = []
        self.bouSide2Y[:] = []
        self.bouSide1X_[:] = []
        self.bouSide1Y_[:] = []
        self.bouSide2X_[:] = []
        self.bouSide2Y_[:] = []
        # print 'callDouble-decker S...'
        lenPosArr_ = 0.0
        lenPosArr_ = len(self.posArr)
        if (lenPosArr_ > 2):
            lenPosArr_ = 2
        for i in range(lenPosArr_):
        # for i in range(len(self.posArr)):
            # x--azimuth y--elevation
            if (self.posArr[i].orientation.x < 0):
                self.posArr[i].orientation.x = self.posArr[i].orientation.x + 360
            if (self.posArr[i].orientation.z < 0):
                self.posArr[i].orientation.z = self.posArr[i].orientation.z + 360
                    # -----------for static version--------------#
                    # Cable of Velodyne 32 point the body who hold the LiDAR, in the previous static experiment, the
                    # bus is on the other side of the cable. Thus, -1 is needed to multiply into the transform
            # self.bouSide1X.append(-1 * (self.posArr[i].orientation.y * (-0.55556) + 50.0) * np.cos((self.posArr[i].orientation.x - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # start point azimuth
            # self.bouSide1Y.append(     (self.posArr[i].orientation.y * (-0.55556) + 50.0) * np.sin(-1 * (self.posArr[i].orientation.x - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # start point elevation
            # self.bouSide2X.append(-1 * (self.posArr[i].orientation.w * (-0.55556) + 50.0) * np.cos(-1 * (self.posArr[i].orientation.z - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # end point azimuth
            # self.bouSide2Y.append(     (self.posArr[i].orientation.w * (-0.55556) + 50.0) * np.sin(-1 * (self.posArr[i].orientation.z - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # end point elevation
                    # -----------for dynamic version--------------#
                    # Cable of Velodyne 32 point the body who hold the LiDAR, in the previous static experiment, the
            # bus is on the other side of the cable
            self.bouSide1X.append(1 * (self.posArr[i].orientation.y * (-0.55556) + 50.0) * np.cos((self.posArr[i].orientation.x - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # start point azimuth
            self.bouSide1Y.append(1 * (self.posArr[i].orientation.y * (-0.55556) + 50.0) * np.sin((self.posArr[i].orientation.x - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # start point elevation

            # print 'start point azimuth in skyplot===', (math.atan2(self.bouSide1Y[-1] , self.bouSide1X[-1])) * 180.0 / (math.pi),'self.calAngNS',self.calAngNS

            self.bouSide2X.append(1 * (self.posArr[i].orientation.w * (-0.55556) + 50.0) * np.cos((self.posArr[i].orientation.z - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # end point azimuth
            self.bouSide2Y.append(1 * (self.posArr[i].orientation.w * (-0.55556) + 50.0) * np.sin((self.posArr[i].orientation.z - 90.0 + self.calAngNS) * 3.14159 / 180.0))  # end point elevation

            uniVecSid1X_= self.bouSide1X[i] / (math.sqrt(self.bouSide1X[i] * self.bouSide1X[i] + self.bouSide1Y[i] * self.bouSide1Y[i]))
            uniVecSid1Y_ = self.bouSide1Y[i] / (math.sqrt(self.bouSide1X[i] * self.bouSide1X[i] + self.bouSide1Y[i] * self.bouSide1Y[i]))
            uniVecSid2X_ = self.bouSide2X[i] / (math.sqrt(self.bouSide2X[i] * self.bouSide2X[i] + self.bouSide2Y[i] * self.bouSide2Y[i]))
            uniVecSid2Y_ = self.bouSide2Y[i] / (math.sqrt(self.bouSide2X[i] * self.bouSide2X[i] + self.bouSide2Y[i] * self.bouSide2Y[i]))
            self.bouSide1X_.append(uniVecSid1X_ * self.maxCircRadi)
            self.bouSide1Y_.append(uniVecSid1Y_ * self.maxCircRadi)
            self.bouSide2X_.append(uniVecSid2X_ * self.maxCircRadi)
            self.bouSide2Y_.append(uniVecSid2Y_ * self.maxCircRadi)

    def CallGNSS_1(self, data):  # GNSS data
        self.GNSS_1 = GNSS_Raw_Array()
        self.GNSS_1 = data
        self.azim_[:] = []
        self.elev_[:] = []
        self.satIdx[:] = []
        self.satecirclRad = 3
        for index_1 in range(len(self.GNSS_1.GNSS_Raws)):  # index all the satellite information in one epoch
            self.azim_.append( (self.GNSS_1.GNSS_Raws[index_1].elevation*(-0.55556)+50.0)* np.cos(-1*(self.GNSS_1.GNSS_Raws[index_1].azimuth-90.0)*3.14159/180.0))
            self.elev_.append( (self.GNSS_1.GNSS_Raws[index_1].elevation*(-0.55556)+50.0)* np.sin(-1*(self.GNSS_1.GNSS_Raws[index_1].azimuth-90.0)*3.14159/180.0))
            self.satIdx.append(self.GNSS_1.GNSS_Raws[index_1].prn_satellites_index)
        self.GNSS_update = 1
        # self.drawAll()

    def PointCloudCall(self,data):
        pointCloud_ = PointCloud2()
        pointCloud_ = data
        if(self.GNSS_update ==1): # make sure that GNSS is updated
            self.drawAll()
            self.GNSS_update =0.0


    def drawAll(self):
        self.drawBase()
        for sate_index in range(len(self.azim_)):  # draw the satellite into the Skyplot (template: circle+G30)
            if(len(self.azim_) > 0):
                self.axes_2.add_artist(Circle((self.azim_[sate_index], self.elev_[sate_index]), self.satecirclRad, color=self.satColorLos))  # self.circle = Circle((0, 0), 50)
                self.axes_2.text(self.azim_[sate_index], self.elev_[sate_index], str(int(self.satIdx[sate_index])), fontdict={'size': self.fontSize, 'color': self.satColor})  # draw the 'E'

            if ((len(self.satIdx) > 0) and ((self.satIdx[sate_index]) in self.excludedSatList)):
                self.axes_2.add_artist(Circle((self.azim_[sate_index], self.elev_[sate_index]), self.satecirclRad,
                                              color=self.satColorNlos))  # self.circle = Circle((0, 0), 50)
                self.axes_2.text(self.azim_[sate_index], self.elev_[sate_index], str(int(self.satIdx[sate_index])),
                                 fontdict={'size': self.fontSize, 'color': self.satColor})  # draw the 'E'
            else:
                self.axes_2.add_artist(Circle((self.azim_[sate_index], self.elev_[sate_index]), self.satecirclRad,
                                              color=self.satColorLos))  # self.circle = Circle((0, 0), 50)
                self.axes_2.text(self.azim_[sate_index], self.elev_[sate_index], str(int(self.satIdx[sate_index])),
                                 fontdict={'size': self.fontSize, 'color': self.satColor})  # draw the 'E'

        self.axes_2.plot([self.bouSide1X, self.bouSide2X], [self.bouSide1Y, self.bouSide2Y], linewidth='2', color='fuchsia')  # draw a line from (0,0) to (50,0)
        self.axes_2.plot([self.bouSide1X, self.bouSide1X_], [self.bouSide1Y, self.bouSide1Y_], linewidth='1',color='fuchsia')  # draw a line from (0,0) to (50,0)
        self.axes_2.plot([self.bouSide2X, self.bouSide2X_], [self.bouSide2Y, self.bouSide2Y_], linewidth='1',color='fuchsia')  # draw a line from (0,0) to (50,0)
        self.canvas.draw()
        self.cleanS()

    def cleanS(self):
        self.axes_2.clear()




class nlosExclusionS_(QMainWindow): # by timer   For paper NLOS exclusion caused by double-decekr bus

    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('evaluation')
        rospy.Subscriber('/novatel_data/bestpos', BESTPOS, self.callcptBestPos_)

    def save_plot(self):
        pass

    def on_about(self):
        pass

    def on_pick(self, event):
        pass

    def on_draw(self):
        # self.axes_1.clear()
        # self.axes_1.grid(True)
        self.canvas.draw()

    def create_main_frame(self):
        self.main_frame = QWidget()
        self.dpi = 100
        self.fig = Figure((5.0, 20.0), dpi=self.dpi) # 5.0 4.0
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        # self.axes = self.fig.add_subplot(111, projection='3d')
        self.axesCurv_1 = self.fig.add_subplot(311)
        self.axesCurv_2 = self.fig.add_subplot(312)
        self.axesCurv_3 = self.fig.add_subplot(313)
        # self.axes.get_zaxis().set_visible(False)


        self.canvas.mpl_connect('pick_event', self.on_pick)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        vbox = QVBoxLayout()
        vbox.addWidget(self.canvas)
        vbox.addWidget(self.mpl_toolbar)
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)

    def PercentCal(self,data_):
        # self.axes_1.grid(True)
        error_ = []
        error_ = data_
        percent_1 = 0.0 # >15
        percent_2 = 0.0 # >25
        percent_3 = 0.0 # >40
        for perCal in range(len(error_)):
            if(error_[perCal]<=15):
                percent_1 = percent_1 + 1
            if (error_[perCal] <= 30):
                percent_2 = percent_2 + 1
            if (error_[perCal] >= 40):
                percent_3 = percent_3 + 1
        percent_1 = percent_1 / len(error_)
        percent_2 = percent_2 / len(error_)
        percent_3 = percent_3 / len(error_)
        print 'percent_1=',percent_1,'percent_2=',percent_2,'percent_3=',percent_3

if __name__ == '__main__':
    app = QApplication(sys.argv)
    rospy.init_node('puSkyplot', anonymous=True)
    puDouDeckBusS_=puDouDeckBusS()
    print 'this is skyplot from Positioning and navigation laboratory'
    puDouDeckBusS_.show()
    app.exec_()