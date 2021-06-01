#!/usr/bin/env python
# license removed for brevity
"""
    GNSS positioning calculation
    Welson Wen, Ph.D.
    https://sites.google.com/view/weisongwen/about-me
"""
from numpy import * # numpy needed
import matplotlib as mpl #plot needed
mpl.use("TkAgg") # Use TKAgg to show figures:set this to show plot
import matplotlib.pyplot as plt #plotting
import pandas as pd # pandas needed renamed as pd
import numpy as np #numpy needed renamed as np
import geometry_msgs.msg as gm #ros geometry message
from geometry_msgs.msg import Quaternion, Point, Pose, Twist,PoseArray # commonly used message type
from sensor_msgs.msg   import NavSatFix # standard message type for GNSSs
from nlosExclusion.msg import GNSS_Raw_Array,GNSS_Raw # customerized ros message type
from matplotlib.patches import Ellipse, Circle # draw circle needs library
import csv # csv reading needed library
import datetime #time format (datetime)
import time #time format (time)
import llh2ecef # llh to ecef
import ecef2llh #ecef coordinate to llh coordinate
from nlosExclusion.msg import Satellite_Info # customized ros message type Satellite_Info containing satellites exclusion numbers
import rospy
from novatel_msgs.msg import BESTPOS

class GNSSPosCal():
    def __init__(self):
        self.calMode = 'LSGPS' #  'LSGPS' 'LSGNSS'
        self.GNSSTim = 0
        self.dop = 0  # dop
        self.toSv = 0
        self.iterations_=0
        self.appInfo = 0
        self.prn = []  # prn
        self.snr = []  # snr
        self.pseudResid = {}  # pseudorange residual
        self.visi = [] # visibility
        self.azimuth = []
        self.elevation =[]
        self.ecef_=[]  # calculation result
        self.llh_ =[]
        # self.GroTruth = [22.303953,114.181925,14.0] # for experiment1: small NLOS reception
        self.GroTruth = [22.303756,114.18215,14.0] # for experiment2: Big NLOS reception
        self.error_ = 0

        self.GPSNum = []  # create a list to save GPS satellites numbers
        self.BeiNum = []  # create a list to save Beidou satellites numbers


    def LSGPSPosCalStep(self,GNSS_one_epoch_, init_x, init_y, init_z, init_b):  # least square
        GNSS_one_epoch = GNSS_Raw_Array()
        GNSS_one_epoch = GNSS_one_epoch_
        rec_ini_pos_x = float(init_x)  # initial position x (m)
        rec_ini_pos_y = float(init_y)  # initial position y (m)
        rec_ini_pos_z = float(init_z)  # initial position z (m)
        rec_clo_bia_b = float(init_b)  # initial distance bias caused by clock bias (m)
        devia_xyz = 1.0  # initial set receiver position estimation error (m)
        gues_pseud = []  # guessed pseudorange
        pseud_error = []  # pseudorange error
        G_matrix_ = array([2, 2, 2, 1], dtype='float')  # creat G matrix to save transform parameters
        self.dop = self.DopCalculation(GNSS_one_epoch)
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            if((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                    GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37))):
                if(self.iterations_<1): # save one time only
                    self.GNSSTim = GNSS_one_epoch.GNSS_Raws[index_1].GNSS_time
                    self.toSv = GNSS_one_epoch.GNSS_Raws[index_1].total_sv
                    self.prn.append(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index)
                    self.snr.append(GNSS_one_epoch.GNSS_Raws[index_1].snr)
                    self.visi.append(GNSS_one_epoch.GNSS_Raws[index_1].visable)
                    self.azimuth.append(GNSS_one_epoch.GNSS_Raws[index_1].azimuth)
                    self.elevation.append(GNSS_one_epoch.GNSS_Raws[index_1].elevation)
                sx_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x)  # satellite position x
                sy_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y)  # satellite position y
                sz_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z)  # satellite position z
                sx_1 = (sx_1 - rec_ini_pos_x) * (sx_1 - rec_ini_pos_x)  # satellite to receiver distance in x idrection
                sy_1 = (sy_1 - rec_ini_pos_y) * (sy_1 - rec_ini_pos_y)  # satellite to receiver distance in y idrection
                sz_1 = (sz_1 - rec_ini_pos_z) * (sz_1 - rec_ini_pos_z)  # satellite to receiver distance in z idrection
                sat2rec_dis = sqrt(sx_1 + sy_1 + sz_1)  # guessed pseudorange
                gues_pseud.append(sat2rec_dis)  # save guessed pseudorange
                pseud_error_element = float(GNSS_one_epoch.GNSS_Raws[index_1].pseudorange) - float(
                    sat2rec_dis) + float(rec_clo_bia_b)  # pseudorange error
                pseud_error.append(pseud_error_element)  # save pseudorange error
                G_row = []  # G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x - rec_ini_pos_x) / float(
                    sat2rec_dis) * -1)  # x for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y - rec_ini_pos_y) / float(
                    sat2rec_dis) * -1)  # y for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z - rec_ini_pos_z) / float(
                    sat2rec_dis) * -1)  # z for G matrix row
                element_float = 1.0  # last element for each row
                G_row.append(element_float)  # save last element for each row
                G_matrix_ = np.row_stack((G_matrix_, G_row))  # add each row to G_matrix
                del G_row[:]  # relief G_row
        # get pseudorange error
        pseud_error_mat = np.array(pseud_error)  # from list to array
        pseud_error_mat = pseud_error_mat.transpose()  # transpose
        # get G matrix
        G_matrix_ = np.delete(G_matrix_, [0], axis=0)  # delete the first row of G matrix

        delta_p = np.dot((G_matrix_.transpose()), G_matrix_)  # G(T) * G
        delta_p_2 = np.linalg.inv(delta_p)  # inverse matrix of G(T) * G
        delta_p = np.dot(delta_p_2, (G_matrix_.transpose()))  # multiply (inverse matrix of G(T) * G) and G(T)
        delta_p = np.dot(delta_p, pseud_error_mat)  # multiply with pseud_error_mat
        rec_ini_pos_x = rec_ini_pos_x + float(delta_p[0])  # update receiver position in x direction
        rec_ini_pos_y = rec_ini_pos_y + float(delta_p[1])  # update receiver position in y idrection
        rec_ini_pos_z = rec_ini_pos_z + float(delta_p[2])  # update receiver position in z idrection
        rec_clo_bia_b = rec_clo_bia_b + float(delta_p[3])  # update receiver clock bias in meters
        devia_x = float(delta_p[0])  # save delta x
        devia_y = float(delta_p[1])  # save delta y
        devia_z = float(delta_p[2])  # save delta z
        devia_b = float(delta_p[3])  # save delta bias
        devia_xyz = sqrt(devia_x * devia_x + devia_y * devia_y + devia_z * devia_z)  # get total bias
        # print 'delta_p',delta_p
        # print 'position estimation x=',rec_ini_pos_x
        # print 'position estimation y=', rec_ini_pos_y
        # print 'position estimation Z=', rec_ini_pos_z
        # print 'position estimation b=', rec_clo_bia_b
        # print 'position estimation devia_xyz=', devia_xyz
        del gues_pseud[:]  # relief gues_pseud[] list
        del pseud_error[:]  # relief pseud_error[] list
        return float(rec_ini_pos_x), float(rec_ini_pos_y), float(rec_ini_pos_z), float(rec_clo_bia_b), float(devia_xyz)

    '''
            GPS:     1:32
            GLONASS: 32 + 1:24
            Galileo: 57 + 1:30
            Beidou:  87 + 1:37
            QZSS:    124 + 1:4
    '''

    def LSGNSSPosCalStep(self,GNSS_one_epoch_, init_x, init_y, init_z, init_b_GPS,
                                                           init_b_Beidou):  # least square for hybrid GNSS positioning (GPS + Beidou)
        GNSS_one_epoch = GNSS_Raw_Array()
        GNSS_one_epoch = GNSS_one_epoch_
        rec_ini_pos_x = float(init_x)  # initial position x (m)
        rec_ini_pos_y = float(init_y)  # initial position y (m)
        rec_ini_pos_z = float(init_z)  # initial position z (m)
        rec_clo_bia_b_GPS = float(init_b_GPS)  # initial distance bias caused by clock bias of GPS (m)
        rec_clo_bia_b_Beidou = float(init_b_Beidou)  # initial distance bias caused by clock bias of Beidou (m)
        devia_xyz = 1.0  # initial set receiver position estimation error (m)
        gues_pseud = []  # guessed pseudorange
        pseud_error = []  # pseudorange error
        G_matrix_ = array([2, 2, 2, 1, 1], dtype='float')  # creat G matrix to save transform parameters
        self.dop = self.DopCalculation(GNSS_one_epoch)
        # get guessed pseudorange and pseudorange error
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                    GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37)) or
                    ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                            GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32))):
                if (self.iterations_ < 1): # save one time only
                    self.GNSSTim = GNSS_one_epoch.GNSS_Raws[index_1].GNSS_time
                    self.toSv = GNSS_one_epoch.GNSS_Raws[index_1].total_sv
                    self.prn.append(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index)
                    self.snr.append(GNSS_one_epoch.GNSS_Raws[index_1].snr)
                    self.visi.append(GNSS_one_epoch.GNSS_Raws[index_1].visable)
                    self.azimuth.append(GNSS_one_epoch.GNSS_Raws[index_1].azimuth)
                    self.elevation.append(GNSS_one_epoch.GNSS_Raws[index_1].elevation)
                sx_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x)  # satellite position x
                sy_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y)  # satellite position y
                sz_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z)  # satellite position z
                # print 'satellite index',GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index
                sx_1 = (sx_1 - rec_ini_pos_x) * (
                        sx_1 - rec_ini_pos_x)  # satellite to receiver distance in x idrection
                sy_1 = (sy_1 - rec_ini_pos_y) * (
                        sy_1 - rec_ini_pos_y)  # satellite to receiver distance in y idrection
                sz_1 = (sz_1 - rec_ini_pos_z) * (
                        sz_1 - rec_ini_pos_z)  # satellite to receiver distance in z idrection
                sat2rec_dis = 0.0  # initialize variable
                sat2rec_dis = sqrt(sx_1 + sy_1 + sz_1)  # guessed pseudorange
                gues_pseud.append(sat2rec_dis)  # save guessed pseudorange
                G_row = []  # G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x - rec_ini_pos_x) / float(
                    sat2rec_dis) * -1)  # x for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y - rec_ini_pos_y) / float(
                    sat2rec_dis) * -1)  # y for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z - rec_ini_pos_z) / float(
                    sat2rec_dis) * -1)  # z for G matrix row
                if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                        GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37))):
                    element_float_GPS = 0.0  # GPS element for each row
                    element_float_Beidou = 1.0  # Beidou element for each row
                    G_row.append(element_float_GPS)  # save last two element for each row
                    G_row.append(element_float_Beidou)  # save last two element for each row
                    pseud_error_element = float(GNSS_one_epoch.GNSS_Raws[index_1].pseudorange) - float(
                        sat2rec_dis) + float(rec_clo_bia_b_Beidou)  # Beidou pseudorange error
                    pseud_error.append(pseud_error_element)  # save pseudorange error
                if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                        GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32)):
                    element_float_GPS = 1.0  # GPS element for each row
                    element_float_Beidou = 0.0  # Beidou element for each row
                    G_row.append(element_float_GPS)  # save last two element for each row
                    G_row.append(element_float_Beidou)  # save last two element for each row
                    pseud_error_element = float(GNSS_one_epoch.GNSS_Raws[index_1].pseudorange) - float(
                        sat2rec_dis) + float(rec_clo_bia_b_GPS)  # GPS pseudorange error
                    pseud_error.append(pseud_error_element)  # save pseudorange error
                # print 'length of G_row',len(G_row),'GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index',GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index
                G_matrix_ = np.row_stack((G_matrix_, G_row))  # add each row to G_matrix
                del G_row[:]  # relief G_row
        # get pseudorange error
        pseud_error_mat = np.array(pseud_error)  # from list to array
        pseud_error_mat = pseud_error_mat.transpose()  # transpose
        # get G matrix
        G_matrix_ = np.delete(G_matrix_, [0], axis=0)  # delete the first row of G matrix
        # get cofactor matrix
        # cofactorMat_ = np.array(self.cofactorMatrixCal(GNSS_one_epoch_))
        # cofactorMat_ = np.diag(cofactorMat_)  # diag matrix
        # print 'cofactors', self.cofactorMatrixCal(GNSS_one_epoch_), cofactorMat_

        delta_p = np.dot((G_matrix_.transpose()), G_matrix_)  # G(T) * G
        delta_p_2 = np.linalg.inv(delta_p)  # inverse matrix of G(T) * G
        delta_p = np.dot(delta_p_2, (G_matrix_.transpose()))  # multiply (inverse matrix of G(T) * G) and G(T)
        delta_p = np.dot(delta_p, pseud_error_mat)  # multiply with pseud_error_mat
        rec_ini_pos_x = rec_ini_pos_x + float(delta_p[0])  # update receiver position in x direction
        rec_ini_pos_y = rec_ini_pos_y + float(delta_p[1])  # update receiver position in y idrection
        rec_ini_pos_z = rec_ini_pos_z + float(delta_p[2])  # update receiver position in z idrection
        rec_clo_bia_b_GPS = rec_clo_bia_b_GPS + float(delta_p[3])  # update receiver clock bias of GPS in meters
        rec_clo_bia_b_Beidou = rec_clo_bia_b_Beidou + float(
            delta_p[4])  # update receiver clock bias of Beidou in meters
        devia_x = float(delta_p[0])  # save delta x
        devia_y = float(delta_p[1])  # save delta y
        devia_z = float(delta_p[2])  # save delta z
        devia_b_GPS = float(delta_p[3])  # save delta bias of GPS
        devia_b_Beidou = float(delta_p[4])  # save delta bias of Beidou
        devia_xyz = sqrt(devia_x * devia_x + devia_y * devia_y + devia_z * devia_z)  # get total bias
        # print 'delta_p',delta_p
        # print 'position estimation x=',rec_ini_pos_x
        # print 'position estimation y=', rec_ini_pos_y
        # print 'position estimation Z=', rec_ini_pos_z
        # print 'position estimation b=', rec_clo_bia_b
        # print 'position estimation devia_xyz=', devia_xyz
        del gues_pseud[:]  # relief gues_pseud[] list
        del pseud_error[:]  # relief pseud_error[] list
        return float(rec_ini_pos_x), float(rec_ini_pos_y), float(rec_ini_pos_z), float(rec_clo_bia_b_GPS), float(
            rec_clo_bia_b_Beidou), float(devia_xyz)

    def WLSGNSSPosCalStep(self,GNSS_one_epoch_, init_x, init_y, init_z, init_b_GPS,
                                                           init_b_Beidou):  # least square for hybrid GNSS positioning (GPS + Beidou)
        GNSS_one_epoch = GNSS_Raw_Array()
        GNSS_one_epoch = GNSS_one_epoch_
        rec_ini_pos_x = float(init_x)  # initial position x (m)
        rec_ini_pos_y = float(init_y)  # initial position y (m)
        rec_ini_pos_z = float(init_z)  # initial position z (m)
        rec_clo_bia_b_GPS = float(init_b_GPS)  # initial distance bias caused by clock bias of GPS (m)
        rec_clo_bia_b_Beidou = float(init_b_Beidou)  # initial distance bias caused by clock bias of Beidou (m)
        devia_xyz = 1.0  # initial set receiver position estimation error (m)
        gues_pseud = []  # guessed pseudorange
        pseud_error = []  # pseudorange error
        G_matrix_ = array([2, 2, 2, 1, 1], dtype='float')  # creat G matrix to save transform parameters
        self.dop = self.DopCalculation(GNSS_one_epoch)
        # get guessed pseudorange and pseudorange error
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                    GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37)) or
                    ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                            GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32))):
                if (self.iterations_ < 1): # save one time only
                    self.GNSSTim = GNSS_one_epoch.GNSS_Raws[index_1].GNSS_time
                    self.toSv = GNSS_one_epoch.GNSS_Raws[index_1].total_sv
                    self.prn.append(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index)
                    self.snr.append(GNSS_one_epoch.GNSS_Raws[index_1].snr)
                    self.visi.append(GNSS_one_epoch.GNSS_Raws[index_1].visable)
                    self.azimuth.append(GNSS_one_epoch.GNSS_Raws[index_1].azimuth)
                    self.elevation.append(GNSS_one_epoch.GNSS_Raws[index_1].elevation)
                sx_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x)  # satellite position x
                sy_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y)  # satellite position y
                sz_1 = float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z)  # satellite position z
                # print 'satellite index',GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index
                sx_1 = (sx_1 - rec_ini_pos_x) * (
                        sx_1 - rec_ini_pos_x)  # satellite to receiver distance in x idrection
                sy_1 = (sy_1 - rec_ini_pos_y) * (
                        sy_1 - rec_ini_pos_y)  # satellite to receiver distance in y idrection
                sz_1 = (sz_1 - rec_ini_pos_z) * (
                        sz_1 - rec_ini_pos_z)  # satellite to receiver distance in z idrection
                sat2rec_dis = 0.0  # initialize variable
                sat2rec_dis = sqrt(sx_1 + sy_1 + sz_1)  # guessed pseudorange
                gues_pseud.append(sat2rec_dis)  # save guessed pseudorange
                G_row = []  # G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x - rec_ini_pos_x) / float(
                    sat2rec_dis) * -1)  # x for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y - rec_ini_pos_y) / float(
                    sat2rec_dis) * -1)  # y for G matrix row
                G_row.append(float(GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z - rec_ini_pos_z) / float(
                    sat2rec_dis) * -1)  # z for G matrix row
                if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                        GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37))):
                    element_float_GPS = 0.0  # GPS element for each row
                    element_float_Beidou = 1.0  # Beidou element for each row
                    G_row.append(element_float_GPS)  # save last two element for each row
                    G_row.append(element_float_Beidou)  # save last two element for each row
                    pseud_error_element = float(GNSS_one_epoch.GNSS_Raws[index_1].pseudorange) - float(
                        sat2rec_dis) + float(rec_clo_bia_b_Beidou)  # Beidou pseudorange error
                    pseud_error.append(pseud_error_element)  # save pseudorange error
                if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                        GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32)):
                    element_float_GPS = 1.0  # GPS element for each row
                    element_float_Beidou = 0.0  # Beidou element for each row
                    G_row.append(element_float_GPS)  # save last two element for each row
                    G_row.append(element_float_Beidou)  # save last two element for each row
                    pseud_error_element = float(GNSS_one_epoch.GNSS_Raws[index_1].pseudorange) - float(
                        sat2rec_dis) + float(rec_clo_bia_b_GPS)  # GPS pseudorange error
                    pseud_error.append(pseud_error_element)  # save pseudorange error
                # print 'length of G_row',len(G_row),'GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index',GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index
                G_matrix_ = np.row_stack((G_matrix_, G_row))  # add each row to G_matrix
                del G_row[:]  # relief G_row
        # get pseudorange error
        pseud_error_mat = np.array(pseud_error)  # from list to array
        pseud_error_mat = pseud_error_mat.transpose()  # transpose
        # get G matrix
        G_matrix_ = np.delete(G_matrix_, [0], axis=0)  # delete the first row of G matrix
        # get cofactor matrix
        cofactorMat_ = np.array(self.cofactorMatrixCal(GNSS_one_epoch))
        cofactorMat_ = np.diag(cofactorMat_) # diag matrix
        # print 'cofactors',self.cofactorMatrixCal(GNSS_one_epoch_),cofactorMat_,len(GNSS_one_epoch.GNSS_Raws)

        delta_p = np.dot((G_matrix_.transpose()), cofactorMat_)  # G(T) * G
        delta_p = np.dot(delta_p, G_matrix_)  # G(T) * G
        delta_p_2 = np.linalg.inv(delta_p)  # inverse matrix of G(T) * G
        delta_p = np.dot(delta_p_2, (G_matrix_.transpose()))  # multiply (inverse matrix of G(T) * G) and G(T)
        delta_p = np.dot(delta_p, cofactorMat_)  # multiply (inverse matrix of G(T) * G) and G(T)
        delta_p = np.dot(delta_p, pseud_error_mat)  # multiply with pseud_error_mat
        rec_ini_pos_x = rec_ini_pos_x + float(delta_p[0])  # update receiver position in x direction
        rec_ini_pos_y = rec_ini_pos_y + float(delta_p[1])  # update receiver position in y idrection
        rec_ini_pos_z = rec_ini_pos_z + float(delta_p[2])  # update receiver position in z idrection
        rec_clo_bia_b_GPS = rec_clo_bia_b_GPS + float(delta_p[3])  # update receiver clock bias of GPS in meters
        rec_clo_bia_b_Beidou = rec_clo_bia_b_Beidou + float(
            delta_p[4])  # update receiver clock bias of Beidou in meters
        devia_x = float(delta_p[0])  # save delta x
        devia_y = float(delta_p[1])  # save delta y
        devia_z = float(delta_p[2])  # save delta z
        devia_b_GPS = float(delta_p[3])  # save delta bias of GPS
        devia_b_Beidou = float(delta_p[4])  # save delta bias of Beidou
        devia_xyz = sqrt(devia_x * devia_x + devia_y * devia_y + devia_z * devia_z)  # get total bias
        # print 'delta_p',delta_p
        # print 'position estimation x=',rec_ini_pos_x
        # print 'position estimation y=', rec_ini_pos_y
        # print 'position estimation Z=', rec_ini_pos_z
        # print 'position estimation b=', rec_clo_bia_b
        # print 'position estimation devia_xyz=', devia_xyz
        del gues_pseud[:]  # relief gues_pseud[] list
        del pseud_error[:]  # relief pseud_error[] list
        return float(rec_ini_pos_x), float(rec_ini_pos_y), float(rec_ini_pos_z), float(rec_clo_bia_b_GPS), float(
            rec_clo_bia_b_Beidou), float(devia_xyz)

    def iterPosCal(self,GNSS_one_epoch_,calMode):
        iterations = 0
        itera_x = 0
        itera_y = 0
        itera_z = 0
        self.calMode = calMode
        self.getSatNum(GNSS_one_epoch_) # get satellites numbers
        if(self.calMode=='LSGNSS') and (self.GPSNum>=4) and (self.BeiNum>1):
            itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou, itera_bias_ = self.LSGNSSPosCalStep(
                GNSS_one_epoch_, 0.1, 0.1, 0.1, 1.0, 1.0)  # first iteration from (0.1, 0.1, 0.1, 1.0)
            self.iterations_ = self.iterations_+1
            while (itera_bias_ > 1e-4) and iterations < 10:  # threshold for iterations:value and times
                itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou, itera_bias_ = self.LSGNSSPosCalStep(
                    GNSS_one_epoch_, itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou)  # iteration
                self.iterations_ = self.iterations_ + 1
                iterations = iterations + 1  # add one iteration
        elif (self.calMode == 'LSGPS') and (self.GPSNum>=4):
            itera_x, itera_y, itera_z, itera_b, itera_bias_ = self.LSGPSPosCalStep(
                GNSS_one_epoch_, 0.1,
                0.1, 0.1,
                1.0)  # first iteration from (0.1, 0.1, 0.1, 1.0)
            self.iterations_ = self.iterations_ + 1
            while (itera_bias_ > 1e-4) and iterations < 10:  # threshold for iterations:value and times
                itera_x, itera_y, itera_z, itera_b, itera_bias_ = self.LSGPSPosCalStep(
                    GNSS_one_epoch_,itera_x, itera_y,itera_z,itera_b)
                self.iterations_ = self.iterations_ + 1
                iterations = iterations + 1  # add one iteration
        elif (self.calMode == 'WLSGNSS') and (self.GPSNum >= 4) and (self.BeiNum > 1):
            itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou, itera_bias_ = self.WLSGNSSPosCalStep(
                GNSS_one_epoch_, 0.1, 0.1, 0.1, 1.0, 1.0)  # first iteration from (0.1, 0.1, 0.1, 1.0)
            self.iterations_ = self.iterations_ + 1
            while (itera_bias_ > 1e-4) and iterations < 10:  # threshold for iterations:value and times
                itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou, itera_bias_ = self.WLSGNSSPosCalStep(
                    GNSS_one_epoch_, itera_x, itera_y, itera_z, itera_b_GPS, itera_b_Beidou)  # iteration
                self.iterations_ = self.iterations_ + 1
                iterations = iterations + 1  # add one iteration
                print 'itera_b_GPS',itera_b_GPS,'itera_b_Beidou-----------------------------------------',itera_b_Beidou
        
        self.ecef_.append(float(itera_x))
        self.ecef_.append(float(itera_y))
        self.ecef_.append(float(itera_z))
        self.pseudoResCal(GNSS_one_epoch_)
        self.PosError()
        self.llh_ = ecef2llh.xyz2llh(self.ecef_)  # ecef to llh coordinates
        iterations = 0.0  # initialize iterations variable

        # GNSSPosCal_    = puGNSSPosCal.GNSSPosCal()
        # GNSSPosCal_.iterGNSSPosCal(self.GNSSArr)

    def DopCalculation(self,GNSS_one_epoch_):  # get Dop in one epoch
        GNSS_one_epoch = GNSS_Raw_Array()
        GNSS_one_epoch = GNSS_one_epoch_
        H_matrix_ = array([2, 2, 2, 1], dtype='float')  # creat H matrix to save transform parameters
        Hdop_ = 0.0  # create an variable to save Hdop
        elemin = 15.0  # minimun elevation angle
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            if (GNSS_one_epoch.GNSS_Raws[index_1].elevation <= elemin):
                # print 'satellite elevation less than 15 degree=',GNSS_one_epoch.GNSS_Raws[index_1].elevation
                continue
            cosel = float(cos(GNSS_one_epoch.GNSS_Raws[index_1].elevation))
            sinel = float(sin(GNSS_one_epoch.GNSS_Raws[index_1].elevation))
            H_row = []  # H matrix row
            H_row.append(float(cosel * sin(GNSS_one_epoch.GNSS_Raws[index_1].azimuth)))
            H_row.append(float(cosel * cos(GNSS_one_epoch.GNSS_Raws[index_1].azimuth)))
            H_row.append(float(sinel))
            H_row.append(1.0)
            H_matrix_ = np.row_stack((H_matrix_, H_row))  # add each row to H_matrix
            del H_row[:]  # relief H_row
        # get H matrix
        H_matrix_ = np.delete(H_matrix_, [0], axis=0)  # delete the first row of H matrix
        # print 'H_matrix_',H_matrix_
        Q_matrix_ = np.dot((H_matrix_.transpose()), H_matrix_)  # H(T) * G
        Q_matrix_ = np.linalg.inv(Q_matrix_)  # inverse matrix of H(T) * G
        Hdop = float(sqrt(Q_matrix_[0, 0] + Q_matrix_[1, 1]))
        # print 'Q_matrix_', Q_matrix_, 'Hdop', Hdop
        return float(Hdop)  # return result

    def getSatNum(self,GNSS_one_epoch): # get number of GPS and Beidou satellites in one epoch and save all the satellite number in list
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32)): # GPS satellites index range
                self.GPSNum.append(float(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index))
            if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37))): # Beidou satellites index range
                self.BeiNum.append(float(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index))

    def PosError(self):
        xyzTru_ = llh2ecef.llh2xyz(self.GroTruth)
        self.error_ = math.sqrt((self.ecef_[0]-xyzTru_[0]) * (self.ecef_[0]-xyzTru_[0]) + (self.ecef_[1]-xyzTru_[1]) * (self.ecef_[1]-xyzTru_[1]) + (self.ecef_[2]-xyzTru_[2]) * (self.ecef_[2]-xyzTru_[2]))
        if(self.error_ > 100):
            self.error_ = 100

    def pseudoResCal(self,GNSS_one_epoch):
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):  # index all the satellite information in one epoch
            pseudoRes_ = 0.0
            res_x = GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_x - self.ecef_[0]
            res_y = GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_y - self.ecef_[1]
            res_z = GNSS_one_epoch.GNSS_Raws[index_1].sat_pos_z - self.ecef_[2]
            pseudoRes_ = int ((math.sqrt(res_x * res_x + res_y * res_y + res_z * res_z)) - GNSS_one_epoch.GNSS_Raws[index_1].pseudorange)
            satIdx_ = float(GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index)
            self.pseudResid[str(satIdx_)] = pseudoRes_
        # print 'self.pseudResid=',self.pseudResid

    def cofactorMatrixCal(self,GNSS_one_epoch):
        snr_1 = 50.0 # T = 50
        snr_A = 30.0 # A = 30
        snr_a = 30.0 # a = 30
        snr_0 = 10.0 # F = 10
        cofactor_ = [] # cofactor of satellite
        for index_1 in range(len(GNSS_one_epoch.GNSS_Raws)):
            if ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 88) and (
                    GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= (87 + 37)) or
                    ((GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index >= 1) and (
                            GNSS_one_epoch.GNSS_Raws[index_1].prn_satellites_index <= 32))):
                snr_R = GNSS_one_epoch.GNSS_Raws[index_1].snr
                elR = GNSS_one_epoch.GNSS_Raws[index_1].elevation
                q_R_1 = 1 / (( sin(elR * pi/180.0 )) ** 2)
                # q_R_1 = 1/q_R_1
                q_R_2 = 10 ** (-(snr_R - snr_1) / snr_a)
                q_R_3 = (((snr_A / (10 ** (-(snr_0 - snr_1) / snr_a)) - 1) / (snr_0 - snr_1)) * (snr_R - snr_1) + 1)
                # q_R = float(1 / (( sin(elR * pi/180.0 )) ** 2) * (10 ** (-(snr_R - snr_1) / snr_a) * (
                #         (snr_A / (10 ** (-(snr_0 - snr_1) / snr_a)) - 1) / (snr_0 - snr_1) * (snr_R - snr_1) + 1)))
                q_R = q_R_1* (q_R_2 * q_R_3)
                cofactor_.append(float(1.0/q_R))
        return cofactor_