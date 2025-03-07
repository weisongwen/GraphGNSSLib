/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 * 
 * Main fucntions: add the raw GNSS data into the rosbag file by alligning the rostopic from SPAN-CPT.
 * input: rosbag, rinex file gnss raw measurements.
 * output: rosbag with raw gnss data 
 * Date: 2020/11/28
 *******************************************************/

// std inputs and outputs, fstream
#include <iostream>
#include <string>  
#include <fstream>
#include<sstream>
#include <stdlib.h>
#include <iomanip>

// math
#include <math.h>
//time 
#include <time.h>
//algorithm 
#include <algorithm>

// google eigen
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include<Eigen/Core>

// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>
// ros
#include <ros/ros.h>
/* Reference from NovAtel GNSS/INS */
#include <novatel_oem7_msgs/INSPVAX.h> // novatel_msgs/INSPVAX
#include <novatel_oem7_msgs/BESTPOS.h> // novatel_msgs/INSPVAX
#include "gnss_tools.h"
#include <nlosExclusion/GNSS_Raw_Array.h>
#include <geometry_msgs/Point32.h>
#include <stdio.h>
#include <queue>
#include <map>
#include <queue>
#include <mutex>
#include <thread>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include "../tic_toc.h"
// allign 
#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

#include <sensor_msgs/NavSatFix.h>

#include <novatel_oem7_msgs/INSPVAX.h> // novatel_msgs/INSPVAX
#include <novatel_oem7_msgs/BESTPOS.h> // novatel_msgs/INSPVAX

// rtklib
#include <stdarg.h>
#include "../../RTKLIB/src/rtklib.h"

class rosbag_generator
{
    // ros::NodeHandle nh;

    /* ros subscriber */
    ros::Publisher pub_WLSENU, pub_FGOENU, pub_global_path, pubGNSSRaw, pubStationGNSSRaw, pubDopplerVel;
    
    std::map<int, nlosExclusion::GNSS_Raw_Array> gnss_raw_map;
    std::map<int, nav_msgs::Odometry> doppler_map;
    std::map<int, nlosExclusion::GNSS_Raw_Array> station_gnss_raw_map;

    GNSS_Tools m_GNSS_Tools; // utilities  

    /* subscriber */
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> gnss_raw_array_sub;
    std::unique_ptr<message_filters::Subscriber<nav_msgs::Odometry>> doppler_sub;
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> station_gnss_raw_array_sub;
    std::unique_ptr<message_filters::TimeSynchronizer<nlosExclusion::GNSS_Raw_Array, nlosExclusion::GNSS_Raw_Array, nav_msgs::Odometry>> syncdoppler2GNSSRaw;

    /* thread lock for data safe */
    std::mutex m_gnss_raw_mux;

    int gnss_frame = 0;
    int curGPSSec = 0;

    bool hasNewData = false;

    /* thread for data processing */
    std::thread publishGNSSTopicThread;

    /* parameters to be get */
    int startGPSSec = 0;
    int endGPSSec   = 0; 

    Eigen::Matrix<double, 3,1> ENU_ref;

    ros::Subscriber span_BP_sub;

    bool finishGNSSReader = false;

public:
    rosbag_generator(ros::NodeHandle& nh)
    {      
        /* thread for factor graph optimization */
        // publishGNSSTopicThread = std::thread(&rosbag_generator::GNSSDataPublisherFun, this);
        
        /* publisher */
        // pub_WLSENU = nh.advertise<nav_msgs::Odometry>("WLSGoGPS", 100); // 

        /* subscriber of three topics  */
        gnss_raw_array_sub.reset(new message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>(nh, "/gnss_preprocessor_node/GNSSPsrCarRov1", 10000));
        doppler_sub.reset(new message_filters::Subscriber<nav_msgs::Odometry>(nh, "/gnss_preprocessor_node/GNSSDopVelRov1", 10000));
        /* GNSS measurements from station (reference end) */
        station_gnss_raw_array_sub.reset(new message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>(nh, "/gnss_preprocessor_node/GNSSPsrCarStation1", 10000)); // measurements from station
        syncdoppler2GNSSRaw.reset(new message_filters::TimeSynchronizer<nlosExclusion::GNSS_Raw_Array, nlosExclusion::GNSS_Raw_Array, nav_msgs::Odometry>(*gnss_raw_array_sub, *station_gnss_raw_array_sub, *doppler_sub, 32));
        syncdoppler2GNSSRaw->registerCallback(boost::bind(&rosbag_generator::gnssraw_doppler_msg_callback,this, _1, _2, _3));

        /* subscribe to span-cpt message */
        span_BP_sub =nh.subscribe("/novatel_data/inspvax", 500, &rosbag_generator::span_bp_callback, this);

        /* publish the path from factor graph optimization */
        pub_global_path = nh.advertise<nav_msgs::Path>("/FGOGlobalPath", 100); // 

        /* publish the raw gnss data*/
        pubGNSSRaw = nh.advertise<nlosExclusion::GNSS_Raw_Array>("/rosbag_generator_node/GNSSPsrCarRov1", 100); //

        /* publish the raw gnss data from station*/ 
        pubStationGNSSRaw = nh.advertise<nlosExclusion::GNSS_Raw_Array>("/rosbag_generator_node/GNSSPsrCarStation1", 100); // 

        /* publish the Doppler gnss data from station*/ 
        pubDopplerVel = nh.advertise<nav_msgs::Odometry>("/rosbag_generator_node/GNSSDopVelRov1", 100); // 
        
        /* reference point for ENU calculation */
        ENU_ref<< ref_lon, ref_lat, ref_alt;

        /* get parameters */
        nh.param("startGPSSec",   startGPSSec, 2);
        nh.param("endGPSSec",     endGPSSec, 2);
        // nh.param("soltype",soltype, 2);

    }

    /* check the valid epoch based on gps time span*/
    bool checkValidEpoch(double gps_sec)
    {
        if((gps_sec >= start_gps_sec) && (gps_sec <=end_gps_sec))
        {
            return true;
        }
        else return false;
    }

    /**
    * @brief span_cpt callback
    * @param span_cpt bestpos msg
    * @return void
    @ 
    */
    void span_bp_callback(const novatel_oem7_msgs::INSPVAXConstPtr& fix_msg)
    {
        int gpsSec = fix_msg->nov_header.gps_week_milliseconds;
        gpsSec = gpsSec / 1000;
        std::cout<<"GPS seconds from Rosbag" << gpsSec << "\n";
        if(finishGNSSReader)
        {
            std::map<int, nlosExclusion::GNSS_Raw_Array>::iterator GNSSiter_pr;
            GNSSiter_pr = gnss_raw_map.begin();
            int length = gnss_raw_map.size();

            /* find rover gnss data and publish */
            nlosExclusion::GNSS_Raw_Array roverGNSSRawData;
            GNSSiter_pr = gnss_raw_map.find(gpsSec);
            if(GNSSiter_pr!=gnss_raw_map.end())
            {
                roverGNSSRawData = GNSSiter_pr->second;
                pubGNSSRaw.publish(roverGNSSRawData);
            }

            /* find station gnss data and publish */
            nlosExclusion::GNSS_Raw_Array stationGNSSRawData;
            GNSSiter_pr = station_gnss_raw_map.find(gpsSec);
            if(GNSSiter_pr!=station_gnss_raw_map.end())
            {
                stationGNSSRawData = GNSSiter_pr->second;
                pubStationGNSSRaw.publish(stationGNSSRawData);
            }

            /* find rover doppler vel data and publish */
            std::map<int,nav_msgs::Odometry>::iterator gnssDopplerVelIter;;
            nav_msgs::Odometry roverDopplerVelData;
            gnssDopplerVelIter = doppler_map.find(gpsSec);
            if(gnssDopplerVelIter!=doppler_map.end())
            {
                roverDopplerVelData = gnssDopplerVelIter->second;
                pubDopplerVel.publish(roverDopplerVelData);
            }
        }
    }

    /**
   * @brief gnss raw msg and doppler msg callback
   * @param gnss raw msg and doppler msg
   * @return void
   @ 
   */
   void gnssraw_doppler_msg_callback(const nlosExclusion::GNSS_Raw_ArrayConstPtr& gnss_msg, const nlosExclusion::GNSS_Raw_ArrayConstPtr& station_gnss_msg, const nav_msgs::OdometryConstPtr& doppler_msg)
    {
        m_gnss_raw_mux.lock();

        
        if(finishGNSSReader) 
        {
            m_gnss_raw_mux.unlock();
            return;
        }

        hasNewData = true;
        gnss_frame++;
        double time0 = gnss_msg->GNSS_Raws[0].GNSS_time;
        double time1 = station_gnss_msg->GNSS_Raws[0].GNSS_time;
        double time_frame = doppler_msg->pose.pose.position.x;

        curGPSSec = gnss_msg->GNSS_Raws[0].GNSS_time;

        // std::cout<<"gnss time0 " <<time0 <<std::endl; 
        // std::cout<<"doppler time_frame " <<time_frame <<std::endl;
        // std::cout<<"station time1 " <<time1 <<std::endl;

        /* save the  */
        if(checkValidEpoch(time_frame) && m_GNSS_Tools.checkRepeating(*gnss_msg))
        {
            if(gnss_msg->GNSS_Raws.size())
            {
                doppler_map[int(time_frame)] = *doppler_msg;
                gnss_raw_map[int(time_frame)] = *gnss_msg;
                station_gnss_raw_map[int(time_frame)] = *station_gnss_msg;


                Eigen::MatrixXd eWLSSolutionECEF = m_GNSS_Tools.WeightedLeastSquare(
                                        m_GNSS_Tools.getAllPositions(*gnss_msg),
                                        m_GNSS_Tools.getAllMeasurements(*gnss_msg),
                                        *gnss_msg, "WLS");
                Eigen::Matrix<double ,3,1> WLSENU;
                WLSENU = m_GNSS_Tools.ecef2enu(ENU_ref, eWLSSolutionECEF);
                LOG(INFO) << "WLSENU -> "<< std::endl << WLSENU;
            }
        }

        if(curGPSSec>end_gps_sec)
        {
            finishGNSSReader = true;
            std::cout<< " you can play the bag file now!  " <<doppler_map.size()<<std::endl;
        }
        
        /* release the lock */
        m_gnss_raw_mux.unlock();
    }

    ~rosbag_generator()
    {
    }
   
};

int main(int argc, char **argv)
{
    FLAGS_logtostderr = 1;  // output to console
    google::InitGoogleLogging(argv[0]); // init the google logging
    google::ParseCommandLineFlags(&argc, &argv, true); // parseCommandLineFlags 
    ros::init(argc, argv, "rosbag_generator_node"); 
    ros::NodeHandle nh;
    ROS_INFO("\033[1;32m----> rosbag_generator_node Started.\033[0m"); 
    // ...
    rosbag_generator rosbag_generator_(nh);
    ros::spin();
    // while(ros::ok())
    // {
    //     ros::spinOnce();
    // }
    
    return 0;
}
