/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 * 
 * Main fucntions: pseudorange/Doppler/Carrier-phase fusion using factor graph optimization 
 * input: pseudorange (psr_), Doppler velocity (doppler_) and carrier-phase (car_) from GPS/BeiDou.
 * output: position of the GNSS receiver, ambiguity of carrier-phase 
 * SPECIFICALLY, this code is an example usage of the DynamicAutoDiffCostFunction, the size of the parameters can be a variable when we use the DynamicAutoDiffCostFunction. 
 * CURRENTLY, the Qaa for LAMBDA algorithm is derived by the covariance estimation of Ceres-solver.
 * Date: 2020/12/06
 * Currently, the covariance estimation step is time consuming
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

// rtklib
#include <stdarg.h>
#include "../../RTKLIB/src/rtklib.h"

 /* factor graph related head file */
#include "../../include/gnss_estimator/psr_doppler_car_rtk_dyna.h"

class psr_doppler_car_rtk
{
    ros::NodeHandle nh;

    /* ros subscriber */
    ros::Publisher pub_WLSENU, pub_FGOENU, pub_global_path, pub_fgo_llh;
    std::map<double, nlosExclusion::GNSS_Raw_Array> gnss_raw_map;
    std::map<double, nav_msgs::Odometry> doppler_map;
    std::map<double, nlosExclusion::GNSS_Raw_Array> station_gnss_raw_map;

    GNSS_Tools m_GNSS_Tools; // utilities  

    /* subscriber */
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> gnss_raw_array_sub;
    std::unique_ptr<message_filters::Subscriber<nav_msgs::Odometry>> doppler_sub;
    std::unique_ptr<message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>> station_gnss_raw_array_sub;
    std::unique_ptr<message_filters::TimeSynchronizer<nlosExclusion::GNSS_Raw_Array, nlosExclusion::GNSS_Raw_Array, nav_msgs::Odometry>> syncdoppler2GNSSRaw;

    /* thread lock for data safe */
    std::mutex m_gnss_raw_mux;

    /* thread for data processing */
    std::thread optimizationThread;

    int gnss_frame = 0;
    Eigen::Matrix<double, 3,1> ENU_ref;
    int slidingWindowSize = 10000; // ? epoch measurements 150 100000
    bool hasNewData = false;

    /* latest state in ENU */
    Eigen::Matrix<double ,3,1> FGOENULatest;

    /* path in ENU */
    nav_msgs::Path fgo_path;

private:
    // std::unique_ptr<factor_graph> factor_graph_ptr_; // factor graph ptr
    FactorGraph factor_graph;
    
    bool initializeFactorGraph(ceres::Solver::Options& options) 
    {
        /* initialize the factor graph size */
        factor_graph.setWindowSize(slidingWindowSize);

        /* set up ceres-solver options */
        factor_graph.setupSolverOptions(options);

        /* set up loss functions (Huber, Cauchy)*/
        factor_graph.setupLossFunction("Huber");
    }


public:
    psr_doppler_car_rtk()
    {      
        /* thread for factor graph optimization */
        optimizationThread = std::thread(&psr_doppler_car_rtk::solveOptimization, this);
        
        /* publisher */
        pub_WLSENU = nh.advertise<nav_msgs::Odometry>("WLSGoGPS", 100); // 
        pub_FGOENU = nh.advertise<nav_msgs::Odometry>("FGO", 100); //  
        pub_fgo_llh = nh.advertise<sensor_msgs::NavSatFix>("fgo_llh", 100);

        /* subscriber of three topics  */
        gnss_raw_array_sub.reset(new message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>(nh, "/gnss_preprocessor_node/GNSSPsrCarRov1", 10000));
        doppler_sub.reset(new message_filters::Subscriber<nav_msgs::Odometry>(nh, "/gnss_preprocessor_node/GNSSDopVelRov1", 10000));
        /* GNSS measurements from station (reference end) */
        station_gnss_raw_array_sub.reset(new message_filters::Subscriber<nlosExclusion::GNSS_Raw_Array>(nh, "/gnss_preprocessor_node/GNSSPsrCarStation1", 256)); // measurements from station
        syncdoppler2GNSSRaw.reset(new message_filters::TimeSynchronizer<nlosExclusion::GNSS_Raw_Array, nlosExclusion::GNSS_Raw_Array, nav_msgs::Odometry>(*gnss_raw_array_sub, *station_gnss_raw_array_sub, *doppler_sub, 32));
        syncdoppler2GNSSRaw->registerCallback(boost::bind(&psr_doppler_car_rtk::gnssraw_doppler_msg_callback,this, _1, _2, _3));

        /* publish the path from factor graph optimization */
        pub_global_path = nh.advertise<nav_msgs::Path>("/FGOGlobalPath", 100); // 
        
        /* reference point for ENU calculation */
        ENU_ref<< ref_lon, ref_lat, ref_alt;

    }


    /**
	 * @brief perform factor graph optimization
	 * @param none
	 * @return none
	 */
    void solveOptimization()
    {
        /* clear data stream */
        factor_graph.clearDataStream();
        
        while(1)
        {
            m_gnss_raw_mux.lock();

            if(factor_graph.getDataStreamSize()>3 && hasNewData)
            {
                /* define the problem */
                ceres::Problem problem;
                ceres::Solver::Options options;
                ceres::Solver::Summary summary;

                /* start clock for factor graph optimization */
                TicToc OptTime;

                /* initialize factor graph */
                initializeFactorGraph(options);

                /* get data stream size */
                factor_graph.getDataStreamSize();

                /* setup memory for position state and ambiguity state */
                factor_graph.setStateMemoryAndInitialize();

                /* initialize factor graph parameters */
                factor_graph.initializeFactorGraphParas();

                /* initialize the previous optimzied states */
                factor_graph.initializeOldGraph();

                /* initialize the newly added states */
                factor_graph.initializeNewlyAddedGraph();

                /* fix the first parameter block */
                factor_graph.fixFirstState(false, problem);

                /* add Doppler FACTORS to factor graph */
                factor_graph.addDopplerFactors(problem);

                /* add double-differenced pseudorange/Carrier-phase FACTORS */
                factor_graph.addDDPsrCarFactors(problem);

                /* solve the factor graph (float solution) */
                factor_graph.solveFactorGraphFloatSolution(problem, options, summary);

                std::cout << "float solution OptTime-> "<< OptTime.toc()<< std::endl;

                /* save graph state to variables */
                factor_graph.saveGraphStateToVector();

                /* get covariance matrix of latest epoch via Ceres-solver*/
                factor_graph.getCovarianceMatrixOfLatestEpoch(problem);

                /* solve the ambiguity resolution of current epoch */
                // factor_graph.solveAmbiguityResolutionViaCovPropogation();

                /* solve the ambiguity resolution of current epoch */
                factor_graph.solveAmbiguityResolutionViaCovCeresSolver();

                /* set reference point for ENU calculation */
                factor_graph.setupReferencePoint();

                /* get the path in factor graph */
                FGOENULatest = factor_graph.getLatestStateENU();

                /* print the lastest float state in factor graph */
                factor_graph.printLatestFloatStateENU();

                /* print the lastest fixed state in factor graph */
                factor_graph.printLatestFixedStateENU();

                /* publish the path from FGO */
                fgo_path = factor_graph.getPathENU(fgo_path);
                pub_global_path.publish(fgo_path);

                /* remove the data outside sliding window */
                factor_graph.removeStatesOutsideSlidingWindow();

                hasNewData = false;

                /** */
                std::cout << "final OptTime-> "<< OptTime.toc()<< std::endl;
            }
            m_gnss_raw_mux.unlock();
            std::chrono::milliseconds dura(10); // this thread sleep for 10 ms
            std::this_thread::sleep_for(dura);
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
        hasNewData = true;
        gnss_frame++;
        double time0 = gnss_msg->GNSS_Raws[0].GNSS_time;
        double time1 = station_gnss_msg->GNSS_Raws[0].GNSS_time;
        double time_frame = doppler_msg->pose.pose.position.x;

        // std::cout<<"gnss time0 " <<time0 <<std::endl; 
        // std::cout<<"doppler time_frame " <<time_frame <<std::endl;
        // std::cout<<"station time1 " <<time1 <<std::endl;

        /* save the  */
        if(checkValidEpoch(time_frame) && m_GNSS_Tools.checkRepeating(*gnss_msg))
        {
            if(gnss_msg->GNSS_Raws.size())
            {
                doppler_map[time_frame] = *doppler_msg;
                gnss_raw_map[time_frame] = *gnss_msg;
                station_gnss_raw_map[time_frame] = *station_gnss_msg;

                factor_graph.input_gnss_raw_data(*gnss_msg, time_frame);
                factor_graph.input_station_data(*station_gnss_msg, time_frame);
                factor_graph.input_doppler_data(*doppler_msg, time_frame);

                Eigen::MatrixXd eWLSSolutionECEF = m_GNSS_Tools.WeightedLeastSquare(
                                        m_GNSS_Tools.getAllPositions(*gnss_msg),
                                        m_GNSS_Tools.getAllMeasurements(*gnss_msg),
                                        *gnss_msg, "WLS");
                Eigen::Matrix<double ,3,1> WLSENU;
                WLSENU = m_GNSS_Tools.ecef2enu(ENU_ref, eWLSSolutionECEF);
                LOG(INFO) << "WLSENU -> "<< std::endl << WLSENU;

                nav_msgs::Odometry odometry;
                odometry.header.frame_id = "map";
                odometry.child_frame_id = "map";
                odometry.pose.pose.position.x = WLSENU(0);
                odometry.pose.pose.position.y = WLSENU(1);
                odometry.pose.pose.position.z = WLSENU(2);
                pub_WLSENU.publish(odometry);
            }
        }
        
        /* release the lock */
        m_gnss_raw_mux.unlock();
    }

    ~psr_doppler_car_rtk()
    {
    }
   
};

int main(int argc, char **argv)
{
    FLAGS_logtostderr = 1;  // output to console
    google::InitGoogleLogging(argv[0]); // init the google logging
    google::ParseCommandLineFlags(&argc, &argv, true); // parseCommandLineFlags 
    ros::init(argc, argv, "psr_doppler_car_rtk_node"); 
    ROS_INFO("\033[1;32m----> psr_doppler_car_rtk_dyna_node Started.\033[0m"); 
    // ...
    psr_doppler_car_rtk psr_doppler_car_rtk;
    ros::spin();
    return 0;
}
