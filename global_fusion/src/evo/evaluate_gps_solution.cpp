/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
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
#include "../../RTKLIB/src/rtklib.h"

using namespace std;

/* ublox TST data 20190428 part of UrbanNav*/
// std::string sol_folder = "/home/wws/GraphGNSSLib/FGO_trajectoryllh_psr_dop_fusion.csv"; // solution from FGO  gnss_ublox_wls.csv  FGO_trajectoryllh_psr_dop_fusion
// std::string gt_sol_folder = "/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/gps_solution_TST/groundTruth_TST.csv";

// FILE* posError = fopen("/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/gps_solution_TST/error.csv", "w+");
// FILE* trajectory = fopen("/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/gps_solution_TST/trajectory.csv", "w+");


/* ublox midundao data 20200312 collected by guohao and Zhangbo*/
// std::string sol_folder = "/home/wws/GraphGNSSLib/FGO_trajectoryllh_psr_dop_fusion.csv"; // solution from FGO  gnss_ublox_wls.csv  FGO_trajectoryllh_psr_dop_fusion
// std::string gt_sol_folder = "/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/Midundao_ublox_M8T/Route1/groundTruth_MDD.csv";
// FILE* posError = fopen("/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/Midundao_ublox_M8T/Route1/error.csv", "w+");
// FILE* trajectory = fopen("/home/wws/GraphGNSSLib/src/GraphGNSSLib/global_fusion/src/rtk_estimator/Midundao_ublox_M8T/Route1/trajectory.csv", "w+");

class evaluate_gps_solution
{
    ros::NodeHandle nh;

    Eigen::Matrix<double, 3,1> ENU_ref;
    
    /* get reference solution */
    std::thread optimizationThread;

    GNSS_Tools m_GNSS_Tools; // utilities

public:
    struct Solution
    {
        double timestamp; // float timestamp
        double gps_week; // gps week 
        double gps_sec; // gps second
        double longitude; 
        double latitude;
        double altitude;
    };

    std::map<int, Solution> m_solutions_Vec;
    std::map<int, Solution> m_gt_solutions_Vec;

    std::string sol_folder;
    std::string gt_sol_folder;

    std::string trajectory_path;
    std::string error_path;

public:
    evaluate_gps_solution(ros::NodeHandle& nh)
    {
        ros::param::get("sol_folder", sol_folder);
        ros::param::get("gt_sol_folder", gt_sol_folder);

        ros::param::get("trajectory_path", trajectory_path);
        ros::param::get("error_path", error_path);


        optimizationThread = std::thread(&evaluate_gps_solution::optimization, this);

        //   ENU_ref<< 114.179000972, 22.3011535667, 0; // 20190428 data
        //   ENU_ref<<-2414266.9200, 5386768.9870, 2407460.0310;
        /* Dynamic data collected in TST*/
        // ENU_ref<< 114.179000972,22.3011535667,6.42821512092;

        /* static data collected in TST by ivan*/
        // ENU_ref<< 114.177707462258,22.2999154035483 ,4.89580645161292;

        /* Kowloon Tong data, evaluation for GNSS-GNC (dynamic, loop) */
        // ENU_ref<< 114.179526047,22.3304807923 ,11.7615366792;
        // ENULlhRef.resize(3,1);
        ENU_ref<< ref_lon, ref_lat, ref_alt;

        getSolutionfromCSV();
        getGroundTruthSolutionfromCSV();
        evaluateSol();


    }
    ~evaluate_gps_solution()
    {
        optimizationThread.detach();
    }

    void optimization()
    {
        
    }
    
    /* get solution from CSV file */
    void getSolutionfromCSV()
    {
        // while(1)
        {
            std::chrono::milliseconds dura(1000); // this thread sleep for any ms
            std::this_thread::sleep_for(dura);
            // load image list
            FILE* solFile;

            solFile = std::fopen((sol_folder).c_str() , "r");
            if(solFile == NULL){
                printf("cannot find file: solution File \n", sol_folder.c_str());
                ROS_BREAK();
            }
            double Timestamp;
            double latitude;
            double longitude;
            double altitude;
            char line[1024];
            Solution Solution_;
            map<int, Solution>::iterator iterSol;

            if(!m_solutions_Vec.size())
            {
                while ((fscanf(solFile, "%[^\n]", line)) != EOF)
                {
                fgetc(solFile);    // Reads in '\n' character and moves file
                                    // stream past delimiting character
                printf("Line = %s \n", line);
                std::stringstream ss(line); // split into three string
                vector<string> result;
                while (ss.good())
                {
                    string substr;
                    getline(ss, substr, ',');
                    result.push_back(substr);
                    std::cout << std::setprecision(17);
                }
                Solution_.gps_week =strtod((result[0]).c_str(), NULL);
                Solution_.gps_sec =strtod((result[1]).c_str(), NULL);

                Solution_.latitude = strtod((result[2]).c_str(), NULL);
                Solution_.longitude = strtod((result[3]).c_str(), NULL);
                Solution_.altitude = strtod((result[4]).c_str(), NULL);
                std::cout << std::setprecision(17);
                m_solutions_Vec[int(Solution_.gps_sec)] = Solution_;
                // std::cout<<"m_NMEAVector.size() = "<<m_NMEAVector.size()<<std::endl;
                }
            }
            std::fclose(solFile);
            // std::cout<<"m_NMEAVector.size() = "<<m_NMEAVector.size()<<std::endl;
        }
    }

    /* get ground truth solution from CSV file */
    void getGroundTruthSolutionfromCSV()
    {
        // while(1)
        {
            std::chrono::milliseconds dura(1000); // this thread sleep for any ms
            std::this_thread::sleep_for(dura);
            // load image list
            FILE* gt_solFile;
            
            gt_solFile = std::fopen((gt_sol_folder).c_str() , "r");
            if(gt_solFile == NULL){
                printf("cannot find file: solution File \n", sol_folder.c_str());
                ROS_BREAK();
            }
            double Timestamp;
            double latitude;
            double longitude;
            double altitude;
            char line[1024];
            Solution Solution_;
            map<int, Solution>::iterator iterSol;

            if(!m_gt_solutions_Vec.size())
            {
                while ((fscanf(gt_solFile, "%[^\n]", line)) != EOF)
                {
                fgetc(gt_solFile);    // Reads in '\n' character and moves file
                                    // stream past delimiting character
                printf("Line = %s \n", line);
                std::stringstream ss(line); // split into three string
                vector<string> result;
                while (ss.good())
                {
                    string substr;
                    getline(ss, substr, ',');
                    result.push_back(substr);
                    std::cout << std::setprecision(17);
                }
                Solution_.gps_week =strtod((result[0]).c_str(), NULL);
                Solution_.gps_sec =strtod((result[1]).c_str(), NULL);

                Solution_.latitude = strtod((result[2]).c_str(), NULL);
                Solution_.longitude = strtod((result[3]).c_str(), NULL);
                Solution_.altitude = strtod((result[4]).c_str(), NULL);
                std::cout << std::setprecision(17);
                m_gt_solutions_Vec[int(Solution_.gps_sec)] = Solution_;
                // std::cout<<"m_NMEAVector.size() = "<<m_NMEAVector.size()<<std::endl;
                }
            }
            std::fclose(gt_solFile);
            // std::cout<<"m_NMEAVector.size() = "<<m_NMEAVector.size()<<std::endl;
        }
    }

    double getMSE(Eigen::MatrixXd est, Eigen::MatrixXd refSat, std::string format)
    {
        double MSE = 0;
        if(format=="2D_Error")
        {
            MSE = sqrt(pow((est(0) - refSat(0)),2) + pow((est(1) - refSat(1)),2));
            return MSE;
        }
        else if(format=="3D_Error")
        {
            MSE = sqrt(pow((est(0) - refSat(0)),2) + pow((est(1) - refSat(1)),2) + pow((est(2) - refSat(2)),2));
            return MSE;
        }
        
    }

    void evaluateSol()
    {
        std::ofstream foutErrorPath(error_path, std::ios::ate);
        foutErrorPath.setf(std::ios::fixed, std::ios::floatfield);
        
        std::ofstream foutPathPath(trajectory_path, std::ios::ate);
        foutPathPath.setf(std::ios::fixed, std::ios::floatfield);

        std::map<int, Solution>::iterator sol_iter, gt_sol_iter;
        int sol_epoch_size = m_solutions_Vec.size();
        sol_iter = m_solutions_Vec.begin();
        gt_sol_iter = m_gt_solutions_Vec.begin();
        for(int i = 0; i < sol_epoch_size; i++, sol_iter++)
        {
            int gps_sec = int(sol_iter->first);
            gt_sol_iter = m_gt_solutions_Vec.find(gps_sec);
            if(gt_sol_iter != m_gt_solutions_Vec.end())
            {

                Eigen::Matrix<double, 3, 1> gt_sol_llh, gt_sol_ecef, gt_sol_enu;
                gt_sol_llh<< gt_sol_iter->second.longitude, gt_sol_iter->second.latitude, gt_sol_iter->second.altitude;
                gt_sol_ecef = m_GNSS_Tools.llh2ecef(gt_sol_llh);
                gt_sol_enu = m_GNSS_Tools.ecef2enu(ENU_ref, gt_sol_ecef);

                Eigen::Matrix<double, 3, 1> sol_llh, sol_ecef, sol_enu;
                sol_llh << sol_iter->second.longitude, sol_iter->second.latitude, sol_iter->second.altitude;
                sol_ecef = m_GNSS_Tools.llh2ecef(sol_llh);
                sol_enu = m_GNSS_Tools.ecef2enu(ENU_ref, sol_ecef);

                double error_2d = getMSE(sol_enu, gt_sol_enu, "2D_Error");
                double error_3d = getMSE(sol_enu, gt_sol_enu, "3D_Error");
                std::cout<<"get ground truth solution at apoch "<< gps_sec << "       error-> " << error_2d<<std::endl;

                
                foutErrorPath.precision(10);             
                foutErrorPath<<gps_sec<<",";
                foutErrorPath<<error_2d<<",";
                foutErrorPath<<error_3d<<",";
                foutErrorPath<<std::endl;

                
                foutPathPath.precision(0);    
                foutPathPath<<gps_sec<<",";  
                foutPathPath.precision(10);       
                foutPathPath<<gt_sol_enu(0)<<",";
                foutPathPath<<gt_sol_enu(1)<<",";
                foutPathPath<<sol_enu(0)<<",";
                foutPathPath<<sol_enu(1)<<",";
                foutPathPath<<std::endl;

            }
        } 

    }


};

int main(int argc, char **argv)
{
    FLAGS_logtostderr = 1;  // output to console
    google::InitGoogleLogging(argv[0]); // init the google logging
    google::ParseCommandLineFlags(&argc, &argv, true); // parseCommandLineFlags 
    ros::init(argc, argv, "evo");  
    // ...
    ros::NodeHandle nh;
    evaluate_gps_solution evaluate_gps_solution(nh);
    ros::spin();
    return 0;
}
