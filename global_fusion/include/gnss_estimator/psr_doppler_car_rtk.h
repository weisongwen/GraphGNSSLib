/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 *******************************************************/

#define D2R 3.1415926/180.0
#include <nlosExclusion/GNSS_Raw_Array.h>
// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>
#include "../gnss_tools.h"

#include "pseudorange_factor.h"
#include "doppler_factor.hpp"
#include "carrier_phase_factor.h"

#define enanle_Gaussian 1
#define state_size 3 // x,y,z, (clk_gps,clk_BeiDou are removed by DD)
#define ar_size 15  // size of variable of ambiguity 15
#define com_size state_size+ar_size // combine the state and the ambiguity in a array

#define enable_doppler_factor 1
#define doppler_vel_var 0.001 //0.01
#define enable_pseudorange_factor 0
#define enable_DD_factor 1 // enable DD measurement factor
#define enable_static_MM 0 // enable zero velocity motion model?

#define enable_ambiguity_resolution 1 // enable Ambiguity Resolution
#define use_fixed_cov 0 // use fixed covariance for DD measurement?

GNSS_Tools m_GNSS_Tools; // utilities

class FactorGraph{
public:
    /* continuous data stream */
    std::map<double, nlosExclusion::GNSS_Raw_Array> gnss_raw_map;
    std::map<double, nav_msgs::Odometry> doppler_map;
    std::map<double, nlosExclusion::GNSS_Raw_Array> station_gnss_raw_map;

    /* Ceres solver object */
    // ceres::Problem problem;
    // ceres::Solver::Options options;
    // ceres::Solver::Summary summary;
    ceres::LossFunction *loss_function;

    /* size of factor graph */
    int sizeOfFactorGraph = 10000000;

    /* position state array of factor graph */
    std::vector<double*> state_array;

    /* ambiguity state array of factor graph */
    std::vector<double*> ar_state_array;

    /* array save the num of ambiguity unknowns for each epoch */
    std::vector<int*> ar_state_num;

    /* gps second array */
    std::vector<int> gps_sec_array;
    
    /* position state saved in vector pair */
    std::vector<std::pair<double, Eigen::Vector3d>> Ps;

    /* ambiguity state saved in vector pair <time, <PRN, ambiguity>> */
    std::vector<std::pair<double, std::vector<std::pair<int, double>>>> AR;

    /* size of the last factor graph optimization */
    int lastFactorGraphSize = 0;

    /* fixed variance of doppler measurements */
    double var = 0.6;

    /* reference point for ENU calculation */
    Eigen::MatrixXd ENULlhRef;

    /* measurements size */
    int measSize = 0; 

    /* parameters */
    int numOfPsrFactors = 0;
    int numOfDopplerFactors = 0;
    int numOfStates =0;

    /* latest GNSS-RTK solution with LAMBDA */
    Eigen::Matrix<double, 3,1> fixedStateGNSSRTK;
    int fixed_cnt = 0;
    

public:
    /* input gnss raw (pseudorange/carrier-phase) data  */
    bool input_gnss_raw_data(nlosExclusion::GNSS_Raw_Array GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            gnss_raw_map[timestamp] = GNSS_data;
            return true;
        }
    }

    /* input Doppler data  */
    bool input_doppler_data(nav_msgs::Odometry dopplerData, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            doppler_map[timestamp] = dopplerData;
            return true;
        }
    }

    /* input GNSS data from station  */
    bool input_station_data(nlosExclusion::GNSS_Raw_Array GNSS_data, double timestamp)
    {
        if(timestamp<0) return false;
        else 
        {
            station_gnss_raw_map[timestamp] = GNSS_data;
            return true;
        }
    }

    /* input gnss doppler data */
    bool setWindowSize(int windowSize)
    {
        sizeOfFactorGraph = windowSize;
        return true;
    }

    /* clear data stream */
    bool clearDataStream()
    {
        gnss_raw_map.clear();
        doppler_map.clear();
        station_gnss_raw_map.clear();
        return true;
    }

    /* set up ceres-solver options */
    bool setupSolverOptions(ceres::Solver::Options& options)
    {
        options.use_nonmonotonic_steps = true;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.trust_region_strategy_type = ceres::TrustRegionStrategyType::DOGLEG;
        options.dogleg_type = ceres::DoglegType::SUBSPACE_DOGLEG;
        options.num_threads = 8;
        options.max_num_iterations = 258;
        return true;
    }

    /* set up Loss functions options */
    bool setupLossFunction(std::string loss)
    {
        if(loss=="Huber")
            loss_function = new ceres::HuberLoss(1.0);
        else 
        {
            loss_function = new ceres::CauchyLoss(1.0);
        }
        return true;
    }

    /* get data stream size */  
    int getDataStreamSize()
    {
        measSize = gnss_raw_map.size();
        return measSize;
    }

    bool initializeFactorGraphParas()
    {
        numOfPsrFactors = 0;
        numOfDopplerFactors = 0;
        numOfStates =0;
    }

    /* setup state size */
    bool setupPositionStateMemory()
    {
        state_array.reserve(measSize);
        ar_state_num.reserve(measSize);
        int length = measSize;
        LOG(INFO) << "length" << length << std::endl;

        for(int i = 0; i < length;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            state_array[i]  = new double[state_size]; //
            ar_state_num[i] = new int[1];
        }

        return true;
    }

    /* setup ambiguity state size 
     * this is slightly complicate, the ambiguity number is an variable
    */
    bool setARStateMemoryAndInitialize()
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_sm; // station gnss measurements map iterator
        int length = measSize;
        gps_sec_array.reserve(length);
        iter_sm = station_gnss_raw_map.begin();
        for(int k = 0;  k < length; k++,iter_sm++) // tranverse the whole station gnss measurements map
        {
            nlosExclusion::GNSS_Raw_Array st_gnss_data = (iter_sm->second);
            int sv_cnt = st_gnss_data.GNSS_Raws.size();
            double t = iter_sm->first;
            
            /* find user end gnss data with closest time */
            nlosExclusion::GNSS_Raw_Array closest_gnss_data;
            findClosestEpoch(t, gnss_raw_map, closest_gnss_data);


            gps_sec_array[k] = int(closest_gnss_data.GNSS_Raws[0].GNSS_time);

            /* get the dd measurements between
            * 1. st_gnss_data from station
            * 2. closest_gnss_data from user end
            */

            /* tranverse station gnss data at a epoch */
            int carrier_phase_index = 0;
            for(int q = 0; q < sv_cnt; q++)
            {
                double sat_id = st_gnss_data.GNSS_Raws[q].prn_satellites_index;
                /*u_master_sv: user to master 
                 *u_iSV: user to ith satellite
                 *r_master_sv: reference to master satellite
                */
                nlosExclusion::GNSS_Raw u_master_sv, u_iSV, r_master_sv;
                
                /* find the master satellite from the user end */
                if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV))
                {
                    if(u_iSV.elevation==0)
                    {
                        LOG(INFO) << "satellite with zero elevation angle---";
                        LOG(INFO) << "satellite with pseudorange---"<<u_iSV.pseudorange;
                    }
                    /* find the satellite from station gnss with the same id with master satellite */
                    findSatellitewithSameId(u_master_sv.prn_satellites_index, st_gnss_data, r_master_sv);

                    DDMeasurement DD_measurement;
                    DD_measurement.u_master_SV = u_master_sv;
                    DD_measurement.u_iSV = u_iSV;

                    DD_measurement.r_master_SV = r_master_sv;
                    DD_measurement.r_iSV = st_gnss_data.GNSS_Raws[q];

                    Eigen::Vector3d base_pose(station_x, station_y, station_z);
                    
                    if(checkCarrierPhaseConsistency(DD_measurement)) 
                    // (carrier_phase_index<10))
                    { 
                        carrier_phase_index++;
                    }
                    else
                    {
                        // LOG(INFO)<<"no carrier-phase measurement";
                    }
                }
            }
            ar_state_num[k][0] = carrier_phase_index;
            
        }

        /* allocate memory to the ambiguity state */
        ar_state_array.reserve(length);

        for(int i = 0; i < length;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            ar_state_array[i] = new double[ar_state_num[i][0]]; //
            for(int j = 0; j < ar_state_num[i][0]; j++)
            {
                ar_state_array[i][j] = 0;
            }
        }

        return true;
    }


    /* setup const ambiguity state size 
    */
    bool setConstARStateMemoryAndInitialize()
    {
        /* allocate memory to the ambiguity state */
        int length = measSize;
        gps_sec_array.reserve(length);
        ar_state_array.reserve(length);

        for(int i = 0; i < length;i++)
        {
            /* ECEF_x, ECEF_y, ECEF_z */
            ar_state_array[i] = new double[ar_size]; //
            for(int j = 0; j < ar_size; j++)
            {
                ar_state_array[i][j] = 0;
            }
        }

        return true;
    }

    /* initialize the previous optimzied states */
    bool initializeOldGraph()
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr;
        iter_pr = gnss_raw_map.begin();
        int length = measSize;

        /* tranverse the stateArray */
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            nlosExclusion::GNSS_Raw_Array gnss_data = (iter_pr->second);
            double time = gnss_data.GNSS_Raws[0].GNSS_time;

            /* if the state vector is NOT empty, update */
                for(int i = 0; i < Ps.size(); i++)
                {
                    if(time == Ps[i].first)
                    {
                        state_array[m][0] = Ps[i].second[0];
                        state_array[m][1] = Ps[i].second[1];
                        state_array[m][2] = Ps[i].second[2];
                        // state_array[m][3] = Clocks[i].second[0];
                        // state_array[m][4] = Clocks[i].second[1];
                    }
                } 
            
        }
        return true;
    }

    /* initialize the newly added state using WLS*/
    bool initializeNewlyAddedGraph()
    {
        int length = measSize;
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter;
        iter = gnss_raw_map.begin();
        for(int i = 0; i <length; i++,iter++)
        {
            if(i>=(lastFactorGraphSize-1))
            {
                nlosExclusion::GNSS_Raw_Array gnss_data = (iter->second);
                Eigen::MatrixXd eWLSSolutionECEF = m_GNSS_Tools.WeightedLeastSquare(
                                    m_GNSS_Tools.getAllPositions(gnss_data),
                                    m_GNSS_Tools.getAllMeasurements(gnss_data),
                                    gnss_data, "WLS");
                state_array[i][0] = eWLSSolutionECEF(0);
                state_array[i][1] = eWLSSolutionECEF(1); 
                state_array[i][2] = eWLSSolutionECEF(2);
            }
        }
        return true;
    }

    /* add parameter blocks */
    bool addParameterBlocksToGraph(ceres::Problem& problem)
    {
        int length = measSize;
        for(int i = 0; i <length; i++)
        {
            /* add parameter block for position state (ECEF_x, ECEF_y, ECEF_z) */
            problem.AddParameterBlock(state_array[i], state_size);

            /* add parameter block for ambiguity state */
            problem.AddParameterBlock(ar_state_array[i], ar_size);
        }
        return true;
    }

    /* fix the first parameter block */
    bool fixFirstState(bool flag, ceres::Problem& problem)
    {
        /* fixed the first state only after first optimization */
        if(flag && Ps.size()>10)
        {
            problem.SetParameterBlockConstant(state_array[0]);
            problem.SetParameterBlockConstant(ar_state_array[0]);
        }    
        return true;
    }

    /* add Doppler FACTORS */
    bool addDopplerFactors(ceres::Problem& problem)
    {
        /* process doppler measurements */
        std::map<double, nav_msgs::Odometry>::iterator iterdopp, iterdoppNext;
        int i = 0;
        for(iterdopp = doppler_map.begin(); iterdopp != doppler_map.end();iterdopp++, i++)
        {
            /* add doppler measurements */
            iterdoppNext = iterdopp;
            iterdoppNext ++;
            if(iterdoppNext != doppler_map.end())
            {
                double delta_t = iterdoppNext->first - iterdopp->first;
                double v_x_i = iterdopp->second.twist.twist.linear.x;
                double v_y_i = iterdopp->second.twist.twist.linear.y;
                double v_z_i = iterdopp->second.twist.twist.linear.z;
                #if enable_static_MM
                v_x_i = 0.0;
                v_y_i = 0.0;
                v_z_i = 0.0;
                #endif
                double dop_var_scal = 0.06;
                double var_x = dop_var_scal * sqrt(iterdopp->second.twist.covariance[0]);
                double var_y = dop_var_scal * sqrt(iterdopp->second.twist.covariance[1]);
                double var_z = dop_var_scal * sqrt(iterdopp->second.twist.covariance[2]);
                double var = doppler_vel_var; // 
                Eigen::Vector3d var_vec(var,var,var);
                ceres::CostFunction* doppler_function = new ceres::AutoDiffCostFunction<dopplerFactor, 3 
                                                        , state_size,state_size>(new 
                                                        dopplerFactor(v_x_i, v_y_i, v_z_i, delta_t,  var_vec));
                problem.AddResidualBlock(doppler_function, loss_function, state_array[i],state_array[i+1]);
            }
        }
        return true;
    }


    /* add double-differenced pseudorange/Carrier-phase FACTORS */
    bool addDDPsrCarFactors(ceres::Problem& problem)
    {
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_sm; // station gnss measurements map iterator
        int length = measSize;
        iter_sm = station_gnss_raw_map.begin();
        for(int k = 0;  k < length; k++,iter_sm++) // tranverse the whole station gnss measurements map
        {
            nlosExclusion::GNSS_Raw_Array st_gnss_data = (iter_sm->second);
            int sv_cnt = st_gnss_data.GNSS_Raws.size();
            double t = iter_sm->first;
            
            /* find user end gnss data with closest time */
            nlosExclusion::GNSS_Raw_Array closest_gnss_data;
            findClosestEpoch(t, gnss_raw_map, closest_gnss_data);

            gps_sec_array[k] = int(closest_gnss_data.GNSS_Raws[0].GNSS_time);

            /* get the dd measurements between
            * 1. st_gnss_data from station
            * 2. closest_gnss_data from user end
            */

            /* tranverse station gnss data at a epoch */
            int carrier_phase_index = 0;
            for(int q = 0; q < sv_cnt; q++)
            {
                double sat_id = st_gnss_data.GNSS_Raws[q].prn_satellites_index;
                /*u_master_sv: user to master 
                 *u_iSV: user to ith satellite
                 *r_master_sv: reference to master satellite
                */
                nlosExclusion::GNSS_Raw u_master_sv, u_iSV, r_master_sv;
                
                /* find the master satellite from the user end */
                if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV))
                {
                    if(u_iSV.elevation==0)
                    {
                        LOG(INFO) << "satellite with zero elevation angle---";
                        LOG(INFO) << "satellite with pseudorange---"<<u_iSV.pseudorange;
                    }
                    /* find the satellite from station gnss iwth the same id with master satellite */
                    findSatellitewithSameId(u_master_sv.prn_satellites_index, st_gnss_data, r_master_sv);

                    DDMeasurement DD_measurement;
                    DD_measurement.u_master_SV = u_master_sv;
                    DD_measurement.u_iSV = u_iSV;

                    DD_measurement.r_master_SV = r_master_sv;
                    DD_measurement.r_iSV = st_gnss_data.GNSS_Raws[q];

                    Eigen::Vector3d base_pose(station_x, station_y, station_z);
                    ceres::CostFunction* dd_pr_function = new ceres::AutoDiffCostFunction<DDpseudorangeVSVFactor, 1 
                                                        , state_size>(new 
                                                        DDpseudorangeVSVFactor(DD_measurement, base_pose));
                    auto ID = problem.AddResidualBlock(dd_pr_function, loss_function, state_array[k]);

                    if(checkCarrierPhaseConsistency(DD_measurement)) 
                    // (carrier_phase_index<10))
                    { 
                        ceres::CostFunction* dd_cp_function = new ceres::AutoDiffCostFunction<DDCarrierPhaseFactor_DDBias, 1 
                                                        , state_size,ar_size>(new 
                                                        DDCarrierPhaseFactor_DDBias(DD_measurement, base_pose,carrier_phase_index));
                        auto IDs = problem.AddResidualBlock(dd_cp_function, loss_function, state_array[k],ar_state_array[k]); 

                        carrier_phase_index++;
                        // LOG(INFO)<<"add double-difference carrier phase factor->"<<carrier_phase_index;
                    }
                    else
                    {
                        // LOG(INFO)<<"no carrier-phase measurement";
                    }
                }
            }

            // LOG(INFO) << "sv_cnt/carrier_phase_index-> "<<sv_cnt<<"  "<<carrier_phase_index;
            ar_state_num[k][0] = carrier_phase_index;
            
        }
    }


    /* solve the factor graph */
    bool solveFactorGraphFloatSolution(ceres::Problem& problem,ceres::Solver::Options& options, ceres::Solver::Summary& summary)
    {
        /* solve the problem*/
        ceres::Solve(options, &problem, &summary);
        return true;
    }

    /* solve the ambiguity resolution */
    bool solveAmbiguityResolutionFixedSolution()
    {
        /* construct the ambiguity resolution problem*/
        int length = measSize;
        /* get the Jacobian matrix for covariance propogation */
        Eigen::MatrixXd jacobian_matrix; 
        Eigen::MatrixXd weighting_matrix;

        /* number of ambiguity variables */
        int AR_cnt = 0;
        for(int i = 0; i < ar_size; i++)
        {
            if(fabs(ar_state_array[length-1][i])>0)
            AR_cnt++;
        }
        LOG(INFO) << "AR_cnt-> " << AR_cnt;
        LOG(INFO) << "ar_state_num[k][0]-> "<< ar_state_num[length-1][0];
        AR_cnt = ar_state_num[length-1][0];
        jacobian_matrix.resize(40, 3 + AR_cnt);
        jacobian_matrix.setIdentity();
        weighting_matrix.resize(40, 40);
        weighting_matrix.setIdentity();

        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_cov; // station gnss measurements map iterator
        iter_cov = station_gnss_raw_map.end();
        iter_cov--;

        nlosExclusion::GNSS_Raw_Array st_gnss_data = (iter_cov->second);
        int sv_cnt = st_gnss_data.GNSS_Raws.size();
        double t = iter_cov->first;

        /* find user end gnss data with closest time */
        nlosExclusion::GNSS_Raw_Array closest_gnss_data;
        findClosestEpoch(t, gnss_raw_map, closest_gnss_data);

        /* get the dd measurements between
        * 1. st_gnss_data from station
        * 2. closest_gnss_data from user end
        */

        /* tranverse station gnss map */
        int carrier_phase_index = 0;
        int jac_row = 0;
        Eigen::Vector3d u_pose(state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]);
        for(int q = 0; q < sv_cnt; q++)
        {
            double sat_id = st_gnss_data.GNSS_Raws[q].prn_satellites_index;
            nlosExclusion::GNSS_Raw u_master_sv, u_iSV, r_master_sv;
            
            /* find the master satellite from the user end */
            if(findMasterSatellite(sat_id, closest_gnss_data, u_master_sv, u_iSV))
            {
                /* find the satellite from station gnss iwth the same id with master satellite */
                findSatellitewithSameId(u_master_sv.prn_satellites_index, st_gnss_data, r_master_sv);

                DDMeasurement DD_measurement;
                DD_measurement.u_master_SV = u_master_sv;
                DD_measurement.u_iSV = u_iSV;

                DD_measurement.r_master_SV = r_master_sv;
                DD_measurement.r_iSV = st_gnss_data.GNSS_Raws[q];

                Eigen::Vector3d base_pose(station_x, station_y, station_z);

                /* get the row for jocabian of pr DD*/
                m_GNSS_Tools.getPrDDJacobian(u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, weighting_matrix);
                
                jac_row++;
                if(checkCarrierPhaseConsistency(DD_measurement) && (carrier_phase_index<AR_cnt)) 
                // (carrier_phase_index<10))
                {
                        /* get the row for jocabian of cp DD*/
                        m_GNSS_Tools.getCpDDJacobian(u_pose, base_pose, DD_measurement, jacobian_matrix,jac_row, carrier_phase_index, weighting_matrix);

                    carrier_phase_index++;
                    // LOG(INFO)<<"add double-difference carrier phase factor";
                    jac_row++;
                }
                else
                {
                    // LOG(INFO)<<"no carrier-phase measurement";
                }
            }
        }

        Eigen::MatrixXd jacobian_matrix_;
        Eigen::MatrixXd weighting_matrix_;
        weighting_matrix_.resize(jac_row,jac_row);
        jacobian_matrix_.resize(jac_row, jacobian_matrix.cols());
        for(int i = 0; i < jac_row; i++)
            for(int j = 0; j < jacobian_matrix.cols(); j++)
            {
                jacobian_matrix_(i,j) = jacobian_matrix(i,j);
            }

        for(int i = 0; i <jac_row; i++)
            for(int j = 0; j <jac_row; j++)
            {
                weighting_matrix_(i,j) = weighting_matrix(i,j);
            }

        Eigen::MatrixXd jacobian_mat;
        jacobian_mat = (jacobian_matrix_.transpose() * weighting_matrix_ *jacobian_matrix_).inverse();
        // std::cout<<"jacobian_mat-> \n"<<std::setprecision(4)<<jacobian_mat<<std::endl;
        // std::cout<<"weighting_matrix_-> \n"<<std::setprecision(4)<<weighting_matrix_<<std::endl;
        Eigen::MatrixXd Qa, Qb; // Qa: covariance of phase bias
        Qa.resize(AR_cnt, AR_cnt);
        Qb.resize(3, 3);
        for(int i = 3; i < AR_cnt+3; i++)
            for(int j = 3; j < AR_cnt+3; j++)
            {
                Qa(i-3,j-3) = jacobian_mat(i,j);
            }
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            {
                Qb(i,j) = jacobian_mat(i,j);
            }
        // std::cout<<"Qa-> \n"<<std::setprecision(4)<<Qa<<std::endl;

        double *a;
        double *Q;
        double *F;
        double S[2];
        double n  = AR_cnt;
        double m = 2;
        F = mat(n,2); Q = mat(n,n); a= mat(n,1);

        for(int i = 0; i < n; i++)
        {
            a[i] = ar_state_array[length-1][i];
        }
        std::vector<double> cov_vector;
        // fi_cov_ar = fi_cov_ar.inverse();
        for(int i = 0; i < n; i++) // cols
            for(int j = 0; j < n; j++) // rows
            {
                cov_vector.push_back(Qa(j,i));
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
            
        for(int i=0; i<cov_vector.size();i++)
        {
            Q[i] = cov_vector[i];
        }

        lambda(n, m, a, Q, F, S);
        Eigen::MatrixXd Q_ba, Q_ab;
        Q_ba.resize(3,AR_cnt);
        Q_ab.resize(AR_cnt,3);
        for(int i = 0; i < 3; i++) // rows
            for(int j = 3; j < (3+ AR_cnt); j++) // cols
            {
                Q_ba(i,j-3) = jacobian_mat(i,j);
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
        for(int i = 3; i < (3 + AR_cnt); i++) // rows
            for(int j = 0; j < 3; j++) // cols
            {
                Q_ab(i-3,j) = jacobian_mat(i,j);
                // LOG(INFO) << "Q[i+j]-> \n" << Q[i+j];
                // if(i=i) 
            }
        Eigen::MatrixXd integer_b, float_b;
        integer_b.resize(AR_cnt, 1);
        float_b.resize(AR_cnt, 1);
        for(int i=0; i<AR_cnt; i++)
        {
            float_b(i) = ar_state_array[length-1][i];
            integer_b(i) = F[i];
        }
        // std::cout<<"Q_ba-> \n"<<std::setprecision(4)<<Q_ba<<std::endl;
        Eigen::Matrix<double, 3,1> a_state;
        a_state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2]; 
        a_state = a_state + Q_ba * Qa.inverse() * (float_b - integer_b);
        // a_state = a_state - Q_ab * Qb.inverse() * (float_b - integer_b);

        Eigen::Matrix<double ,3,1> ENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], state_array[length-1][1], state_array[length-1][2];
        double fix_flag = 0;
        if((S[1]/S[0])> 3)
        {
            state = a_state; // 
            fix_flag = 1;
            // LOG(INFO) << "<<<<<Fixed!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>";
            fixed_cnt++;
        }
        fixedStateGNSSRTK = state;

        return true;
    }

    /* save graph state to vector for next solving */
    bool saveGraphStateToVector()
    {
        /* save the size of current factor graph */
        lastFactorGraphSize = measSize;

        /* get time from data stream */
        std::map<double, nlosExclusion::GNSS_Raw_Array>::iterator iter_pr;
        iter_pr = gnss_raw_map.begin();
        int length = measSize;

        /* tranverse the stateArray */
        for(int m = 0;  m < length; m++,iter_pr++) // 
        {
            nlosExclusion::GNSS_Raw_Array gnss_data = (iter_pr->second);
            double time = gnss_data.GNSS_Raws[0].GNSS_time;
            double prn = gnss_data.GNSS_Raws[0].prn_satellites_index;

            /* if the state vector is empty, override */
            if(Ps.size()==0)
            {
                Ps.push_back(std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2])));
                // Clocks.push_back(std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4])));
            }
            /* if the state vector is NOT empty, update */
            else
            {
                bool findTimeKey = false;
                for(int i = 0; i < Ps.size(); i++)
                {
                    if(time == Ps[i].first)
                    {
                        Ps[i] = std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2]));
                        // Clocks[i] = std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4]));
                        findTimeKey = true;
                    }
                }
                /* new time frame, add to state vector*/
                if(findTimeKey==false)
                {
                    Ps.push_back(std::make_pair(time, Eigen::Vector3d(state_array[m][0],state_array[m][1],state_array[m][2])));
                    // Clocks.push_back(std::make_pair(time, Eigen::Vector2d(state_array[m][3],state_array[m][4])));
                } 
            }
            
        }
        return true;
    }
    
    /* set up the reference point for ENU calculation */
    bool setupReferencePoint()
    {
        /* reference point for ENU calculation */
        ENULlhRef.resize(3,1);
        ENULlhRef<< ref_lon, ref_lat, ref_alt;
        return true;
    }

    /* get the latest state in ENU */
    Eigen::Matrix<double ,3,1> getLatestStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);    
        return FGOENU;
    }

     /* print the latest state in ENU */
    bool printLatestFloatStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        state<< state_array[length-1][0], 
                state_array[length-1][1], 
                state_array[length-1][2];
        FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);    
        std::cout << "FGOENU-> "<< FGOENU<< std::endl;  
        return true;
    }

    /* print the latest state in ENU */
    bool printLatestFixedStateENU()
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> fixedFGOENU;
        fixedFGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, fixedStateGNSSRTK);    
        std::cout << "fixedFGOENU-> "<< fixedFGOENU<< std::endl;  
        return true;
    }

    /* get the path of FGO in ENU */
    nav_msgs::Path getPathENU(nav_msgs::Path& fgo_path)
    {
        int length = measSize;
        Eigen::Matrix<double ,3,1> FGOENU;
        Eigen::Matrix<double, 3,1> state;
        fgo_path.poses.clear();
        fgo_path.header.frame_id = "map";
        for(int i = 0; i < length;i++)
        {
            state<< state_array[i][0], 
                    state_array[i][1], 
                    state_array[i][2];
            FGOENU = m_GNSS_Tools.ecef2enu(ENULlhRef, state);  
            geometry_msgs::PoseStamped pose_stamped;
            pose_stamped.header.stamp = ros::Time::now();
            pose_stamped.header.frame_id = "map";
            pose_stamped.pose.position.x = FGOENU(0);
            pose_stamped.pose.position.y = FGOENU(1);
            pose_stamped.pose.position.z = 10;
            fgo_path.poses.push_back(pose_stamped);
            // std::cout << "pose_stamped- FGO-> "<< std::endl<< pose_stamped;
        }
              
        return fgo_path;
    }

   /**
   * @brief maintain sliding window slidingWindowSize
   * @param gnss raw msg and doppler msg
   * @return void
   @ 
   */

    void removeStatesOutsideSlidingWindow()
    {
        int numElementsToRemove = 0;
        
        /* sliding window gnss raw pseudorange*/
        numElementsToRemove = gnss_raw_map.size() - sizeOfFactorGraph;
        if(numElementsToRemove<=0) return;
        auto i = gnss_raw_map.begin();
        while (i != gnss_raw_map.end() && numElementsToRemove > 0)
        {
            i = gnss_raw_map.erase(i);
            --numElementsToRemove;
        }

        /* sliding window station gnss raw pseudorange*/
        numElementsToRemove = station_gnss_raw_map.size() - sizeOfFactorGraph;
        if(numElementsToRemove<=0) return;
        auto s = station_gnss_raw_map.begin();
        while (s != station_gnss_raw_map.end() && numElementsToRemove > 0)
        {
            s = station_gnss_raw_map.erase(s);
            --numElementsToRemove;
        }

        /* sliding window gnss raw doppler */
        numElementsToRemove = doppler_map.size() - sizeOfFactorGraph;
        if(numElementsToRemove<=0) return;
        auto j = doppler_map.begin();
        while (j != doppler_map.end() && numElementsToRemove > 0)
        {
            j = doppler_map.erase(j);
            --numElementsToRemove;
        }
    }


    /* free memory */
    bool freeStateMemory()
    {
        int length = measSize;
        for(int i = 0; i < length;i++)
        {
            free(state_array[i]);
        }     
        return true;
    }

};


/* check the valid epoch based on gps time span*/
bool checkValidEpoch(double gps_sec)
{
    if((gps_sec >= start_gps_sec) && (gps_sec <=end_gps_sec))
    {
        return true;
    }
    else return false;
}