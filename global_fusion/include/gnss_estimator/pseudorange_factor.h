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

/* pseudorange factor*/
struct pseudorangeFactor
{
    pseudorangeFactor(std::string sat_sys, double s_g_x, double s_g_y, double s_g_z, double pseudorange, double var)
                :sat_sys(sat_sys),s_g_x(s_g_x), s_g_y(s_g_y), s_g_z(s_g_z), pseudorange(pseudorange),var(var){}

    template <typename T>
    bool operator()(const T* state, T* residuals) const
    {
        T est_pseudorange; 
        T delta_x = pow((state[0] - s_g_x),2);
        T delta_y = pow((state[1] - s_g_y),2);
        T delta_z = pow((state[2] - s_g_z),2);

        est_pseudorange = sqrt(delta_x+ delta_y + delta_z);

        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;
        est_pseudorange = est_pseudorange + OMGE_ * (s_g_x*state[1]-s_g_y*state[0])/CLIGHT_;
        
        if(sat_sys == "GPS") 
        {
            est_pseudorange = est_pseudorange - state[3];
        }
        
        else if(sat_sys == "BeiDou") 
        {
            est_pseudorange = est_pseudorange - state[4];
        }

        residuals[0] = (est_pseudorange - T(pseudorange)) / T(var);
        // std::cout << "residuals[0]-> "<< residuals[0]<<std::endl;

        return true;
    }

    double s_g_x, s_g_y, s_g_z, pseudorange, var;
    std::string sat_sys; // satellite system

};


/* DDpseudorangeFactor factor :: different satellite system different sat pose
    * M: master satellite
    * i: satellite to be double-difference based on master satellite
    * r: reference satellite
    * u: use end satellite (GNSS receiver)
*/
struct DDpseudorangeVSVFactor
{
    DDpseudorangeVSVFactor(DDMeasurement dd_measurement, Eigen::Vector3d base_pos)
                :dd_measurement(dd_measurement),base_pos(base_pos){}

    template <typename T>
    bool operator()(const T* state, T* residuals) const
    {

        /**/
        Eigen::Vector3d pose_r = base_pos;
        GNSS_Tools m_GNSS_Tools1; // utilities
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;


        /* satellite position*/
        Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
        double var_u2m = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.u_master_SV);

        Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
        double var_u2i = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.u_iSV);
        // LOG(INFO)<<"var_u2i-> " <<var_u2i;

        Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
        double var_r2m = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.r_master_SV);
        // LOG(INFO)<<"var_r2m-> " <<var_r2m;

        Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
        double var_r2i = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.r_iSV);
        // LOG(INFO)<<"var_r2i-> " <<var_r2i;
        
        T est_p_r2m = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_m));
        // est_p_r2m = est_p_r2m + OMGE_ * (rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT_;
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;
        
        T est_p_r2i = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_i)); 
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        T est_p_u2m;
        T delta_x1 = pow((state[0] - u_pose_m(0)),2);
        T delta_y1 = pow((state[1] - u_pose_m(1)),2);
        T delta_z1 = pow((state[2] - u_pose_m(2)),2);
        est_p_u2m = sqrt(delta_x1+ delta_y1 + delta_z1);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*state[1]-u_pose_m(1)*state[0])/CLIGHT_;

        T est_p_u2i;
        T delta_x2 = pow((state[0] - u_pose_i(0)),2);
        T delta_y2 = pow((state[1] - u_pose_i(1)),2);
        T delta_z2 = pow((state[2] - u_pose_i(2)),2);
        est_p_u2i = sqrt(delta_x2+ delta_y2 + delta_z2);

        /* relatistic effects*/
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*state[1]-u_pose_i(1)*state[0])/CLIGHT_;
        
        /* expected DD measurments */
        T est_DD_pr = (est_p_u2i - est_p_r2i) - (est_p_u2m - est_p_r2m);

        T p_r2m = T(dd_measurement.r_master_SV.raw_pseudorange);
        T p_r2i = T(dd_measurement.r_iSV.raw_pseudorange);
        T p_u2m = T(dd_measurement.u_master_SV.raw_pseudorange);
        T p_u2i = T(dd_measurement.u_iSV.raw_pseudorange);
        
        /*observation of DD measurement*/
        T DD_pr = (p_u2i - p_r2i) - (p_u2m - p_r2m);

        // double var = 1000;
        double var = var_u2m + var_u2i + var_r2m + var_r2i;
        var = var / 4.0;
        // residuals[0] = (est_DD_pr - DD_pr) /T(var);
        #if use_fixed_cov
        residuals[0] = (est_DD_pr - DD_pr) /T(pow(0.4, 2)); 
        #else
        // residuals[0] = (est_DD_pr - DD_pr) /T(var_u2i);
        
        residuals[0] = (est_DD_pr - DD_pr) /T(var);
        // std::cout<< "residuals[0]-> " << residuals[0] << std::endl;
        #endif

        return true;
    }

    DDMeasurement dd_measurement;
    Eigen::Vector3d base_pos;

};



/* DDpseudorangeVSVConstraint  :: different satellite system different sat pose
    * M: master satellite
    * i: satellite to be double-difference based on master satellite
    * r: reference satellite
    * u: use end satellite (GNSS receiver)
*/
struct DDpseudorangeVSVConstraint
{
    typedef ceres::DynamicAutoDiffCostFunction<DDpseudorangeVSVConstraint, 10>
      DDpseudorangeVSVDynaCostFunction;
    DDpseudorangeVSVConstraint(DDMeasurement dd_measurement, Eigen::Vector3d base_pos, int keyIndex)
                :dd_measurement(dd_measurement),base_pos(base_pos),keyIndex(keyIndex){}

    template <typename T>
    bool operator()(T const* const* state, T* residuals) const
    {
        /**/
        Eigen::Vector3d pose_r = base_pos;
        GNSS_Tools m_GNSS_Tools1; // utilities
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;


        /* satellite position*/
        Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
        double var_u2m = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.u_master_SV);

        Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
        double var_u2i = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.u_iSV);
        // LOG(INFO)<<"var_u2i-> " <<var_u2i;

        Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
        double var_r2m = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.r_master_SV);
        // LOG(INFO)<<"var_r2m-> " <<var_r2m;

        Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
        double var_r2i = m_GNSS_Tools1.getVarofpr_ele_SNR(dd_measurement.r_iSV);
        // LOG(INFO)<<"var_r2i-> " <<var_r2i;
        
        T est_p_r2m = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_m));
        // est_p_r2m = est_p_r2m + OMGE_ * (rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT_;
        est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;
        
        T est_p_r2i = T(m_GNSS_Tools1.getDistanceFrom2Points(pose_r, r_pose_i)); 
        est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

        T est_p_u2m;
        T delta_x1 = pow((state[keyIndex][0] - u_pose_m(0)),2);
        T delta_y1 = pow((state[keyIndex][1] - u_pose_m(1)),2);
        T delta_z1 = pow((state[keyIndex][2] - u_pose_m(2)),2);
        est_p_u2m = sqrt(delta_x1+ delta_y1 + delta_z1);
        est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*state[keyIndex][1]-u_pose_m(1)*state[keyIndex][0])/CLIGHT_;

        T est_p_u2i;
        T delta_x2 = pow((state[keyIndex][0] - u_pose_i(0)),2);
        T delta_y2 = pow((state[keyIndex][1] - u_pose_i(1)),2);
        T delta_z2 = pow((state[keyIndex][2] - u_pose_i(2)),2);
        est_p_u2i = sqrt(delta_x2+ delta_y2 + delta_z2);

        /* relatistic effects*/
        est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*state[keyIndex][1]-u_pose_i(1)*state[keyIndex][0])/CLIGHT_;
        
        /* expected DD measurments */
        T est_DD_pr = (est_p_u2i - est_p_r2i) - (est_p_u2m - est_p_r2m);

        T p_r2m = T(dd_measurement.r_master_SV.raw_pseudorange);
        T p_r2i = T(dd_measurement.r_iSV.raw_pseudorange);
        T p_u2m = T(dd_measurement.u_master_SV.raw_pseudorange);
        T p_u2i = T(dd_measurement.u_iSV.raw_pseudorange);
        
        /*observation of DD measurement*/
        T DD_pr = (p_u2i - p_r2i) - (p_u2m - p_r2m);

        // double var = 1000;
        double var = var_u2m + var_u2i + var_r2m + var_r2i;
        var = var / 4.0;
        // residuals[0] = (est_DD_pr - DD_pr) /T(var);
        #if use_fixed_cov
        residuals[0] = (est_DD_pr - DD_pr) /T(pow(0.4, 2)); 
        #else
        // residuals[0] = (est_DD_pr - DD_pr) /T(var_u2i);
        residuals[0] = (est_DD_pr - DD_pr) /T(var);
        #endif

        return true;
    }

    // Factory method to create a CostFunction from a DDpseudorangeVSVConstraint to
    // conveniently add to a ceres problem.
    static DDpseudorangeVSVDynaCostFunction* Create(DDMeasurement dd_measurement, Eigen::Vector3d base_pos, int keyIndex, std::vector<double*>* state_array, std::vector<double*>* pose_parameter_blocks, std::vector<int*>* ar_state_num) {
        
        DDpseudorangeVSVConstraint* constraint = new DDpseudorangeVSVConstraint(
            dd_measurement, base_pos, keyIndex);
        
        DDpseudorangeVSVDynaCostFunction* cost_function = new DDpseudorangeVSVDynaCostFunction(constraint);
        
        pose_parameter_blocks->clear();
        // double a[5] = {1,2,3,4,5};
        // parameter_blocks->push_back(a);
        // parameter_blocks->push_back(&((*state_array)[keyIndex]));
        // parameter_blocks->push_back(state_array[keyIndex]);
        
        for(int i = 0; i <(keyIndex+1); i++)
        {
            pose_parameter_blocks->push_back((*state_array)[i]);
            cost_function->AddParameterBlock(3 + (*ar_state_num)[i][0]);
        }
        // std::cout << "parameter_blocks.size()-> " << parameter_blocks->size() << std::endl;
        // cost_function->AddParameterBlock(1);
        // cost_function->AddParameterBlock(5);
        
        cost_function->SetNumResiduals(1);
        return (cost_function);
    }

    DDMeasurement dd_measurement;
    Eigen::Vector3d base_pos;
    int keyIndex;

};

