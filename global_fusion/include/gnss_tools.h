/*******************************************************
 * Copyright (C) 2019, Intelligent Positioning and Navigation Lab, Hong Kong Polytechnic University
 * 
 * This file is part of GraphGNSSLib.
 * Licensed under the GNU General Public License v3.0;
 * you may not use this file except in compliance with the License.
 *
 * Author: Weisong Wen (weisong.wen@connect.polyu.hk)
 *******************************************************/

#ifndef GNSS_Tools_HPP
#define GNSS_Tools_HPP
#include <nlosExclusion/GNSS_Raw_Array.h>
// google implements commandline flags processing.
#include <gflags/gflags.h>
// google loging tools
#include <glog/logging.h>

using namespace Eigen;
#define pi_ 3.1415926
#define minGPSCnt 4
#define minBeidouCnt 1
#define D2R 3.1415926/180.0
#define R2D 180.0/3.1415926

#define useEleVar 1

#define use_fixed_cov_ar 1


struct DDMeasurement
{
    /* the master satellite is found from user end*/
    nlosExclusion::GNSS_Raw u_master_SV;
    nlosExclusion::GNSS_Raw u_iSV;

    nlosExclusion::GNSS_Raw r_master_SV;
    nlosExclusion::GNSS_Raw r_iSV;

    double var_pr;
    double var_cp;
};


/**
 * @brief GNSS Tools
 * @note  GNSS related functions
 */
class GNSS_Tools {
public:
  GNSS_Tools() {
    double reserve = 0.01;
  }

public:

double getMSE(Eigen::MatrixXd est, Eigen::MatrixXd refSat, std::string format)
{
  double MSE = 0;
  if(format=="2D_Error")
  {
    MSE = sqrt(pow((est(0) - refSat(0)),2) + pow((est(1) - refSat(1)),2));
    return MSE;
  }
  if(format=="3D_Error")
  {
    MSE = sqrt(pow((est(0) - refSat(0)),2) + pow((est(1) - refSat(1)),2) + pow((est(2) - refSat(2)),2));
    return MSE;
  }
}

Eigen::Vector3d getXYZ_error(Eigen::MatrixXd est, Eigen::MatrixXd refSat, std::string format)
{
  Eigen::Vector3d xyz_error;
  if(format=="xyz_Error")
  {
    xyz_error(0) = est(0) - refSat(0);
    xyz_error(1) = est(1) - refSat(1);
    xyz_error(2) = est(2) - refSat(2);
    return xyz_error;
  }
  else 
    return xyz_error;
}

 Eigen::MatrixXd getAllPositions(nlosExclusion::GNSS_Raw_Array GNSS_data)
{
  Eigen::MatrixXd eAllSVPositions; // satellite positions   
  eAllSVPositions.resize(GNSS_data.GNSS_Raws.size(), 4);
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++) // for weighted least square
  {
    eAllSVPositions(i,0) = GNSS_data.GNSS_Raws[i].prn_satellites_index;
    eAllSVPositions(i,1) = GNSS_data.GNSS_Raws[i].sat_pos_x;
    eAllSVPositions(i,2) = GNSS_data.GNSS_Raws[i].sat_pos_y;
    eAllSVPositions(i,3) = GNSS_data.GNSS_Raws[i].sat_pos_z;
  }
  return eAllSVPositions;
}

Eigen::MatrixXd getAllMeasurements(nlosExclusion::GNSS_Raw_Array GNSS_data)
{
  Eigen::MatrixXd eAllMeasurement; // pseudorange measurements 
  eAllMeasurement.resize(GNSS_data.GNSS_Raws.size(), 3);
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++) // for weighted least square
  {
    eAllMeasurement(i,0) = GNSS_data.GNSS_Raws[i].prn_satellites_index;
    eAllMeasurement(i,1) = GNSS_data.GNSS_Raws[i].snr;
    eAllMeasurement(i,2) = GNSS_data.GNSS_Raws[i].pseudorange;
  }
  return eAllMeasurement;
}
  /*
author: WEN Weisong, visiting Ph.D student in Univeristy of California, Berkeley. (weisong.wen@berkeley.edu)
function: check GNSS availability
input: GNSS data
output: bool 
*/
bool checkAvailability(nlosExclusion::GNSS_Raw_Array GNSS_data) //
{
  bool result = true;
  if(validateSV(getGPSCnt(GNSS_data), getBeiDouCnt(GNSS_data)) && (checkRepeating(GNSS_data)))
  {
    return true;
  }
  else
  {
    return false;
  }
}

int getBeiDouCnt(nlosExclusion::GNSS_Raw_Array GNSS_data)
{
  int cnt = 0; 
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++)
  {
    if(PRNisBeidou(GNSS_data.GNSS_Raws[i].prn_satellites_index))
    {
      cnt++;
    }
  }
  return cnt;
}

int getGPSCnt(nlosExclusion::GNSS_Raw_Array GNSS_data)
{
  int cnt = 0; 
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++)
  {
    if(PRNisGPS(GNSS_data.GNSS_Raws[i].prn_satellites_index))
    {
      cnt++;
    }
  }
  return cnt;
}

bool validateSV(int gpsCnt, int BeidouCnt)
{
  if((gpsCnt >= minGPSCnt)&&(BeidouCnt >= minBeidouCnt))
  {
    return true;
  }
  else
  {
    if((gpsCnt + BeidouCnt) > 5)
    {
      return true;
    }
    LOG(ERROR) << "satellite number is not enough" <<"GPS cnt= "
              <<gpsCnt<< "  beidou cnt "<< BeidouCnt;
    return false;
  }
}

bool checkRepeating(nlosExclusion::GNSS_Raw_Array GNSS_data) 
{
  std::vector<int> satList;
  nlosExclusion::GNSS_Raw_Array GNSS_dataTmp;
  satList.clear();
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++) // 
  {
    if(std::find(satList.begin(), satList.end(), GNSS_data.GNSS_Raws[i].prn_satellites_index) != satList.end())
    {
      LOG(ERROR) <<"find repeted satellite!!!!!!!!!!!!! " << GNSS_data.GNSS_Raws[i].prn_satellites_index <<" at epoch "<<GNSS_data.GNSS_Raws[i].GNSS_time;
      return false; // find repeted satellite ID
    }
    satList.push_back(GNSS_data.GNSS_Raws[i].prn_satellites_index);
    
  }
  return true;
}

void removeRepeatedSV(nlosExclusion::GNSS_Raw_Array &GNSS_data) 
{
  std::vector<int> satList;
  nlosExclusion::GNSS_Raw_Array GNSS_dataTmp;
  satList.clear();
  for(int i =0; i < GNSS_data.GNSS_Raws.size(); i++) // 
  {
    if(std::find(satList.begin(), satList.end(), GNSS_data.GNSS_Raws[i].prn_satellites_index) != satList.end())
    {
      // LOG(ERROR) <<"find repeted satellite " << GNSS_data.GNSS_Raws[i].prn_satellites_index <<" at epoch "<<GNSS_data.GNSS_Raws[i].GNSS_time / 1e-3;
    }
    else
    {
      GNSS_dataTmp.GNSS_Raws.push_back(GNSS_data.GNSS_Raws[i]);
    }
    satList.push_back(GNSS_data.GNSS_Raws[i].prn_satellites_index);
    
  }
  GNSS_data = GNSS_dataTmp;
}

  /*
author: WEN Weisong, visiting Ph.D student in Univeristy of California, Berkeley. (weisong.wen@berkeley.edu)
function: llh to ecef
input: llh (Matrix3d)
output: ecef (Matrix3d)
*/
Eigen::MatrixXd llh2ecef(Eigen::MatrixXd data) // transform the llh to ecef
{
  Eigen::MatrixXd ecef; // the ecef for output
  ecef.resize(3, 1);
  double a = 6378137.0;
  double b = 6356752.314;
  double n, Rx, Ry, Rz;
  double lon = (double)data(0) * 3.1415926 / 180.0; // lon to radis
  double lat = (double)data(1) * 3.1415926 / 180.0; // lat to radis
  double alt = (double)data(2); // altitude
  n = a * a / sqrt(a * a * cos(lat) * cos(lat) + b * b * sin(lat) * sin(lat));
  Rx = (n + alt) * cos(lat) * cos(lon);
  Ry = (n + alt) * cos(lat) * sin(lon);
  Rz = (b * b / (a * a) * n + alt) * sin(lat);
  ecef(0) = Rx; // return value in ecef
  ecef(1) = Ry; // return value in ecef
  ecef(2) = Rz; // return value in ecef
  return ecef;

  /**************for test purpose*************************
  Eigen::MatrixXd llh;
  llh.resize(3, 1);
  Eigen::MatrixXd ecef;
  ecef.resize(3, 1);
  llh(0) = 114.1772621294604;
  llh(1) = 22.29842880200087;
  llh(2) = 58;
  ecef = llh2ecef(llh);
  cout << "ecef ->: " << ecef << "\n";
  */
}

/*
author: WEN Weisong, visiting Ph.D student in Univeristy of California, Berkeley. (weisong.wen@berkeley.edu)
function: ecef to llh
input: ecef (Matrix3d)
output: llh (Matrix3d)
*/
Eigen::MatrixXd ecef2llh(Eigen::MatrixXd data) // transform the ecef to llh
{
  Eigen::MatrixXd llh; // the ecef for output
  double pi = 3.1415926; // pi
  llh.resize(3, 1);
  double x = data(0); // obtain ecef 
  double y = data(1);
  double z = data(2);
  double x2 = pow(x, 2);
  double y2 = pow(y, 2);
  double z2 = pow(z, 2);

  double a = 6378137.0000; //earth radius in meters
  double b = 6356752.3142; // earth semiminor in meters
  double e = sqrt(1 - (b / a) * (b / a));
  double b2 = b*b;
  double e2 = e*e;
  double  ep = e*(a / b);
  double  r = sqrt(x2 + y2);
  double  r2 = r*r;
  double  E2 = a * a - b*b;
  double F = 54 * b2*z2;
  double G = r2 + (1 - e2)*z2 - e2*E2;
  double c = (e2*e2*F*r2) / (G*G*G);
  double s = (1 + c + sqrt(c*c + 2 * c));
  s = pow(s, 1 / 3);
  double P = F / (3 * ((s + 1 / s + 1)*(s + 1 / s + 1)) * G*G);
  double Q = sqrt(1 + 2 * e2*e2*P);
  double ro = -(P*e2*r) / (1 + Q) + sqrt((a*a / 2)*(1 + 1 / Q) - (P*(1 - e2)*z2) / (Q*(1 + Q)) - P*r2 / 2);
  double tmp = (r - e2*ro)*(r - e2*ro);
  double U = sqrt(tmp + z2);
  double V = sqrt(tmp + (1 - e2)*z2);
  double zo = (b2*z) / (a*V);

  double height = U*(1 - b2 / (a*V));

  double lat = atan((z + ep*ep*zo) / r);

  double temp = atan(y / x);
  double long_;
  if (x >= 0)
    long_ = temp;
  else if ((x < 0) && (y >= 0))
    long_ = pi + temp;
  else
    long_ = temp - pi;
  llh(0) = (long_)*(180 / pi);
  llh(1) = (lat)*(180 / pi);
  llh(2) = height;
  return llh;

  /**************for test purpose*************************
  Eigen::MatrixXd ecef;
  ecef.resize(3, 1);
  Eigen::MatrixXd llh;
  llh.resize(3, 1);
  ecef(0) = -2418080.9387265667;
  ecef(1) = 5386190.3905763263;
  ecef(2) = 2405041.9305451373;
  llh = ecef2llh(ecef);
  cout << "llh ->: " << llh << "\n";
  */
}

/*
author: WEN Weisong, visiting Ph.D student in Univeristy of California, Berkeley. (weisong.wen@berkeley.edu)
function: ecef to enu
input: original llh, and current ecef (Matrix3d)
output: enu (Matrix3d)
*/
Eigen::MatrixXd ecef2enu(Eigen::MatrixXd originllh, Eigen::MatrixXd ecef) // transform the ecef to enu 
{
  double pi = 3.1415926; // pi 
  double DEG2RAD = pi / 180.0;
  double RAD2DEG = 180.0 / pi;

  Eigen::MatrixXd enu; // the enu for output
  enu.resize(3, 1); // resize to 3X1
  Eigen::MatrixXd oxyz; // the original position 
  oxyz.resize(3, 1); // resize to 3X1

  double x, y, z; // save the x y z in ecef
  x = ecef(0);
  y = ecef(1);
  z = ecef(2);

  double ox, oy, oz; // save original reference position in ecef
  oxyz = llh2ecef(originllh);
  ox = oxyz(0); // obtain x in ecef 
  oy = oxyz(1); // obtain y in ecef
  oz = oxyz(2); // obtain z in ecef

  double dx, dy, dz;
  dx = x - ox;
  dy = y - oy;
  dz = z - oz;

  double lonDeg, latDeg, _; // save the origin lon alt in llh
  lonDeg = originllh(0);
  latDeg = originllh(1);
  double lon = lonDeg * DEG2RAD;
  double lat = latDeg * DEG2RAD;

  //save ENU
  enu(0) = -sin(lon) * dx + cos(lon) * dy;
  enu(1) = -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz;
  enu(2) = cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz;
  return enu;

  /**************for test purpose*****suqare distance is about 37.4 meters********************
  Eigen::MatrixXd llh;  //original
  llh.resize(3, 1);
  llh(0) = 114.1775072541416;
  llh(1) = 22.29817969722738;
  llh(2) = 58;
  Eigen::MatrixXd ecef;
  ecef.resize(3, 1);
  ecef(0) = -2418080.9387265667;
  ecef(1) = 5386190.3905763263;
  ecef(2) = 2405041.9305451373;
  Eigen::MatrixXd enu;
  enu.resize(3, 1);
  enu = ecef2enu(llh, ecef);
  cout << "enu ->: " << enu << "\n";
  */
}

/*
author: WEN Weisong, visiting Ph.D student in Univeristy of California, Berkeley. (weisong.wen@berkeley.edu)
function: ecef to enu
input: original llh, and current ecef (Matrix3d)
output: enu (Matrix3d)
*/
Eigen::MatrixXd enu2ecef(Eigen::MatrixXd originllh, Eigen::MatrixXd enu) // transform the ecef to enu 
{
  // enu to ecef
  double  e = enu(0);
  double  n = enu(1);
  double  u = enu(2);
  double lon = (double)originllh(0) * D2R;
  double lat = (double)originllh(1) * D2R;
  Eigen::MatrixXd oxyz; // the original position 
  oxyz.resize(3, 1); // resize to 3X1
  oxyz = llh2ecef(originllh);
  double ox = oxyz(0);
  double oy = oxyz(1);
  double oz = oxyz(2);

  oxyz(0) = ox - sin(lon) * e - cos(lon) * sin(lat) * n + cos(lon) * cos(lat) * u;
  oxyz(1) = oy + cos(lon) * e - sin(lon) * sin(lat) * n + cos(lat) * sin(lon) * u;
  oxyz(2) = oz + cos(lat) * n + sin(lat) * u;
  return oxyz;
}

/**
   * @brief  least square for signle point positioning
   * @param eAllSVPositions ((n,4) prn, sx, sy, sz, )     eAllSVPositions ((n,3) PRN CNO Pseudorange)
   * @return eWLSSolution 5 unknowns with two clock bias variables
   @ 
  */
  Eigen::MatrixXd LeastSquare(Eigen::MatrixXd eAllSVPositions, Eigen::MatrixXd eAllMeasurement){
  
    Eigen::MatrixXd eWLSSolution;
    eWLSSolution.resize(5, 1);

    /**after read the obs file, one measure is not right**/
    int validNumMeasure=0;
    std::vector<int> validMeasure;
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          validNumMeasure++;
          validMeasure.push_back(int(eAllMeasurement(idx, 0)));
        }
      }
    }

    Eigen::MatrixXd validMeasurement; // for WLS 
    validMeasurement.resize(validNumMeasure,eAllMeasurement.cols());
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllMeasurement.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            validMeasurement(idx, kdx) = eAllMeasurement(idx, kdx);
            
          }
        }
      }
    }



    int iNumSV = validMeasurement.rows();

    /*Find the received SV and Sort based on the order of Measurement matrix*/
    Eigen::MatrixXd eExistingSVPositions; // for WLS
    eExistingSVPositions.resize(iNumSV, eAllSVPositions.cols());

    for (int idx = 0; idx < validMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(validMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllSVPositions.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            eExistingSVPositions(idx, kdx) = eAllSVPositions(jdx, kdx);
            
          }
        }
      }
    } 
    //for (int idx = 0; idx < eExistingSVPositions.rows(); idx++){
    //  printf("%2d-[%3d] - (%10.2f,%10.2f,%10.2f) %f\n", idx, int(eExistingSVPositions(idx, 0)), eExistingSVPositions(idx, 1), eExistingSVPositions(idx, 2), eExistingSVPositions(idx, 3), eExistingSVPositions(idx, 4)*CLIGHT);
    //}

    //Intialize the result by guessing.
    for (int idx = 0; idx < eWLSSolution.rows(); idx++){
      eWLSSolution(idx, 0) = 0;
    }
    
    // for the case of insufficient satellite
    if (iNumSV < 5){
      return eWLSSolution;
    }

    bool bWLSConverge = false;

    int count = 0;
    while (!bWLSConverge)
    {
      Eigen::MatrixXd eH_Matrix;
      eH_Matrix.resize(iNumSV, eWLSSolution.rows());

      Eigen::MatrixXd eDeltaPr;
      eDeltaPr.resize(iNumSV, 1);

      Eigen::MatrixXd eDeltaPos;
      eDeltaPos.resize(eWLSSolution.rows(), 1);

      for (int idx = 0; idx < iNumSV; idx++){

        int prn = int(validMeasurement(idx, 0));
        double pr = validMeasurement(idx, 2);
        
        // Calculating Geometric Distance
        double rs[3], rr[3], e[3];
        double dGeoDistance;

        rs[0] = eExistingSVPositions(idx, 1);
        rs[1] = eExistingSVPositions(idx, 2);
        rs[2] = eExistingSVPositions(idx, 3);

        rr[0] = eWLSSolution(0);
        rr[1] = eWLSSolution(1);
        rr[2] = eWLSSolution(2);

        // dGeoDistance = geodist(rs, rr, e);
        dGeoDistance = sqrt(pow((rs[0] - rr[0]),2) + pow((rs[1] - rr[1]),2) +pow((rs[2] - rr[2]),2));

        // Making H matrix      
        eH_Matrix(idx, 0) = -(rs[0] - rr[0]) / dGeoDistance;
        eH_Matrix(idx, 1) = -(rs[1] - rr[1]) / dGeoDistance;
        eH_Matrix(idx, 2) = -(rs[2] - rr[2]) / dGeoDistance;

        if (PRNisGPS(prn)){
          eH_Matrix(idx, 3) = 1;
          eH_Matrix(idx, 4) = 0;
        }
        else if (PRNisBeidou(prn))
        {
          eH_Matrix(idx, 3) = 1;
          eH_Matrix(idx, 4) = 1;
        }

        // Making delta pseudorange
        double rcv_clk_bias;
        if (PRNisGPS(prn)){
          rcv_clk_bias = eWLSSolution(3);       
        }
        else if (PRNisBeidou(prn))
        {
          rcv_clk_bias = eWLSSolution(4);
        }
        // double sv_clk_bias = eExistingSVPositions(idx, 4) * CLIGHT;
        eDeltaPr(idx, 0) = pr - dGeoDistance + rcv_clk_bias;
        //printf("%2d - %f %f %f %f \n", prn, pr, dGeoDistance, eDeltaPr(idx, 0), rcv_clk_bias);
      }

      // Least Square Estimation 
      eDeltaPos = (eH_Matrix.transpose() * eH_Matrix).ldlt().solve(eH_Matrix.transpose() *  eDeltaPr);
      //eDeltaPos = (eH_Matrix.transpose() * eH_Matrix).inverse() * eH_Matrix.transpose() *  eDeltaPr;
      //eDeltaPos = eH_Matrix.householderQr().solve(eDeltaPr);

      //for (int idx = 0; idx < eDeltaPos.rows(); idx++)
      //  printf("%f ", eDeltaPos(idx));
      //printf("\n");

      eWLSSolution(0) += eDeltaPos(0);
      eWLSSolution(1) += eDeltaPos(1);
      eWLSSolution(2) += eDeltaPos(2);
      eWLSSolution(3) += eDeltaPos(3);
      eWLSSolution(4) += eDeltaPos(4);

      for (int i = 0; i < 3; ++i){
        //printf("%f\n", fabs(eDeltaPos(i)));
        if (fabs(eDeltaPos(i)) >1e-4)
        {
          bWLSConverge = false;
        }
        else { 
          bWLSConverge = true;
        };
        
      }
      count += 1;
      if (count > 6)
        bWLSConverge = true;
    }
    // printf("WLS -> (%11.2f,%11.2f,%11.2f)\n\n", eWLSSolution(0), eWLSSolution(1), eWLSSolution(2));
    std::cout << std::setprecision(12);
    // cout<< "---------------WLS (ECEF) x, y, z, bias_gps, bias_beidou-----------------  \n"<<eWLSSolution<<endl;

    return eWLSSolution;
  }

  /**
   * @brief  weighted least square for signle point positioning
   * @param eAllSVPositions ((n,4) prn, sx, sy, sz, )     eAllSVPositions ((n,3) PRN CNO Pseudorange)
   * @return eWLSSolution 5 unknowns with two clock bias variables
   @ 
  */
  Eigen::MatrixXd WeightedLeastSquare(Eigen::MatrixXd eAllSVPositions, Eigen::MatrixXd eAllMeasurement, nlosExclusion::GNSS_Raw_Array GNSS_data, std::string method){
  
    Eigen::MatrixXd eWLSSolution;
    eWLSSolution.resize(5, 1);

    MatrixXd weight_matrix; //goGPS weighting
    if(method == "WLS")
    {
      weight_matrix = cofactorMatrixCal_WLS(GNSS_data, method); //goGPS weighting;
    }
    else if(method == "R-WLS")
    {
      weight_matrix = cofactorMatrixCal_R_WLS(GNSS_data, method); //goGPS weighting;
    }
    else if (method == "LS")
    {
      weight_matrix = cofactorMatrixCal_WLS(GNSS_data, method); //goGPS weighting;
      weight_matrix.setIdentity();
    }

    // std::cout << " weight_matrix-> " << weight_matrix<< std::endl;    

    /**after read the obs file, one measure is not right**/
    int validNumMeasure=0;
    std::vector<int> validMeasure;
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          validNumMeasure++;
          validMeasure.push_back(int(eAllMeasurement(idx, 0)));
        }
      }
    }

    Eigen::MatrixXd validMeasurement; // for WLS 
    validMeasurement.resize(validNumMeasure,eAllMeasurement.cols());
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllMeasurement.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            validMeasurement(idx, kdx) = eAllMeasurement(idx, kdx);
            
          }
        }
      }
    }



    int iNumSV = validMeasurement.rows();

    /*Find the received SV and Sort based on the order of Measurement matrix*/
    Eigen::MatrixXd eExistingSVPositions; // for WLS
    eExistingSVPositions.resize(iNumSV, eAllSVPositions.cols());

    for (int idx = 0; idx < validMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(validMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllSVPositions.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            eExistingSVPositions(idx, kdx) = eAllSVPositions(jdx, kdx);
            
          }
        }
      }
    } 
    //for (int idx = 0; idx < eExistingSVPositions.rows(); idx++){
    //  printf("%2d-[%3d] - (%10.2f,%10.2f,%10.2f) %f\n", idx, int(eExistingSVPositions(idx, 0)), eExistingSVPositions(idx, 1), eExistingSVPositions(idx, 2), eExistingSVPositions(idx, 3), eExistingSVPositions(idx, 4)*CLIGHT);
    //}

    //Intialize the result by guessing.
    for (int idx = 0; idx < eWLSSolution.rows(); idx++){
      eWLSSolution(idx, 0) = 0;
    }
    
    // for the case of insufficient satellite
    if (iNumSV < 5){
      return eWLSSolution;
    }

    bool bWLSConverge = false;

    int count = 0;
    Eigen::MatrixXd eH_Matrix_;
    while (!bWLSConverge)
    {
      Eigen::MatrixXd eH_Matrix;
      eH_Matrix.resize(iNumSV, eWLSSolution.rows());

      Eigen::MatrixXd eDeltaPr;
      eDeltaPr.resize(iNumSV, 1);

      Eigen::MatrixXd eDeltaPos;
      eDeltaPos.resize(eWLSSolution.rows(), 1);

      for (int idx = 0; idx < iNumSV; idx++){

        int prn = int(validMeasurement(idx, 0));
        double pr = validMeasurement(idx, 2);
        
        // Calculating Geometric Distance
        double rs[3], rr[3], e[3];
        double dGeoDistance;

        rs[0] = eExistingSVPositions(idx, 1);
        rs[1] = eExistingSVPositions(idx, 2);
        rs[2] = eExistingSVPositions(idx, 3);

        rr[0] = eWLSSolution(0);
        rr[1] = eWLSSolution(1);
        rr[2] = eWLSSolution(2);

        // dGeoDistance = geodist(rs, rr, e);
        dGeoDistance = sqrt(pow((rs[0] - rr[0]),2) + pow((rs[1] - rr[1]),2) +pow((rs[2] - rr[2]),2));
        
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;
        dGeoDistance = dGeoDistance + OMGE_ * (rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT_;

        // Making H matrix      
        eH_Matrix(idx, 0) = -(rs[0] - rr[0]) / dGeoDistance;
        eH_Matrix(idx, 1) = -(rs[1] - rr[1]) / dGeoDistance;
        eH_Matrix(idx, 2) = -(rs[2] - rr[2]) / dGeoDistance;

        if (PRNisGPS(prn)){
          eH_Matrix(idx, 3) = 1;
          eH_Matrix(idx, 4) = 0;
        }
        else
        {
          eH_Matrix(idx, 3) = 0;
          eH_Matrix(idx, 4) = 1;
        }

        // Making delta pseudorange
        double rcv_clk_bias = 0;
        if (PRNisGPS(prn)){
          rcv_clk_bias = eWLSSolution(3);       
        }
        // else if (PRNisBeidou(prn))
        else
        {
          rcv_clk_bias = eWLSSolution(4);
        }
        // double sv_clk_bias = eExistingSVPositions(idx, 4) * CLIGHT;
        eDeltaPr(idx, 0) = pr - dGeoDistance - rcv_clk_bias;
        // printf("%2d - %f %f %f %f \n", prn, pr, dGeoDistance, eDeltaPr(idx, 0), rcv_clk_bias);
      }

      // Least Square Estimation 
      // eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).ldlt().solve(eH_Matrix.transpose() * weight_matrix *  eDeltaPr);
      eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).inverse() * eH_Matrix.transpose() * weight_matrix * eDeltaPr;
      //eDeltaPos = eH_Matrix.householderQr().solve(eDeltaPr);

      //for (int idx = 0; idx < eDeltaPos.rows(); idx++)
      //  printf("%f ", eDeltaPos(idx));
      //printf("\n");

      eWLSSolution(0) += eDeltaPos(0);
      eWLSSolution(1) += eDeltaPos(1);
      eWLSSolution(2) += eDeltaPos(2);
      eWLSSolution(3) += eDeltaPos(3);
      eWLSSolution(4) += eDeltaPos(4);

      eH_Matrix_ = eH_Matrix;

      for (int i = 0; i < 3; ++i){
        //printf("%f\n", fabs(eDeltaPos(i)));
        if (fabs(eDeltaPos(i)) >1e-4)
        {
          bWLSConverge = false;
        }
        else { 
          bWLSConverge = true;
          // std::cout<< "eDeltaPr-> " << eDeltaPr << std::endl;
        };
        
      }
      count += 1;
      if (count > 10)
      {
        bWLSConverge = true;
        // std::cout<<" more than 10 times in iterations"<<std::endl;
        
      }
      
    }
    // printf("WLS -> (%11.2f,%11.2f,%11.2f)\n\n", eWLSSolution(0), eWLSSolution(1), eWLSSolution(2));
    std::cout << std::setprecision(12);
    
    #if 0 // debug result
    std::cout<<"eH_Matrix-> \n" << eH_Matrix_ << std::endl;
    std::cout<<"iterations-> " << count << std::endl;
    std::cout<< "---------------WLS (ECEF) x, y, z, bias_gps, bias_beidou-----------------  \n"<<eWLSSolution<<std::endl;
    #endif

    return eWLSSolution;
  }

  Eigen::MatrixXd getPseudorangeResidual(Eigen::MatrixXd eWLSSolution, nlosExclusion::GNSS_Raw_Array GNSS_data)
  {
    Eigen::MatrixXd eWLSSolutionResidual;
    eWLSSolutionResidual.resize(GNSS_data.GNSS_Raws.size(),1);
    for(int i = 0; i <GNSS_data.GNSS_Raws.size(); i++)
    {
      int prn = int(GNSS_data.GNSS_Raws[i].prn_satellites_index);
      double pr = GNSS_data.GNSS_Raws[i].pseudorange;
      
      // Calculating Geometric Distance
      double rs[3], rr[3], e[3];
      double dGeoDistance;

      rs[0] = GNSS_data.GNSS_Raws[i].sat_pos_x;
      rs[1] = GNSS_data.GNSS_Raws[i].sat_pos_y;
      rs[2] = GNSS_data.GNSS_Raws[i].sat_pos_z;

      rr[0] = eWLSSolution(0);
      rr[1] = eWLSSolution(1);
      rr[2] = eWLSSolution(2);

      dGeoDistance = sqrt(pow((rs[0] - rr[0]),2) + pow((rs[1] - rr[1]),2) +pow((rs[2] - rr[2]),2));

      // Making delta pseudorange
      double rcv_clk_bias;
      if (PRNisGPS(prn)){
        rcv_clk_bias = eWLSSolution(3);       
      }
      else if (PRNisBeidou(prn))
      {
        rcv_clk_bias = eWLSSolution(4);
      }
      // double sv_clk_bias = eExistingSVPositions(idx, 4) * CLIGHT;
      eWLSSolutionResidual(i, 0) = pr - dGeoDistance - rcv_clk_bias;
      // std::cout <<"pr-> "<< pr <<  "  dGeoDistance-> "<<dGeoDistance <<"  rcv_clk_biasG-> "<< eWLSSolution(3) << "  rcv_clk_biasB-> "<< eWLSSolution(4)<< std::endl;
      // std::cout <<"pr-> "<< pr <<  "  pr - dGeoDistance-> "<<pr - dGeoDistance<< std::endl;
    }
    return eWLSSolutionResidual;
  }

/**
   * @brief satellite set validation
   * @param prn
   * @return ture/false
   @ 
   */
  bool PRNisGPS(int prn)
  {
    if (prn <= 32 || prn == 84)
      return true;
    else{
      return false;
    } 
  }

  /**
   * @brief satellite set validation
   * @param prn
   * @return ture/false
   @ 
   */
  bool PRNisGLONASS(int prn)
  {
    if (prn > 32 && prn <= 56)
      return true;
    else{
      return false;
    }
  }

  /**
   * @brief satellite set validation
   * @param prn
   * @return ture/false
   @ 
   */
  bool PRNisBeidou(int prn)
  {
    if ((prn <= 121) && (prn >= 87))
      return true;
    else{
      return false;
    }
  }



  /**
   * @brief covariance estimation
   * @param nlosExclusion::GNSS_Raw_Array GNSS_data
   * @return weight_matrix
   @ 
   */
  Eigen::MatrixXd cofactorMatrixCal_WLS(nlosExclusion::GNSS_Raw_Array GNSS_data, std::string method)
  {
    Eigen::Matrix<double,4,1> parameters;
    parameters<<50.0, 30.0, 30.0, 10.0; // loosely coupled 
    // parameters<<50.0, 30.0, 20.0, 30.0; // loosely coupled 
    double snr_1 = parameters(0); // T = 50
    double snr_A = parameters(1); // A = 30
    double snr_a = parameters(2);// a = 30
    double snr_0 = parameters(3); // F = 10

    VectorXd cofactor_;  // cofactor of satellite
    cofactor_.resize(GNSS_data.GNSS_Raws.size());
    for(int i = 0; i < GNSS_data.GNSS_Raws.size(); i++)
    {
      if( 1 )
      {
        double snr_R = GNSS_data.GNSS_Raws[i].snr;
        double elR = GNSS_data.GNSS_Raws[i].elevation;
        double q_R_1 = 1 / (pow(( sin(elR * 3.1415926/180.0 )),2));
        double q_R_2 = pow(10,(-(snr_R - snr_1) / snr_a));
        double q_R_3 = (((snr_A / (pow(10,(-(snr_0 - snr_1) / snr_a))) - 1) / (snr_0 - snr_1)) * (snr_R - snr_1) + 1);
        double q_R = q_R_1* (q_R_2 * q_R_3);
        cofactor_[i]=(1.0/float(q_R)); // uncertainty: cofactor_[i] larger, larger uncertainty
      }
    }
    // cout<<"cofactor_ -> "<<cofactor_<<endl;

    MatrixXd weight_matrix;
    weight_matrix.resize(GNSS_data.GNSS_Raws.size(),GNSS_data.GNSS_Raws.size());
    weight_matrix.setIdentity();
    for(int k = 0; k < weight_matrix.rows(); k++)
    {
      weight_matrix.row(k) = weight_matrix.row(k) * cofactor_(k);
    }
    if(method == "WLS")
    {
      return weight_matrix;
    }
    else
    {
      weight_matrix.setIdentity();
      return weight_matrix;
    }
  }

    /**
   * @brief covariance estimation for R-WLS
   * @param nlosExclusion::GNSS_Raw_Array GNSS_data
   * @return weight_matrix
   @ 
   */
  Eigen::MatrixXd cofactorMatrixCal_R_WLS(nlosExclusion::GNSS_Raw_Array GNSS_data, std::string method)
  {
    Eigen::Matrix<double,4,1> parameters;
    double factor = 0.75;
    // parameters<<50.0, 30.0, 30.0, 10.0; // original 
    // parameters<<50.0, 30.0, 20.0, 30.0; // remodeling 
    parameters<<50.0, 30.0, 20.0, 35.0; // remodeling 
    double snr_1 = parameters(0); // T = 50
    double snr_A = parameters(1); // A = 30
    double snr_a = parameters(2);// a = 30
    double snr_0 = parameters(3); // F = 10

    // double snr_1 = 60.0; // T = 50
    // double snr_A = 40.0; // A = 30
    // double snr_a = 30.0;// a = 30
    // double snr_0 = 10.0; // F = 10
    VectorXd cofactor_;  // cofactor of satellite
    cofactor_.resize(GNSS_data.GNSS_Raws.size());
    for(int i = 0; i < GNSS_data.GNSS_Raws.size(); i++)
    {
      if( (PRNisGPS(GNSS_data.GNSS_Raws[i].prn_satellites_index)) || (PRNisBeidou(GNSS_data.GNSS_Raws[i].prn_satellites_index)) )
      {
        double snr_R = GNSS_data.GNSS_Raws[i].snr;
        double elR = GNSS_data.GNSS_Raws[i].elevation;
        double q_R_1 = 1 / (pow(( sin(elR * 3.1415926/180.0 )),2));
        double q_R_2 = pow(10,(-(snr_R - snr_1) / snr_a));
        double q_R_3 = (((snr_A / (pow(10,(-(snr_0 - snr_1) / snr_a))) - 1) / (snr_0 - snr_1)) * (snr_R - snr_1) + 1);
        double q_R = q_R_1* (q_R_2 * q_R_3);
        cofactor_[i]=(1.0/float(q_R)); // uncertainty: cofactor_[i] larger, larger uncertainty
        if(GNSS_data.GNSS_Raws[i].visable==0)
        {
          cofactor_[i] = factor * cofactor_[i]; // NLOS with larger uncertainty 
        }
      }
    }
    // std::cout<<"cofactor_ -> "<<cofactor_<<std::endl;

    MatrixXd weight_matrix;
    weight_matrix.resize(GNSS_data.GNSS_Raws.size(),GNSS_data.GNSS_Raws.size());
    weight_matrix.setIdentity();
    for(int k = 0; k < weight_matrix.rows(); k++)
    {
      weight_matrix.row(k) = weight_matrix.row(k) * cofactor_(k);
    }
    if(method == "R-WLS")
    {
      return weight_matrix;
    }
    else
    {
      weight_matrix.setIdentity();
      return weight_matrix;
    }
  }



/**
   * @brief  weighted least square for signle point positioning
   * @param eAllSVPositions ((n,4) prn, sx, sy, sz, )     eAllSVPositions ((n,3) PRN CNO Pseudorange)
   * @return eWLSSolution 5 unknowns with two clock bias variables
   @ 
  */
  Eigen::MatrixXd WeightedLeastSquare_GPS(Eigen::MatrixXd eAllSVPositions, Eigen::MatrixXd eAllMeasurement, nlosExclusion::GNSS_Raw_Array GNSS_data){
  
    Eigen::MatrixXd eWLSSolution;
    eWLSSolution.resize(4, 1);

    MatrixXd weight_matrix = cofactorMatrixCal_WLS(GNSS_data, "WLS");

    /**after read the obs file, one measure is not right**/
    int validNumMeasure=0;
    std::vector<int> validMeasure;
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          validNumMeasure++;
          validMeasure.push_back(int(eAllMeasurement(idx, 0)));
        }
      }
    }

    Eigen::MatrixXd validMeasurement; // for WLS 
    validMeasurement.resize(validNumMeasure,eAllMeasurement.cols());
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllMeasurement.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            validMeasurement(idx, kdx) = eAllMeasurement(idx, kdx);
            
          }
        }
      }
    }



    int iNumSV = validMeasurement.rows();

    /*Find the received SV and Sort based on the order of Measurement matrix*/
    Eigen::MatrixXd eExistingSVPositions; // for WLS
    eExistingSVPositions.resize(iNumSV, eAllSVPositions.cols());

    for (int idx = 0; idx < validMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(validMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllSVPositions.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            eExistingSVPositions(idx, kdx) = eAllSVPositions(jdx, kdx);
            
          }
        }
      }
    } 
    //for (int idx = 0; idx < eExistingSVPositions.rows(); idx++){
    //  printf("%2d-[%3d] - (%10.2f,%10.2f,%10.2f) %f\n", idx, int(eExistingSVPositions(idx, 0)), eExistingSVPositions(idx, 1), eExistingSVPositions(idx, 2), eExistingSVPositions(idx, 3), eExistingSVPositions(idx, 4)*CLIGHT);
    //}

    //Intialize the result by guessing.
    for (int idx = 0; idx < eWLSSolution.rows(); idx++){
      eWLSSolution(idx, 0) = 0;
    }
    
    // for the case of insufficient satellite
    if (iNumSV < 5){
      std::cout<<"satellite number is not enough" <<std::endl;
      return eWLSSolution;
    }

    bool bWLSConverge = false;

    int count = 0;
    while (!bWLSConverge)
    {
      Eigen::MatrixXd eH_Matrix;
      eH_Matrix.resize(iNumSV, eWLSSolution.rows());

      Eigen::MatrixXd eDeltaPr;
      eDeltaPr.resize(iNumSV, 1);

      Eigen::MatrixXd eDeltaPos;
      eDeltaPos.resize(eWLSSolution.rows(), 1);

      for (int idx = 0; idx < iNumSV; idx++){

        int prn = int(validMeasurement(idx, 0));
        double pr = validMeasurement(idx, 2);
        
        // Calculating Geometric Distance
        double rs[3], rr[3], e[3];
        double dGeoDistance;

        rs[0] = eExistingSVPositions(idx, 1);
        rs[1] = eExistingSVPositions(idx, 2);
        rs[2] = eExistingSVPositions(idx, 3);

        rr[0] = eWLSSolution(0);
        rr[1] = eWLSSolution(1);
        rr[2] = eWLSSolution(2);

        // dGeoDistance = geodist(rs, rr, e);
        dGeoDistance = sqrt(pow((rs[0] - rr[0]),2) + pow((rs[1] - rr[1]),2) +pow((rs[2] - rr[2]),2));
        
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;
        dGeoDistance = dGeoDistance + OMGE_ * (rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT_;

        // Making H matrix      
        eH_Matrix(idx, 0) = -(rs[0] - rr[0]) / dGeoDistance;
        eH_Matrix(idx, 1) = -(rs[1] - rr[1]) / dGeoDistance;
        eH_Matrix(idx, 2) = -(rs[2] - rr[2]) / dGeoDistance;

        if (PRNisGPS(prn)){
          eH_Matrix(idx, 3) = 1;
        }
        // Making delta pseudorange
        double rcv_clk_bias;
        if (PRNisGPS(prn)){
          rcv_clk_bias = eWLSSolution(3);       
        }
        // double sv_clk_bias = eExistingSVPositions(idx, 4) * CLIGHT;
        eDeltaPr(idx, 0) = pr - dGeoDistance + rcv_clk_bias;
        //printf("%2d - %f %f %f %f \n", prn, pr, dGeoDistance, eDeltaPr(idx, 0), rcv_clk_bias);
      }

      // Least Square Estimation 
      // eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).ldlt().solve(eH_Matrix.transpose() * weight_matrix *  eDeltaPr);
      eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).inverse() * eH_Matrix.transpose() * weight_matrix * eDeltaPr;
      //eDeltaPos = eH_Matrix.householderQr().solve(eDeltaPr);

      //for (int idx = 0; idx < eDeltaPos.rows(); idx++)
      //  printf("%f ", eDeltaPos(idx));
      //printf("\n");

      eWLSSolution(0) += eDeltaPos(0);
      eWLSSolution(1) += eDeltaPos(1);
      eWLSSolution(2) += eDeltaPos(2);
      eWLSSolution(3) += eDeltaPos(3);

      for (int i = 0; i < 3; ++i){
        //printf("%f\n", fabs(eDeltaPos(i)));
        if (fabs(eDeltaPos(i)) >1e-4)
        {
          bWLSConverge = false;
        }
        else { 
          bWLSConverge = true;
        };
        
      }
      count += 1;
      if (count > 6)
      {
        bWLSConverge = true;
        std::cout<<" more than 6 times in iterations"<<std::endl;
      }
    }
    // printf("WLS -> (%11.2f,%11.2f,%11.2f)\n\n", eWLSSolution(0), eWLSSolution(1), eWLSSolution(2));
    std::cout << std::setprecision(12);
    // std::cout<< "---------------WLS (ECEF) x, y, z, bias_gps, bias_beidou-----------------  \n"<<eWLSSolution<<std::endl;

    return eWLSSolution;
  }



  /**
   * @brief  get covariance matrix
   * @param eAllSVPositions ((n,4) prn, sx, sy, sz, )     eAllSVPositions ((n,3) PRN CNO Pseudorange)
   * @return covariance matrix for the estimated states
   @ 
  */
  Eigen::MatrixXd getCovarianceMatrix(Eigen::MatrixXd eAllSVPositions, Eigen::MatrixXd eAllMeasurement, nlosExclusion::GNSS_Raw_Array GNSS_data, std::string method){
  
    Eigen::MatrixXd eWLSSolution;
    eWLSSolution.resize(5, 1);

    MatrixXd weight_matrix; //goGPS weighting
    if(method == "WLS")
    {
      weight_matrix = cofactorMatrixCal_WLS(GNSS_data, method); //goGPS weighting;
    }
    else if(method == "R-WLS")
    {
      weight_matrix = cofactorMatrixCal_R_WLS(GNSS_data, method); //goGPS weighting;
    }
    else if (method == "LS")
    {
      weight_matrix = cofactorMatrixCal_WLS(GNSS_data, method); //goGPS weighting;
      weight_matrix.setIdentity();
    }

    // std::cout << " weight_matrix-> " << weight_matrix<< std::endl;    

    /**after read the obs file, one measure is not right**/
    int validNumMeasure=0;
    std::vector<int> validMeasure;
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          validNumMeasure++;
          validMeasure.push_back(int(eAllMeasurement(idx, 0)));
        }
      }
    }

    Eigen::MatrixXd validMeasurement; // for WLS 
    validMeasurement.resize(validNumMeasure,eAllMeasurement.cols());
    for (int idx = 0; idx < eAllMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(eAllMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllMeasurement.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            validMeasurement(idx, kdx) = eAllMeasurement(idx, kdx);
            
          }
        }
      }
    }



    int iNumSV = validMeasurement.rows();

    /*Find the received SV and Sort based on the order of Measurement matrix*/
    Eigen::MatrixXd eExistingSVPositions; // for WLS
    eExistingSVPositions.resize(iNumSV, eAllSVPositions.cols());

    for (int idx = 0; idx < validMeasurement.rows(); idx++){
      for (int jdx = 0; jdx < eAllSVPositions.rows(); jdx++){
        if (int(validMeasurement(idx, 0)) == int(eAllSVPositions(jdx, 0))){
          for (int kdx = 0; kdx < eAllSVPositions.cols(); kdx++){
            // std::cout<<"satellite prn -> "<<eAllMeasurement(idx, 0)<<"\n"<<std::endl;
            eExistingSVPositions(idx, kdx) = eAllSVPositions(jdx, kdx);
            
          }
        }
      }
    } 
    //for (int idx = 0; idx < eExistingSVPositions.rows(); idx++){
    //  printf("%2d-[%3d] - (%10.2f,%10.2f,%10.2f) %f\n", idx, int(eExistingSVPositions(idx, 0)), eExistingSVPositions(idx, 1), eExistingSVPositions(idx, 2), eExistingSVPositions(idx, 3), eExistingSVPositions(idx, 4)*CLIGHT);
    //}

    //Intialize the result by guessing.
    for (int idx = 0; idx < eWLSSolution.rows(); idx++){
      eWLSSolution(idx, 0) = 0;
    }
    
    // for the case of insufficient satellite
    // if (iNumSV < 5){
    //   return eWLSSolution;
    // }

    bool bWLSConverge = false;

    int count = 0;
    Eigen::MatrixXd covarianceMatrix;
    covarianceMatrix.resize(5,5);
    covarianceMatrix.setIdentity();
    while (!bWLSConverge && (iNumSV>5))
    {
      Eigen::MatrixXd eH_Matrix;
      eH_Matrix.resize(iNumSV, eWLSSolution.rows());

      Eigen::MatrixXd eDeltaPr;
      eDeltaPr.resize(iNumSV, 1);

      Eigen::MatrixXd eDeltaPos;
      eDeltaPos.resize(eWLSSolution.rows(), 1);

      for (int idx = 0; idx < iNumSV; idx++){

        int prn = int(validMeasurement(idx, 0));
        double pr = validMeasurement(idx, 2);
        
        // Calculating Geometric Distance
        double rs[3], rr[3], e[3];
        double dGeoDistance;

        rs[0] = eExistingSVPositions(idx, 1);
        rs[1] = eExistingSVPositions(idx, 2);
        rs[2] = eExistingSVPositions(idx, 3);

        rr[0] = eWLSSolution(0);
        rr[1] = eWLSSolution(1);
        rr[2] = eWLSSolution(2);

        // dGeoDistance = geodist(rs, rr, e);
        dGeoDistance = sqrt(pow((rs[0] - rr[0]),2) + pow((rs[1] - rr[1]),2) +pow((rs[2] - rr[2]),2));
        
        double OMGE_ = 7.2921151467E-5;
        double CLIGHT_ = 299792458.0;
        dGeoDistance = dGeoDistance + OMGE_ * (rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT_;

        // Making H matrix      
        eH_Matrix(idx, 0) = -(rs[0] - rr[0]) / dGeoDistance;
        eH_Matrix(idx, 1) = -(rs[1] - rr[1]) / dGeoDistance;
        eH_Matrix(idx, 2) = -(rs[2] - rr[2]) / dGeoDistance;

        if (PRNisGPS(prn)){
          eH_Matrix(idx, 3) = 1;
          eH_Matrix(idx, 4) = 0;
        }
        else
        {
          eH_Matrix(idx, 3) = 0;
          eH_Matrix(idx, 4) = 1;
        }

        // Making delta pseudorange
        double rcv_clk_bias = 0;
        if (PRNisGPS(prn)){
          rcv_clk_bias = eWLSSolution(3);       
        }
        // else if (PRNisBeidou(prn))
        else
        {
          rcv_clk_bias = eWLSSolution(4);
        }
        // double sv_clk_bias = eExistingSVPositions(idx, 4) * CLIGHT;
        eDeltaPr(idx, 0) = pr - dGeoDistance - rcv_clk_bias;
        // printf("%2d - %f %f %f %f \n", prn, pr, dGeoDistance, eDeltaPr(idx, 0), rcv_clk_bias);
      }

      // Least Square Estimation 
      // eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).ldlt().solve(eH_Matrix.transpose() * weight_matrix *  eDeltaPr);
      eDeltaPos = (eH_Matrix.transpose() * weight_matrix * eH_Matrix).inverse() * eH_Matrix.transpose() * weight_matrix * eDeltaPr;
      //eDeltaPos = eH_Matrix.householderQr().solve(eDeltaPr);

      //for (int idx = 0; idx < eDeltaPos.rows(); idx++)
      //  printf("%f ", eDeltaPos(idx));
      //printf("\n");

      eWLSSolution(0) += eDeltaPos(0);
      eWLSSolution(1) += eDeltaPos(1);
      eWLSSolution(2) += eDeltaPos(2);
      eWLSSolution(3) += eDeltaPos(3);
      eWLSSolution(4) += eDeltaPos(4);

      // eH_Matrix_.resize(iNumSV, 4);
      // for(int i =0; i < iNumSV; i++)
      //   for(int j =0; j < 4; j++)
      //   {
      //     eH_Matrix_(i,j) = 1.0/ eH_Matrix(i,j);
      //     if(j==3) eH_Matrix_(i,j) = 1;
      //   }
      // eH_Matrix_ = eH_Matrix.block<a, 4>(0,0);
      covarianceMatrix = (eH_Matrix.transpose() * weight_matrix * eH_Matrix);
      // covarianceMatrix = (eH_Matrix.transpose()  * eH_Matrix);
      covarianceMatrix = covarianceMatrix.inverse();
      // for(int i = 0; i < covarianceMatrix.rows(); i++)
      //   for(int j = 0; j < covarianceMatrix.cols(); j++)
      //   {
      //     covarianceMatrix(i,j) = 1.0 / covarianceMatrix(i,j);
      //   }

      for (int i = 0; i < 3; ++i){
        //printf("%f\n", fabs(eDeltaPos(i)));
        if (fabs(eDeltaPos(i)) >1e-4)
        {
          bWLSConverge = false;
        }
        else { 
          bWLSConverge = true;
        };
        
      }
      count += 1;
      if (count > 10)
      {
        bWLSConverge = true;
        // std::cout<<" more than 10 times in iterations"<<std::endl;
        
      }
      
    }
    // printf("WLS -> (%11.2f,%11.2f,%11.2f)\n\n", eWLSSolution(0), eWLSSolution(1), eWLSSolution(2));
    std::cout << std::setprecision(12);
    
    #if 0 // debug result
    std::cout<<"eH_Matrix-> \n" << eH_Matrix_ << std::endl;
    std::cout<<"iterations-> " << count << std::endl;
    std::cout<< "---------------WLS (ECEF) x, y, z, bias_gps, bias_beidou-----------------  \n"<<eWLSSolution<<std::endl;
    #endif

    return covarianceMatrix;
  }

  double getDistanceFrom2Points(Eigen::Vector3d p1, Eigen::Vector3d p2)
  {
      double xMod = pow((p1.x() - p2.x()), 2);
      double yMod = pow((p1.y() - p2.y()), 2);
      double zMod = pow((p1.z() - p2.z()), 2);
      double mod = sqrt(xMod + yMod + zMod);
      return mod;
  }

  /* getFullRankARCovMatrix
  * defi_cov_ar: the input defi_cov_ar is the deficient matrix from ceres-solver
  * n: number of phase-bias 
  * r: reference satellite
  * u: use end satellite (GNSS receiver)
  */
  void getFullRankARCovMatrix(Eigen::MatrixXd defi_cov_ar, Eigen::MatrixXd& fi_cov_ar, int& n)
  {
    bool find_zero = false;
    for(int i = 0; i < defi_cov_ar.rows(); i++)
    {
      if((defi_cov_ar(i,0)==0) && (find_zero==0))
      {
        find_zero = 1;
        n = i ;
      }
    }

    fi_cov_ar.resize(n,n);
    for(int i = 0; i<n; i++)
      for(int j = 0; j<n; j++)
      {
        fi_cov_ar(i,j) = defi_cov_ar(i,j);
      }
  }

  /* get pseudorange double-differenced Jacobian matrix */
  void getPrDDJacobian(Eigen::Vector3d u_pose, Eigen::Vector3d base_pose, DDMeasurement dd_measurement, Eigen::MatrixXd& jacobian_matrix, int jac_row,Eigen::MatrixXd& weighting_matrix)
  {
    Eigen::Vector3d pose_r = base_pose;
    double OMGE_ = 7.2921151467E-5;
    double CLIGHT_ = 299792458.0;

    /* satellite position*/
    Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
    double var_u2m = getVarofPr(dd_measurement.u_master_SV.elevation);

    Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
    double var_u2i = getVarofPr(dd_measurement.u_iSV.elevation);

    Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
    double var_r2m = getVarofPr(dd_measurement.r_master_SV.elevation);

    Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
    double var_r2i = getVarofPr(dd_measurement.r_iSV.elevation);

    double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
    est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

    double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);
    est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

    double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
    est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*u_pose(1)-u_pose_m(1)*u_pose(0))/CLIGHT_;
    

    double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);
    est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*u_pose(1)-u_pose_i(1)*u_pose(0))/CLIGHT_;

    // jacobian_matrix(jac_row, 0) = (u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m;

    // jacobian_matrix(jac_row, 1) = (u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m;

    // jacobian_matrix(jac_row, 2) = (u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m;

    double factor = 1;
    jacobian_matrix(jac_row, 0) = factor * ((u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m);

    jacobian_matrix(jac_row, 1) = factor * ((u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m);

    jacobian_matrix(jac_row, 2) = factor * ((u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m);

    for(int i = 3; i < jacobian_matrix.cols(); i++)
    {
      jacobian_matrix(jac_row, i) = 0;
    }
    
    #if use_fixed_cov_ar
    weighting_matrix(jac_row,jac_row) = 1.0 / (pow(0.4, 2));
    #else
    weighting_matrix(jac_row,jac_row) = 1.0 / ((var_u2m + var_u2i + var_r2m + var_r2i)/4.0); 
    #endif 

  }

  /* get carrier-phase double-differenced Jacobian matrix */
  void getCpDDJacobian(Eigen::Vector3d u_pose, Eigen::Vector3d base_pose, DDMeasurement dd_measurement, Eigen::MatrixXd& jacobian_matrix, int jac_row, int carrier_phase_index, Eigen::MatrixXd& weighting_matrix)
  {
    Eigen::Vector3d pose_r = base_pose;
    double OMGE_ = 7.2921151467E-5;
    double CLIGHT_ = 299792458.0;

    /* satellite position*/
    Eigen::Vector3d u_pose_m(dd_measurement.u_master_SV.sat_pos_x, dd_measurement.u_master_SV.sat_pos_y,dd_measurement.u_master_SV.sat_pos_z);
    double var_u2m = getVarofCp(dd_measurement.u_master_SV.elevation);

    Eigen::Vector3d u_pose_i(dd_measurement.u_iSV.sat_pos_x, dd_measurement.u_iSV.sat_pos_y,dd_measurement.u_iSV.sat_pos_z); 
    double var_u2i = getVarofCp(dd_measurement.u_iSV.elevation);

    Eigen::Vector3d r_pose_m(dd_measurement.r_master_SV.sat_pos_x, dd_measurement.r_master_SV.sat_pos_y,dd_measurement.r_master_SV.sat_pos_z);
    double var_r2m = getVarofCp(dd_measurement.r_master_SV.elevation);

    Eigen::Vector3d r_pose_i(dd_measurement.r_iSV.sat_pos_x, dd_measurement.r_iSV.sat_pos_y,dd_measurement.r_iSV.sat_pos_z); 
    double var_r2i = getVarofCp(dd_measurement.r_iSV.elevation);

    // double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
    // double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);

    // double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
    // double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);

    double est_p_r2m = getDistanceFrom2Points(pose_r, r_pose_m);
    est_p_r2m = est_p_r2m + OMGE_ * (r_pose_m(0)*pose_r(1)-r_pose_m(1)*pose_r(0))/CLIGHT_;

    double est_p_r2i = getDistanceFrom2Points(pose_r, r_pose_i);
    est_p_r2i = est_p_r2i + OMGE_ * (r_pose_i(0)*pose_r(1)-r_pose_i(1)*pose_r(0))/CLIGHT_;

    double est_p_u2m = getDistanceFrom2Points(u_pose, u_pose_m);
    est_p_u2m = est_p_u2m + OMGE_ * (u_pose_m(0)*u_pose(1)-u_pose_m(1)*u_pose(0))/CLIGHT_;
    

    double est_p_u2i = getDistanceFrom2Points(u_pose, u_pose_i);
    est_p_u2i = est_p_u2i + OMGE_ * (u_pose_i(0)*u_pose(1)-u_pose_i(1)*u_pose(0))/CLIGHT_;

    double factor = 1;
    // LOG(INFO)<<"jacobian_matrix.cols()"<<jacobian_matrix.cols();
    // LOG(INFO)<<"jacobian_matrix.rows()"<<jacobian_matrix.rows();
    // LOG(INFO)<< "jac_row-> "<<jac_row;
    // LOG(INFO)<<"carrier_phase_index->" << carrier_phase_index;
    // LOG(INFO)<<"jacobian_matrix(jac_row, 0)->" << jacobian_matrix(jac_row, 0);
    jacobian_matrix(jac_row, 0) = factor * ((u_pose_i(0) - u_pose(0))/(est_p_u2i) - (u_pose_m(0) - u_pose(0))/est_p_u2m);

    jacobian_matrix(jac_row, 1) = factor * ((u_pose_i(1) - u_pose(1))/(est_p_u2i) - (u_pose_m(1) - u_pose(1))/est_p_u2m);

    jacobian_matrix(jac_row, 2) = factor * ((u_pose_i(2) - u_pose(2))/(est_p_u2i) - (u_pose_m(2) - u_pose(2))/est_p_u2m);

    for(int i = 3; i < jacobian_matrix.cols(); i++)
    {
      jacobian_matrix(jac_row, i) = 0; 
    }
    
    jacobian_matrix(jac_row, 3 + carrier_phase_index) = dd_measurement.r_master_SV.lamda;

    // weighting_matrix(jac_row,jac_row) = 1.0/(var_u2m + var_u2i + var_r2m + var_r2i);
    weighting_matrix(jac_row,jac_row) = 1.0/(pow(0.004, 2));

    #if use_fixed_cov_ar
    
    weighting_matrix(jac_row,jac_row) = 1.0/(pow(0.004, 2));
    #else
    weighting_matrix(jac_row,jac_row) = 1.0 / ((var_u2m + var_u2i + var_r2m + var_r2i)/4.0); 
    #endif 
  }

  /* get variance for carrier-phase from a single satellite based on elevation */
  double getVarofCp(double ele)
  {
    double a = 0.01; 
    double b = 0.01;
    double c = 0;
    double d = 0;
    double square_sigma = pow(a,2) + pow(b,2) / pow(sin(ele * D2R), 2) + pow(c,2) + pow(d,2);
    // return 0.004;
    return square_sigma;
  }

  /* get variance for carrier-phase from a single satellite based on elevation/SNR */
  double getVarofCp_ele_SNR(nlosExclusion::GNSS_Raw single_sat_data)
  {
    Eigen::Matrix<double,4,1> parameters;
    parameters<<50.0, 30.0, 30.0, 10.0; // loosely coupled 
    // parameters<<50.0, 30.0, 20.0, 30.0; // loosely coupled 
    double snr_1 = parameters(0); // T = 50
    double snr_A = parameters(1); // A = 30
    double snr_a = parameters(2);// a = 30
    double snr_0 = parameters(3); // F = 10
    double snr_R = single_sat_data.snr;
    double elR = single_sat_data.elevation;
    // if(elR<15) elR = 30;
    if(elR<15)
    {
      LOG(INFO) << "satellite elevation -> " << elR;
    }
    double q_R_1 = 1 / (pow(( sin(elR * 3.1415926/180.0 )),2));
    double q_R_2 = pow(10,(-(snr_R - snr_1) / snr_a));
    double q_R_3 = (((snr_A / (pow(10,(-(snr_0 - snr_1) / snr_a))) - 1) / (snr_0 - snr_1)) * (snr_R - snr_1) + 1);
    double q_R = q_R_1* (q_R_2 * q_R_3);
    double var =(1.0/float(q_R)); // uncertainty: cofactor_[i] larger, larger uncertainty

    double a = 0.01; 
    double b = 0.01;
    double c = 0;
    double d = 0;
    double var_ele = pow(a,2) + pow(b,2) / pow(sin(elR * D2R), 2) + pow(c,2) + pow(d,2);
    // return 0.004;
    #if useEleVar
      return sqrt(var_ele);
    #else 
      0.0001 * sqrt(1/var);
    #endif 
  }

  /* get variance for pseudorange from a single satellite based on elevation */
  double getVarofPr(double ele)
  {
    double p_c_ratio = 1; // 10
    double a = 3 * p_c_ratio; 
    double b = 3 * p_c_ratio;
    double c = 0;
    double d = 0;
    double square_sigma = pow(a,2) + pow(b,2) / pow(sin(ele * D2R), 2) + pow(c,2) + pow(d,2);
    // return 0.4;
    // std::cout<<"var of pseudorange -> "<<var<<std::endl;
    return square_sigma;
  }

    /* get variance for pseudorange from a single satellite based on elevation/SNR */
  double getVarofpr_ele_SNR(nlosExclusion::GNSS_Raw single_sat_data)
  {
    Eigen::Matrix<double,4,1> parameters;
    parameters<<50.0, 30.0, 30.0, 10.0; // loosely coupled 
    // parameters<<50.0, 30.0, 20.0, 30.0; // loosely coupled 
    double snr_1 = parameters(0); // T = 50
    double snr_A = parameters(1); // A = 30
    double snr_a = parameters(2);// a = 30
    double snr_0 = parameters(3); // F = 10
    double snr_R = single_sat_data.snr;
    double elR = single_sat_data.elevation;
    // if(elR<15) elR = 30;
    if(elR<15)
    {
      // LOG(INFO) << "satellite elevation -> " << elR;
    }
    double q_R_1 = 1 / (pow(( sin(elR * 3.1415926/180.0 )),2));
    double q_R_2 = pow(10,(-(snr_R - snr_1) / snr_a));
    double q_R_3 = (((snr_A / (pow(10,(-(snr_0 - snr_1) / snr_a))) - 1) / (snr_0 - snr_1)) * (snr_R - snr_1) + 1);
    double q_R = q_R_1* (q_R_2 * q_R_3);
    double var =(1.0/float(q_R)); // uncertainty: cofactor_[i] larger, larger uncertainty

    double p_c_ratio = 1; // 10
    double a = 3 * p_c_ratio; 
    double b = 3 * p_c_ratio;
    double c = 0;
    double d = 0;
    double var_ele = pow(a,2) + pow(b,2) / pow(sin(elR * D2R), 2) + pow(c,2) + pow(d,2);
    #if useEleVar
      return sqrt(var_ele);
    #else 
      return 1 * sqrt(1/var);
    #endif 
  }


};

#endif // POSE_SYSTEM_HPP
