cmake_minimum_required(VERSION 2.8.3)
project(global_fusion)

set(CMAKE_BUILD_TYPE "Debug")

# set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS "-std=c++14")
# C++ 14 is required
add_definitions(-std=c++14)
set(CMAKE_CXX_STANDARD 14)

#-DEIGEN_USE_MKL_ALL")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(Eigen3)
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

find_package(Gflags)
find_package(Glog)
include_directories(${GFLAGS_INCLUDE_DIRS} ${GLOG_INCLUDE_DIRS})
include_directories(include)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  nlosExclusion
)

find_package(Ceres REQUIRED)

add_subdirectory(./ThirdParty/GeographicLib/)
add_subdirectory(./RTKLIB/)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ./ThirdParty/GeographicLib/include/
)


catkin_package()

# add rtklib
add_library(convgpx RTKLIB/src/convgpx.c)
add_library(convkml RTKLIB/src/convkml.c)
add_library(convrnx RTKLIB/src/convrnx.c)
add_library(datum RTKLIB/src/datum.c)
add_library(download RTKLIB/src/download.c)
add_library(ephemeris RTKLIB/src/ephemeris.c)
add_library(geoid RTKLIB/src/geoid.c)
add_library(gis RTKLIB/src/gis.c)
add_library(ionex RTKLIB/src/ionex.c)
add_library(lambda RTKLIB/src/lambda.cpp)
add_dependencies(lambda nlosExclusion_gencpp)
add_library(options RTKLIB/src/options.c)
add_library(pntpos RTKLIB/src/pntpos.cpp)
add_dependencies(pntpos nlosExclusion_gencpp)
add_library(postpos RTKLIB/src/postpos.cpp)
add_dependencies(postpos nlosExclusion_gencpp)
add_library(ppp_ar RTKLIB/src/ppp_ar.c)
add_library(ppp RTKLIB/src/ppp.c)
add_library(preceph RTKLIB/src/preceph.c)
# add_library(qzslex RTKLIB/src/qzslex.c)
add_library(rcvraw RTKLIB/src/rcvraw.c)
add_library(rinex RTKLIB/src/rinex.c)
add_library(rtcm RTKLIB/src/rtcm.c)
add_library(rtcm2 RTKLIB/src/rtcm2.c)
add_library(rtcm3 RTKLIB/src/rtcm3.c)
add_library(rtcm3e RTKLIB/src/rtcm3e.c)
add_library(rtkcmn RTKLIB/src/rtkcmn.c)
add_library(rtkpos RTKLIB/src/rtkpos.cpp)
add_dependencies(rtkpos nlosExclusion_gencpp)
add_library(rtksvr RTKLIB/src/rtksvr.c)
add_library(sbas RTKLIB/src/sbas.c)
add_library(solution RTKLIB/src/solution.c)
add_library(stream RTKLIB/src/stream.c)
add_library(streamsvr RTKLIB/src/streamsvr.c)
add_library(tides RTKLIB/src/tides.c)
add_library(tle RTKLIB/src/tle.c)

add_executable(testRTKLIBNode 
	src/testRTKLIBNode.cpp)
target_link_libraries(testRTKLIBNode ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 



add_executable(gnss_preprocessor_node 
	src/gnss_preprocessor/gnss_preprocessor.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(gnss_preprocessor_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

add_executable(rosbag_generator_node 
	src/gnss_preprocessor/rosbag_generator.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(rosbag_generator_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

add_executable(psr_spp_node 
	src/gnss_estimator/psr_spp.cpp)
target_link_libraries(psr_spp_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 


add_executable(psr_spp_dyna_auto_node 
	src/gnss_estimator/psr_spp_dyna_auto.cpp)
target_link_libraries(psr_spp_dyna_auto_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 


add_executable(psr_doppler_fusion_node 
	src/gnss_estimator/psr_doppler_fusion.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(psr_doppler_fusion_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

add_executable(psr_doppler_car_rtk_node 
	src/gnss_estimator/psr_doppler_car_rtk.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(psr_doppler_car_rtk_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

add_executable(psr_doppler_car_rtk_dyna_node 
	src/gnss_estimator/psr_doppler_car_rtk_dyna.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(psr_doppler_car_rtk_dyna_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

add_executable(evaluate_gps_solution_node 
	src/evo/evaluate_gps_solution.cpp
	./RTKLIB/src/rtkpos.cpp
	./RTKLIB/src/postpos.cpp
	./RTKLIB/src/ppp.c)
target_link_libraries(evaluate_gps_solution_node ${catkin_LIBRARIES} ${CERES_LIBRARIES} 
						convgpx convkml convrnx datum download ephemeris geoid gis ionex lambda
						options pntpos postpos ppp_ar ppp preceph rcvraw rinex
						rtcm rtcm2 rtcm3 rtcm3e rtkcmn rtkpos rtksvr sbas solution
						stream streamsvr tides tle
						libGeographiccc) 

