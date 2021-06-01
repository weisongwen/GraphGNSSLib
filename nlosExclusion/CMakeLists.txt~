cmake_minimum_required(VERSION 2.8.3)
project(nlosExclusion)

find_package(catkin REQUIRED COMPONENTS
  message_filters
  nodelet
  roscpp
  pcl_conversions
  sensor_msgs
  tf2
  tf2_ros
  tf2_sensor_msgs
  pcl_ros
  rospy
  std_msgs
  message_generation
)
add_message_files(
  FILES
  GNSS_Raw.msg
  GNSS_Raw_Array.msg
  DOP.msg
  Satellite_Info.msg
  exclusionSatNum.msg
  Error.msg
)
generate_messages(
  DEPENDENCIES 
   std_msgs
  )

catkin_package(


  CATKIN_DEPENDS message_runtime roscpp message_filters nodelet sensor_msgs tf2 tf2_ros
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

find_package(PCL 1.7 REQUIRED)
include_directories(include ${PCL_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS} ${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})


install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
install(FILES nodelets.xml
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)

install(DIRECTORY launch
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
)
