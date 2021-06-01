#!/usr/bin/env python  
import roslib
import rospy
import math
import tf
import geometry_msgs.msg
import turtlesim.srv
from nav_msgs.msg import Odometry
from geometry_msgs.msg import Pose, PoseWithCovarianceStamped, Point, Quaternion, Twist, PoseStamped
from nav_msgs.msg import Path

if __name__ == '__main__':
    rospy.init_node('CartoTF2Pose')
    poseCartoPub = rospy.Publisher('/poseCarto', Odometry, queue_size=10)  # customerized odometry from ENU
    cartoPathpub = rospy.Publisher('/cartoPath', Path, queue_size=100)
    cartoPath_ = Path()
    CartoPose_ = PoseStamped()
    CartoPose_.header.frame_id = 'velodyne'
    listener = tf.TransformListener()
    rate = rospy.Rate(int(1))
    print 'hshshsh'
    poseCarto = Odometry()
    while not rospy.is_shutdown():
        try:
            (trans,rot) = listener.lookupTransform('/map', '/base_link', rospy.Time(0))
        except (tf.LookupException, tf.ConnectivityException, tf.ExtrapolationException):
            continue
        print 'carto trans',trans
        print 'carto rotation',rot
        print '------'
        poseCarto.header.frame_id = 'map'
        poseCarto.pose.pose.position.x = trans[0]
        poseCarto.pose.pose.position.y = trans[1]
        poseCarto.pose.pose.position.z = trans[2]
        poseCarto.pose.pose.orientation.x = rot[0]
        poseCarto.pose.pose.orientation.y = rot[1]
        poseCarto.pose.pose.orientation.z = rot[2]
        poseCarto.pose.pose.orientation.w = rot[3]
        CartoPose_.pose = poseCarto.pose.pose
        cartoPath_.poses.append(CartoPose_)
        cartoPath_.header.frame_id = 'velodyne'
        print len(cartoPath_.poses)
        cartoPathpub.publish(cartoPath_)
        poseCartoPub.publish(poseCarto)
        rate.sleep()
