#!/bin/bash
docker container ls -a -f name=graph1 | grep graph1$ > /dev/null

if [ $? == 0 ]
then
	docker container start graph1
	docker exec -it graph1 /bin/bash

else
	xhost +
	SHARED_DOCKER_DIR=/root/graph1/shared_dir
	SHARED_HOST_DIR=$HOME/shared_dir
	if [ ! -d "$SHARED_HOST_DIR" ];then
	    mkdir $SHARED_HOST_DIR
	fi
	#docker run --gpus all -it -v /home/hrz/project/ros/melodic/LvisamTest:/root/LvisamTest  -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=unix$DISPLAY -e GDK_SCALE -e GDK_DPI_SCALE -e ROS_MASTER_URI=http://172.17.0.2:8899 -e ROS_HOSTNAME=172.17.0.2 --cpus 0.5 --cpuset-cpus=2 --name melodic melodic-gpu /bin/bash 
	docker run -it   -v /tmp/.X11-unix:/tmp/.X11-unix -v $SHARED_HOST_DIR:$SHARED_DOCKER_DIR  -e DISPLAY=unix$DISPLAY -e GDK_SCALE -e GDK_DPI_SCALE -e ROS_MASTER_URI=http://172.17.0.2:11311 -e ROS_HOSTNAME=172.17.0.2 --name graph1 ros:GraphGNSSLib /bin/bash 
	sudo chmod 777 $SHARED_HOST_DIR
fi
