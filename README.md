# Loosely-Coupled-ROS-Navigation
This repository contains the code to fuse GPS/IMU/Wheel odometry for the vehicle localization in ROS. It's available for users to their own sensors such as GPS, IMU and wheel odometry depending on vehicles.
![image](https://github.com/Xiang4587/Loosely-Coupled-ROS-Navigation/blob/master/doc/Fusion_in_ROS.png)
## reference
P. D. Froves, "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems, [Book review]," _IEEE_ Arospace and Electronic Systems Magazine, vol. 30, no. 2, pp. 26-27, 2015

## Installation & Prerequies
1. Install ROS (kinetic), rospy, and ros drivers for sensors depending on yours. We perform this package on Autoware platform.\
```
$ sudo apt-get update -y
```
```
$ sudo apt-get install -y python-rospy
```
2. Install the following packages, `numpy`, `scipy`, e.g. with the pip command
```
$ pip install numpy scipy 
```
3. Clone this repo in catkin_ws
```
$ cd catkin_ws/src
```
```
$ git clone git@github.com:Xiang4587/Loosely-Coupled-ROS-Navigation.git
```

## Usage
Run this package
```
$ roscore
```
Online sensing GPS/IMU/wheel odometry data or play rosbag file
```
$ cd catkin_ws/src/Loosly-Coupled-ROS-Navigation/src
$ chmod +x LC_main.py
```
```
$ rosrun Loosly-Coupled-ROS-Navigation LC_main.py
```

## Description
The state of the extended Kalman filter for GPS/IMU/Odometry fusion includes position, velocity, attitude:
[P_n, P_e, P_d, V_n, V_e, V_d, roll, pitch, yaw, bias_g1, bias_g2, bias_g3, bias_w1, bias_w2, bias_w3]
