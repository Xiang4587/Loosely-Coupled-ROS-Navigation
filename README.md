# Loosely-Coupled-ROS-Navigation
This repository contains the code to fuse GPS/IMU/Wheel odometry for the vehicle localization in ROS. It's available for users to their own sensors such as GPS, IMU and wheel odometry depending on vehicles.
![image](https://github.com/Xiang4587/Loosely-Coupled-ROS-Navigation/blob/master/doc/Fusion_in_ROS.png)
## reference
P. D. Froves, "Principles of GNSS, Inertial, and Multisensor Integrated Navigation Systems, [Book review]," _IEEE_ Arospace and Electronic Systems Magazine, vol. 30, no. 2, pp. 26-27, 2015

## Installation & Prerequies
1. Install ROS (kinetic), rospy, and ros drivers for sensors depending on yours. We perform this package on Autoware platform.\
```
sudo apt-get update -y
```
```
sudo apt-get install -y python-rospy
```
2. Install the following packages, `numpy`, `scipy`, e.g. with the pip command
```
pip install numpy scipy 
```
3. Clone this repo in catkin_ws
```
cd catkin_ws/src
```
```
git clone 
```
