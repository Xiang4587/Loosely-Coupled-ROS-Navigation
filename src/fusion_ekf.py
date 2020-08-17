#!/usr/bin/env python

import rospy
from std_msgs.msg import String
from geometry_msgs.msg import PoseStamped
from geometry_msgs.msg import Twist, Point, Quaternion
from sensor_msgs.msg import Imu
from dbw_mkz_msgs.msg import WheelSpeedReport
from nmea_msgs.msg import GpsInfo
from novatel_gps_msgs.msg import Gpgga
import matplotlib.pyplot as plt
import numpy as np
import math
from numpy.linalg import inv
import tf
from tf.transformations import euler_from_quaternion, quaternion_from_euler, euler_matrix, euler_from_matrix
from scipy.linalg import block_diag
from rotations import Quaternion, skew_symmetric

def NomalizeAngle(angle):
    if angle <= -(math.pi):
        angle += 2*(math.pi)
    if angle > math.pi:
        angle -= 2*(math.pi)
    return angle

def so3exp(rho):
    angle = np.linalg.norm(rho)

    # Near phi==0, use first order Taylor expansion
    if np.abs(angle) < 1e-8:
        skew_rho = np.array([[0, -rho[2], rho[1]],
                             [rho[2], 0, -rho[0]],
                             [-rho[1], rho[0], 0]])
        return np.identity(3) + skew_rho

    axis = rho / angle
    skew_axis = np.array([[0, -axis[2], axis[1]],
                          [axis[2], 0, -axis[0]],
                          [-axis[1], axis[0], 0]])
    s = np.sin(angle)
    c = np.cos(angle)
    return c * np.eye(3) + (1 - c) * np.outer(axis, axis) + s * skew_axis

class kalman_class():

    def __init__(self, x, P):

        self.previous_x = np.array(x)
        self.previous_P = np.diag(P)
        self.estimated_x = np.zeros((15,1)) 
        self.estimated_P = np.zeros((15,15)) 
        self.predicted_x = np.zeros((15,1)) 
        self.predicted_P = np.zeros((15,15))
        self.A = np.zeros((15,15))
        self.Phi = np.zeros((15,15))
        self.Qmat = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0], # north
                               [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0], # east
                               [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0], # down
                               [0,0,0,0.1**2,0,0,0,0,0,0,0,0,0,0,0], # vnorth
                               [0,0,0,0,0.1**2,0,0,0,0,0,0,0,0,0,0], # veast
                               [0,0,0,0,0,0.1**2,0,0,0,0,0,0,0,0,0], # vdown
                               [0,0,0,0,0,0,0.0175**2,0,0,0,0,0,0,0,0], # roll
                               [0,0,0,0,0,0,0,0.0175**2,0,0,0,0,0,0,0], # pitch
                               [0,0,0,0,0,0,0,0,0.1571**2,0,0,0,0,0,0], # yaw
                               [0,0,0,0,0,0,0,0,0,4.0e-04,0,0,0,0,0], # bias_a1
                               [0,0,0,0,0,0,0,0,0,0,4.0e-04,0,0,0,0], # bias_a2
                               [0,0,0,0,0,0,0,0,0,0,0,4.0e-04,0,0,0], # bias_a3
                               [0,0,0,0,0,0,0,0,0,0,0,0,4.3633e-04,0,0], # bias_w1
                               [0,0,0,0,0,0,0,0,0,0,0,0,0,4.3633e-04,0], # bias_w2
                               [0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.3633e-04]]) # bias_w3
        
        # Longitude, Latitude, Height constriant
        self.constantH = -61.7
        self.vary_h = 0
        self.br_lc = tf.TransformBroadcaster()
    
    def Predict(self, count, time, T, input):
        earth_gravity = 9.807
        euler = tf.transformations.euler_from_quaternion(input.quatern)
        roll = self.previous_x[6][0]
        pitch = self.previous_x[7][0]
        yaw = self.previous_x[8][0]
        #roll = euler[0]
        #pitch = euler[1]
        #yaw = euler[2]
        print("yaw from IMU:", yaw)
        print("yaw from GPS:", input.gps_yaw)
        cyaw = math.cos(yaw)
        syaw = math.sin(yaw)

        ## ----------------------------
        if count == 0:
            print("Initial roll,pitch,yaw: ", roll, pitch, yaw)
            print("Initial state:", self.previous_x[:3])
            self.predicted_x = self.previous_x
            self.predicted_P = self.previous_P

        if count < 1:  ## old_time:0, new_time:0
            pass
        else:
            if count == 1 or count == 2:        ## old_time:0, new_time:15..... ==> T:15... ; old_time:15..., new_time:15... ==> T:0
                T = 0.1
            
            print("T_adjusted:", T)
            print("Time_rightnow:", time)
            print("previous state: ", self.previous_x[:3])

            # # Eular to rotation matrix(static axis xyz)
            R_prev = euler_matrix(roll, pitch, yaw, 'sxyz') #syxz
            R_prev = np.delete(R_prev, 3, 0)
            R_prev = np.delete(R_prev, 3, 1) #3x3
            
            ## Accelerate propagate
            imu_acc_input = np.array([[input.acce1-self.previous_x[9][0]],
                                      [input.acce2-self.previous_x[10][0]],
                                      [input.acce3-self.previous_x[11][0]]])

            acc = R_prev.dot(imu_acc_input) + np.array([[0],[0],[earth_gravity]])

            # ## Velocity propagate
            self.predicted_x[3] = self.previous_x[3] + acc[0]*T
            self.predicted_x[4] = self.previous_x[4] + acc[1]*T
            self.predicted_x[5] = self.previous_x[5] + acc[2]*T
            
            ## Position propagate
            self.predicted_x[0] = self.previous_x[0] + self.predicted_x[3]*T + 1/2*acc[0]*(T**2)
            self.predicted_x[1] = self.previous_x[1] + self.predicted_x[4]*T + 1/2*acc[1]*(T**2)
            self.predicted_x[2] = self.previous_x[2] + self.predicted_x[5]*T + 1/2*acc[2]*(T**2)
            
            ## Attitude propogate
            omega = np.array([input.rate1-self.previous_x[12][0], input.rate2-self.previous_x[13][0], input.rate3-self.previous_x[14][0]])
            #print("omega:", omega)
            R_predict = R_prev.dot(so3exp(omega*T))
            self.predicted_x[6], self.predicted_x[7], self.predicted_x[8] = euler_from_matrix(R_predict, 'sxyz')
            
            ## Linearize Motion Model and compute Jacobians
            self.A[0][3] = 1
            self.A[1][4] = 1
            self.A[2][5] = 1
            
            self.A[3][6] = (input.acce3 - input.bias_a3)*syaw
            self.A[3][7] = (input.acce3 - input.bias_a3)*cyaw
            self.A[3][8] = -(input.acce1 - input.bias_a1)*syaw - (input.acce2 - input.bias_a2)*cyaw + (input.acce3 - input.bias_a3)*(-syaw*pitch + cyaw*roll)
            self.A[3][9] = -cyaw;
            self.A[3][10] = syaw;
            self.A[3][11] = -(cyaw*pitch + syaw*roll)

            self.A[4][6] = -(input.acce3 - input.bias_a3)*cyaw
            self.A[4][7] = (input.acce3 - input.bias_a3)*syaw
            self.A[4][8] = (input.acce1 - input.bias_a1)*cyaw - (input.acce2 - input.bias_a2)*syaw + (input.acce3 - input.bias_a3)*(cyaw*pitch + syaw*roll)
            self.A[4][9] = -syaw;
            self.A[4][10] = -cyaw;
            self.A[4][11] = -(syaw*pitch - cyaw*roll)
            
            self.A[5][6] = (input.acce2 - input.bias_a2)
            self.A[5][7] = -(input.acce1 - input.bias_a1)
            self.A[5][9] = pitch
            self.A[5][10] = -roll
            self.A[5][11] = -1

            self.A[6][7] = input.rate3 - input.bias_g3;
            self.A[6][12] = -1;
            self.A[6][14] = -pitch;

            self.A[7][6] = -(input.rate3 - input.bias_g3);
            self.A[7][13] = -1;
            self.A[7][14] = roll;

            self.A[8][6] = input.rate2 - input.bias_g2;
            self.A[8][13] = -roll;
            self.A[8][14] = -1;

            self.Phi = np.eye(15) + self.A*T #15x15

            self.predicted_P = self.Phi.dot(self.previous_P.dot(self.Phi.T)) + self.Qmat*T

            ##
            #print("predicted state: ", self.predicted_x[:6])
            print("predicted yaw: ", self.predicted_x[8][0])



    def Update(self, time, T, input):
        GNSS_input = False
        Odom_input = False
        #
        ## GNSS update
        # ==============================================================================
        H_gnss = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                            [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                            [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                            [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                            [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                            [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
                            [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0]])  # 7x15

        H_gnss_redu = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        	                    [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
        	                    [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                                [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0]]) # 6x15
        
        # R will be changed and ZIHR (vdown ~= 0)
        speed1 = math.sqrt(self.previous_x[3]**2 + self.previous_x[4]**2) ## only vector vnorth + vector veast
        cond1 = speed1 - 2 # 2 is the threshold
        cond2 = 10*(abs(input.rate3)-0.1)
        factor1 = (1 + 5*math.exp(-cond1/2)) * (1 + 5*math.exp(-cond2/1))
        r_gnss_yaw = factor1 * 0.0175**2 
        #print("r_gnss_yaw", r_gnss_yaw)
        R_gnss = np.array([[0.5563**2 ,0,0,0,0,0,0],
                            [0,0.5563**2 ,0,0,0,0,0],
                            [0,0,2.2253**2 ,0,0,0,0],
                            [0,0,0,0.1       ,0,0,0],
                            [0,0,0,0,0.1       ,0,0],
                            [0,0,0,0,0,0.1       ,0],
                            [0,0,0,0,0,0,r_gnss_yaw]]) # 7x7

        R_gnss_redu = np.array([[0.5563**2 ,0,0,0,0,0],
                                [0,0.5563**2 ,0,0,0,0],
                                [0,0,2.2253**2 ,0,0,0],
                                [0,0,0,0.1       ,0,0],
                                [0,0,0,0,0.1       ,0],
                                [0,0,0,0,0,0.1       ]]) # 6x6
        
        #print("position_z:", self.predicted_x[2],"constantH:", self.constantH, "vary_h: ", self.vary_h, "postion_z:", (self.predicted_x[2]-self.constantH)/0.1)
        measurement = np.array([[input.position_x], [input.position_y], [input.position_z], 
                                [input.vnorth], [input.veast], [0], # vdown = 0
                                [input.gps_yaw]])
        output_predict = np.array([[self.previous_x[0][0]], [self.previous_x[1][0]], [self.previous_x[2][0]],
                                   [self.previous_x[3][0]], [self.previous_x[4][0]], [self.previous_x[5][0]],
                                   [self.previous_x[8][0]]])

        measurement_redu = np.array([[input.position_x], [input.position_y], [input.position_z - self.vary_h*T], 
                                [input.vnorth], [input.veast], [0]])
        output_predict_redu = np.array([[self.previous_x[0][0]], [self.previous_x[1][0]], [self.previous_x[2][0]],
                           [self.previous_x[3][0]], [self.previous_x[4][0]], [self.previous_x[5][0]]])
       
        
        if input.gpstime >= time-T and input.gpstime <= time:
            GNSS_input = True
        	## Caulculate residual
            if input.gps_qual == 4 and input.num_sats > 15:
                P_z = H_gnss.dot(self.predicted_P.dot(H_gnss.T)) + R_gnss #7x7
                #P_z = H_gnss_redu.dot(self.predicted_P.dot(H_gnss_redu.T)) + R_gnss_redu #6x6

                K = self.predicted_P.dot(H_gnss.T.dot(inv(P_z))) # 15x7
                #K = self.predicted_P.dot(H_gnss_redu.T.dot(inv(P_z))) #15x6

                residual =  measurement - output_predict # 7x1
                #residual =  measurement_redu - output_predict_redu # 6x1

                ## Nomalize yaw angle
                #residual[6][0] = NomalizeAngle(residual[6][0])
                
                self.estimated_x = self.predicted_x + K.dot(residual)

                self.estimated_P = self.predicted_P - K.dot(H_gnss.dot(self.predicted_P))
                #self.estimated_P = self.predicted_P - K.dot(H_gnss_redu.dot(self.predicted_P))
            else:
                measurement_con = np.array([[self.previous_x[0][0]], [self.previous_x[1][0]], [self.previous_x[2][0]]]) # 3x1
                output_predict_con = np.array([[self.predicted_x[0][0]], [self.predicted_x[1][0]], [self.predicted_x[2][0]]]) # 3x1
                R_ins = np.array([[0.5563**2,0,0], #0.5563**2
                              [0,0.5563**2,0],     #0.5563**2
                              [0,0,2.2253**2 ]]) # 3x3
                H_ins = np.array([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                              [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                              [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0]]) #3x15
                P_z_con = H_ins.dot(self.predicted_P.dot(H_ins.T)) + R_ins
                K_con = self.predicted_P.dot(H_ins.T.dot(inv(P_z_con)))

                residual_con = measurement_con - output_predict_con #3x1

                self.estimated_x = self.predicted_x + K_con.dot(residual_con)
                self.estimated_P = self.predicted_P - K_con.dot(H_ins.dot(self.predicted_P))

        else:
            GNSS_input = False
            self.estimated_x = self.predicted_x
            self.estimated_P = self.predicted_P
        

        self.vary_h = (self.previous_x[2][0] - self.constantH)/0.06
        
        #print("Estimate_GNSS:", self.estimated_x[:6])
        print("Estimate_GNSS_yaw:", self.estimated_x[8][0])
        #
        ## Odometry update
        # ==============================================================================
        #if input.odomtime >= time-T and input.odomtime <= time:
        #Odom_input = True
        
        INS_speed = math.sqrt(self.estimated_x[3]**2 + self.estimated_x[4]**2 + self.estimated_x[5]**2)
        #print("odom_speed:", input.odom_speed, "velocity:", INS_speed)
        residual_speed = input.odom_speed - INS_speed # 1x1
        H_speed = np.array([[0,0,0,self.estimated_x[3][0]/INS_speed,self.estimated_x[4][0]/INS_speed,self.estimated_x[5][0]/INS_speed,0,0,0,0,0,0,0,0,0]]) # 1x15
        if GNSS_input == True:
            if input.odom_speed < 2:
                std_speed = 0.5
            else:
                std_speed = 5/input.odom_speed
        else:
            std_speed = 0.01
        R_speed =  std_speed**2
        P_z_speed = H_speed.dot(self.estimated_P.dot(H_speed.T)) + R_speed
        K_speed = self.estimated_P.dot(H_speed.T.dot(inv(P_z_speed))) # 15x1
        self.estimated_x = self.estimated_x + K_speed.dot(residual_speed)
        self.estimated_P = self.estimated_P - K_speed.dot(H_speed.dot(self.estimated_P))
        
        print("Estimate_Odom_yaw:", self.estimated_x[8][0])

            

    def publish_message(self, count, time, caller_obj, TFbroadcaster):
        # publisher
        Pub = rospy.Publisher('/looselycouple/kf', PoseStamped, queue_size=10)
        state_new = PoseStamped()
        # Publish kalma.x and kalman.P
        x = self.estimated_x[0][0]
        y = self.estimated_x[1][0]
        z = self.estimated_x[2][0]

        # extracting Quaternion from the homogeneous
        print("Estimate:", self.estimated_x[:3])
        print("roll:", self.estimated_x[6], "pitch:", self.estimated_x[7], "yaw:", self.estimated_x[8])
        #print("Unormalize:", self.estimated_x[8])
        print("Before changing previous state yaw 11111:", self.estimated_x[6:9])

        if count < 3:
            self.previous_x = self.predicted_x
            self.previous_P = self.predicted_P
        else:
            self.previous_x = self.estimated_x
            self.previous_P = self.estimated_P
        
        quatern = quaternion_from_euler(self.estimated_x[6][0], self.estimated_x[7][0], NomalizeAngle(self.estimated_x[8][0])) ## This function will convert the original angle to half
        print("Before changing previous state yaw 22222:", self.estimated_x[6:9])
        
        time = rospy.Time.from_sec(time)
        # Constructing the message
        state_new.header.stamp = time
        state_new.header.frame_id = 'LC'
        state_new.pose.position = Point(x, y, z)
        state_new.pose.orientation = Quaternion(*quatern)

        # publishing the message
        Pub.publish(state_new)
            
        ## Publish trajectory
        TFbroadcaster.sendTransform((x,y,z),  ## z-34
                        (quatern[0], quatern[1], quatern[2], quatern[3]),
                         time,
                         "LC",
                         "map")  
        
        print("After changing previous state yaw:", self.previous_x[6:9])
        print("-------------------------------------------")

class Kalman_caller(object):
    def __init__(self):
        self.imuSub = rospy.Subscriber("/xsens/imu/data", Imu, self.imu_cb)
        self.gpsSub = rospy.Subscriber("/xsens/gps/pose", PoseStamped, self.gps_cb)
        self.gpsInfoSub = rospy.Subscriber("/xsens/gps/info", GpsInfo, self.gpsInfo_cb)
        self.odoSub = rospy.Subscriber("/vehicle/wheel_speed_report", WheelSpeedReport, self.odo_cb)
        self.novatelSub = rospy.Subscriber("/novatel/gpgga", Gpgga, self.novatel_cb)
        self.n = 0
        # self.odom_covariance = np.empty((4, 4), dtype=int)
        # self.imu_covariance = np.empty((2, 2), dtype=int)
        # self.I_see_something = False
        self.time_stamp = None
        
        ## GPS
        self.gpstime = None
        self.position_x = None
        self.position_y = None
        self.position_z = None
        self.vnorth = None
        self.veast = None
        self.velocity = None
        self.track = None
        self.gps_ori = [0,0,0,0]
        self.gps_pose_ori = None
        self.gps_yaw = None

        ## IMU
        self.IMUtime = None
        self.acce1 = None
        self.acce2 = None
        self.acce3 = None
        self.rate1 = None
        self.rate2 = None
        self.rate3 = None
        self.orientation_cov = np.empty((3,3))
        self.bias_a1 = 0.0004
        self.bias_a2 = 0.0004
        self.bias_a3 = 0.0004
        self.bias_g1 = 0.00043633
        self.bias_g2 = 0.00043633
        self.bias_g3 = 0.00043633
        self.ori_x = 0
        self.ori_y = 0
        self.ori_z = 0
        self.ori_w = 0
        self.quatern = [0,0,0,0]
        
        # Odometry
        self.odomtime = None
        self.odom_time_stamp = None
        self.odom_speed = None
        
        # NoVatel
        self.novateltime = None
        self.num_sats = None
        self.gps_qual = None

    def imu_cb(self, imu_msg):
        self.IMUtime = imu_msg.header.stamp.secs + imu_msg.header.stamp.nsecs*10**(-9)
        x_rot = np.array([[1, 0, 0],
                         [0, math.cos(-math.pi/2), math.sin(-math.pi/2)],
                         [0, -math.sin(-math.pi/2), math.cos(-math.pi/2)]])
        z_rot = np.array([[math.cos(-math.pi/2), math.sin(-math.pi/2), 0],
                         [-math.sin(-math.pi/2), math.cos(-math.pi/2), 0],
                         [0, 0, 1]])
        convert = x_rot.dot(z_rot)
        
        
        self.acce1 = imu_msg.linear_acceleration.x
        self.acce2 = imu_msg.linear_acceleration.y
        self.acce3 = -imu_msg.linear_acceleration.z
        self.rate1 = imu_msg.angular_velocity.x
        self.rate2 = imu_msg.angular_velocity.y
        self.rate3 = imu_msg.angular_velocity.z
        self.ori_x = imu_msg.orientation.x
        self.ori_y = imu_msg.orientation.y
        self.ori_z = imu_msg.orientation.z
        self.ori_w = imu_msg.orientation.w
        self.quatern[0] = self.ori_x
        self.quatern[1] = self.ori_y
        self.quatern[2] = self.ori_z
        self.quatern[3] = self.ori_w


    def gps_cb(self, gps_msg):
        self.gpstime = gps_msg.header.stamp.secs + gps_msg.header.stamp.nsecs*10**(-9)
        #print("GPS:", self.gpstime)
        self.position_x = gps_msg.pose.position.x
        self.position_y = gps_msg.pose.position.y
        self.position_z = gps_msg.pose.position.z
        self.gps_ori[0] = gps_msg.pose.orientation.x
        self.gps_ori[1] = gps_msg.pose.orientation.y
        self.gps_ori[2] = gps_msg.pose.orientation.z
        self.gps_ori[3] = gps_msg.pose.orientation.w
        self.gps_pose_ori = tf.transformations.euler_from_quaternion(self.gps_ori)
        self.gps_yaw = self.gps_pose_ori[2]
       


    def gpsInfo_cb(self, gpsInfo_msg):
        self.vnorth = gpsInfo_msg.vnorth
        self.veast = gpsInfo_msg.veast
        self.velocity = gpsInfo_msg.velocity
        course = gpsInfo_msg.course
        self.track = course
        #self.track = NomalizeAngle(course)  ## 0~2pi
        #print("yaw from GPS:", self.track)
        #print("Orignal radian: ", gpsInfo_msg.course)
        #print("Nomalize:", self.track ,"Nomalize radian: ", NomalizeAngle(self.track))


    def odo_cb(self, odo_msg):
        self.odomtime = odo_msg.header.stamp.secs + odo_msg.header.stamp.nsecs*10**(-9)
        odom_front_left = odo_msg.front_left
        odom_front_right = odo_msg.front_right
        odom_rear_left = odo_msg.rear_left
        odom_rear_right = odo_msg.rear_right
        self.odom_speed = 0.3393*0.5*(odom_rear_left + odom_rear_right) # wheel radius: 0.3393 m

    def novatel_cb(self, novatel_msg):
        self.novateltime = novatel_msg.header.stamp.secs + novatel_msg.header.stamp.nsecs*10**(-9)
        self.num_sats = novatel_msg.num_sats
        self.gps_qual = novatel_msg.gps_qual
