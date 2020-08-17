#!/usr/bin/env python

import rospy
from std_msgs.msg import String
from geometry_msgs.msg import PoseStamped
from sensor_msgs.msg import Imu
import matplotlib.pyplot as plt
import fusion_ekf
import math
import tf

## Specify initial state from loosely couple 
# =====================================================
init_n = 6762.97459455  
init_e = -8181.38892808
init_d = -59.4381523223
init_vn = -0.44428169012944596
init_ve = -1.0356160590378831
init_vd = 7.930784577471797
init_roll = 0.017436878      ## Origin roll, pitch, yaw angle
init_pitch = -0.973799989 
init_yaw = 2.727936519 
init_ba1 = 8.798919297464678e-06
init_ba2 = 3.365448668524238e-07
init_ba3 = 4.076829264994726e-05
init_bg1 = -2.175173428157313e-07
init_bg2 = 8.906784458343999e-06
init_bg3 = -5.259844050797483e-08

x = [[init_n], [init_e], [init_d], 
     [init_vn], [init_ve], [init_vd],
     [init_roll], [init_pitch], [init_yaw],
     [init_ba1], [init_ba2], [init_ba3],
     [init_bg1], [init_bg2], [init_bg3]]

P = [0.551715385464264,
     0.5517116877346624,
     1.9064471206229054,
     0.28788978939418647,
     0.2877699007732788,
     0.014122916385081376,
     0.017649483451635056,
     0.01764004200956554,
     0.03851406597962453,
     0.000400065055115232,
     0.0004000651039640964,
     0.00040006410399308226,
     0.0004364097810290437,
     0.000436409733250568,
     0.0004364049001403606] 

            
def main():
    rospy.init_node('LooselyCoupled')
    tf_broadcast = tf.TransformBroadcaster()
    
    listener = fusion_ekf.Kalman_caller()
    kalman = fusion_ekf.kalman_class(x, P)
    
    Rate = rospy.Rate(8.0) ## Depends on sensors sampling rate ex. GPS: 4 Hz; IMU: 100- Hz; Wheel_speed_report: 100+ Hz 
    try:
        count = 0
        old_time = rospy.Time().now().to_sec()
        while not rospy.is_shutdown():

            # time step for prediction
            new_time = rospy.Time().now().to_sec()
            T = new_time - old_time
            print("old_time:", old_time)
            print("new_time:", new_time)
            print("T:", T)
            print("count:", count)


            kalman.Predict(count, new_time, T, listener)   ## Dead reckoning
            
            # Update (the first time has no input)
            if count > 0 and listener.position_z != None:
                kalman.Update(new_time, T, listener)

            # # Publish LC
            kalman.publish_message(count, old_time, listener, tf_broadcast)

            old_time = new_time
            count = count + 1
            Rate.sleep()

    except rospy.ROSInterruptException:
        pass
    rospy.spin()
    

if __name__=='__main__':
    main()
