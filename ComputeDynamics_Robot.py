'''
Created on Mar 7, 2010

@author: cbpark
'''

import sys

from scipy.linalg import *
from numpy import *
from MyMathLib import *

from RobotBaseClass import *

class SetRobot(RobotBaseClass):
    def __init__(self, EEF_Pos=array([0,0,0]), EEF_Ori=RotY_rad(d2r(-90))):

############################## Enter Robot Data below
         
        self.Name = 'Robot' 
        self.DOF = DOF = 2
        
        # unit mm
        l1 = 1000
        l2 = 1000
        
        scale = 0.001 # mm -> m
 
        # Direction of Joint Rotation
        w = array([
             [0,1,0],
             [0,1,0],
             ])

        # Position of Joint     
        p = scale*array([
                   [0,0,0],
                   [l1,0,0],
                   ])
       
        # Gravity Vector    
        self.gravity = [0, 0, -9.81]
        
        # 
        self.m = [1, 1, ] # link weight in kg
        self.r = array([ # mass center wrt i-th frame 
                        array([0.5,0,0]),
                        array([0.5,0,0]),
                        ]) 
        
        self.joint_max_angle_deg = array([360, 360])# Max. Angle, in degrees
        self.joint_min_angle_deg = -self.joint_max_angle_deg
        self.joint_angle_home_deg = array([0,0]) # Home Position
  
############## Commmon Part
        self.scale = scale
        

        self.p = p
        self.w = w

        
        self.M = []
        
        self.M.append(eye(4))
        self.M[0][0:3, 3] = p[0]

        for i in range(1, DOF, 1):
            self.M.append(eye(4))
            self.M[i][0:3, 3] = p[i] - p[i - 1]
        
        P0 = self.p[DOF-1]+EEF_Pos
        
        R0 = EEF_Ori
        self.P0 = P0

        self.M.append(eye(4))

        self.M[DOF][0:3, 0:3] = R0
        self.M[DOF][0:3, 3] = P0 - p[DOF - 1]
        
        self.M0 = eye(4)
        self.M0[0:3, 0:3] = R0
        self.M0[0:3, 3] = P0
        
        zero_v = [0, 0, 0]
        v = []; self.A = []; self.S = [];
        for i in range(0, DOF, 1):
            v.append(cross(p[i], w[i]))
            
            A = hstack([w[i], v[i]])
            S = hstack([w[i], zero_v]) 
            self.A.append(A)
            self.S.append(S)
        
        self.R0 = RotY_rad(d2r(90))
    
        self.InvK = eye(self.DOF);

        self.I = [];# link inertia matrix wrt i-th frame
        for i in range(DOF):
            self.I.append(zeros((3, 3)))

    def SetEndEffectorKinParameter(self, EEF_Pos=array([0,0,0]), EEF_Ori=RotY_rad(d2r(-90)) ):
        DOF = self.DOF
        P0 = self.p[DOF-1]+EEF_Pos
        R0 = EEF_Ori
        self.P0 = P0

        self.M.append(eye(4))

        dprint(P0)
        dprint(self.p[DOF - 1])

        self.M[DOF][0:3, 0:3] = R0
        self.M[DOF][0:3, 3] = P0 - self.p[DOF - 1]
        
        self.M0 = eye(4)
        self.M0[0:3, 0:3] = R0
        self.M0[0:3, 3] = P0
        
                
if __name__ == "__main__":
    Deg2Rad = pi / 180.0;
    Rad2Deg = 180.0 / pi;
    
    q_zero = array([0, 0]) * Deg2Rad;
    ToolRelPos = array([1.0, 0, 0.0]);
    Robot = SetRobot(ToolRelPos)
    
    q_deg = array([0, 90])
    q = q_deg * Deg2Rad;
    dq = array([0, 0])
    ddq = array([0, 0])
    
    ToolPose = FwdKin(Robot, q)
    TCP_Pos = ToolPose[0:3, 3]
    
    nprint (ToolPose)
    nprint (TCP_Pos)
    
    EE_Force = array([0, 0, 0])
    EE_Torque = array([0, 0, 0])
    F_ext = hstack([EE_Torque, EE_Force])
    nprint (F_ext)
    
    gravity = array([0, 0, -9.81])
    jointTorque = ComputeTorque(Robot, q, dq, ddq, F_ext, gravity)
    nprint (jointTorque)
