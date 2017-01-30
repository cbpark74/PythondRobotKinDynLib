'''
Created on Mar 7, 2010

@author: cbpark
'''

import sys
sys.path.append(r'C:\UserSource\PythonCodes\MyPythonLib')

from scipy.linalg import *
from numpy import *
from MyMathLib import *

from RobotBaseClass import *

class SetRobot(RobotBaseClass):
    def __init__(self, EEF_Pos=array([0,0,0]), EEF_Ori=RotY_rad(d2r(-90))):

############################## Enter Robot Data below
         
        self.Name = 'Robot' 
        self.DOF = DOF = 6
        
        # unit mm
        h1 = 740
        h2 = 1075
        h3 = 250
        d1 = 305
        d2 = 1275 
        d3 = 240
        
        scale = 0.001 # mm -> m
        
        w = array([
             [0,0,1],
             [0,1,0],
             [0,1,0],
             [1,0,0],
             [0,1,0],
             [1,0,0]
             ])
             
        p = scale*array([
                   [0,0,0],
                   [d1,0,h1],
                   [d1,0,h1+h2],
                   [d1,0,h1+h2+h3],
                   [d1+d2,0,h1+h2+h3],
                   [d1+d2+d3,0,h1+h2+h3],
                   ])
        

        self.gravity = [0, 0, -9.81]
        
        # 
        self.m = [200, 200, 200, 100, 100, 200] # mass in kg
        self.r = array([ # mass center wrt i-th frame 
                        array([0,0,0]),
                        array([0,0,0]),
                        array([0,0,0]),
                        array([0,0,0]),
                        array([0,0,0]),
                        array([0,0,0]),
                        ]) 
        
        self.joint_max_angle_deg = array([150,90,120,120,255,120])# Max. Angle, in degrees
        self.joint_min_angle_deg = -self.joint_max_angle_deg
        self.joint_angle_home_deg = array([0,0,0,0,0,0]) # Home Position
  
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
        
        dprint(P0)
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
    q_zero = array([0, 0, 0, 0, 0, 0]) * pi / 180
    ToolRelPos = array([0.0, 0, 0.0])*0.5
    R = eye(3)
    Robot = SetRobot(ToolRelPos)
    Robot.SetEndEffectorKinParameter(ToolRelPos, R)
    
    q_arb = array([10, 10, 10, 90, 0, 90]) * pi / 180
    q_home = array([90, 0, 0, 0, 0, 0]) * pi / 180
    q = q_home
    dq = zeros(Robot.DOF)
    ddq = zeros(Robot.DOF)
    
    T = FwdKin(Robot, q)
    nprint (T)
    
    Force = array([0, 0, 0])
    Torque = array([0, 0, 0])
    F_ext = hstack([Torque, Force])
    nprint (F_ext)
    
    gravity = array([0, 0, -9.81])
    jointTorque = ComputeTorque(Robot, q, dq, ddq, F_ext, gravity)
    nprint (jointTorque)
