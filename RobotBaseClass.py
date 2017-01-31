'''
Created on Apr 8, 2010

@author: cbpark
'''
from scipy.linalg import *
from numpy import *
from MyMathLib import *
from numpy import array as ar

# http://www.gebrauchtroboter.com/rob_angeb/R065_Kuka_KR200_KRC1/daten_arbeitsbereich.html
class RobotBaseClass():
    def __init__(self, EEF_Pos=array([0,0,0]), EEF_Ori=eye(3)):
        pass
    
    def CheckJointLimit(self, joint_angle_deg):
        for i in range(self.DOF):
            if (joint_angle_deg[i] < self.joint_min_angle_deg[i]) | (joint_angle_deg[i] > self.joint_max_angle_deg[i]):
                return False
            
        return True
    def InvKinNum_Safe(self, T, ja):
        ja_backup = ja.copy()
        
        ja, ret = InvKinNum(self,T,ja_backup)
        if ret == False:
            ja = ja_backup.copy()
            return ja, False
            
        jl = self.CheckJointLimit(r2d(ja))
        if jl == False:
            ja = ja_backup.copy()
            return ja, False
        return ja, True
    

        
    def ComputeWorkspace(self):
        ja_home = d2r(self.joint_angle_home_deg)
        ja = ja_home.copy()
        EEF_T0 = FwdKin (self, ja)
        R0, P0 = FromSE3(EEF_T0)
        
        x0 = P0[0]
        y0 = P0[1]
        z0 = P0[2]
        x, y, z = x0, y0, z0
        dx, dy, dz = 0.05,  0.05, 0.05
        
        WSUIB = []
        WSUOB = []
        WSLIB = []
        WSLOB = []

        x = x0
        z = z0
        while(True):
            print 1, x,y,z
            P = array([x,y,z])
            T = ToSE3(R0, P)
            ja, ret = self.InvKinNum_Safe(T,ja_home)
            if ret == False:
                break
            
            x = x0
            while(True):
#                print 2, x,y,z
                P = array([x,y,z])
                T = ToSE3(R0, P)
                ja, ret = self.InvKinNum_Safe(T,ja_home)
                if ret == False:
                    WSUIB.append(array([x_old,y,z]))
                    
                    break
                x_old = x
                x = x - dx

            x = x0
            while(True):
#                print 3, x,y,z
            
                P = array([x,y,z])
                T = ToSE3(R0, P)
                ja, ret = self.InvKinNum_Safe(T,ja_home)
                if ret == False:
                    WSUOB.append(array([x_old,y,z]))
                    break
                x_old = x
                x = x + dx
            
            
            z = z + dz

        x = x0
        z = z0
        ja = ja_home
        
        while(True):
            P = array([x,y,z])
            print 4, x,y,z
           
            # I.K.
            T = ToSE3(R0, P)
            ja, ret = self.InvKinNum_Safe(T,ja_home)
            if ret == False:
                break
            
            x = x0
            while(True):
                print 5, x,y,z

                P = array([x,y,z])
                T = ToSE3(R0, P)
                ja, ret = self.InvKinNum_Safe(T,ja_home)
                if ret == False:
                    WSLIB.append(array([x_old,y,z]))
                    break
                x_old = x
                x = x - dx

            x = x0
            while(True):
                print 6, x,y,z
                
                P = array([x,y,z])
                T = ToSE3(R0, P)
                ja, ret = self.InvKinNum_Safe(T,ja_home)
                if ret == False:
                    WSLOB.append(array([x_old,y,z]))
                    break
                x_old = x
                x = x + dx
            z = z - dz
    
        dprint(WSLIB)
        dprint(WSLOB)
        dprint(WSUIB)
        dprint(WSUOB)
        