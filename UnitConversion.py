'''
Created on Mar 19, 2010

@author: cbpark
'''
from numpy import *

def deg2rad(deg):
    rad = deg*pi/180
    return rad

def rad2deg(rad):
    deg = rad/pi*180
    return deg

def d2r(deg):
    rad = deg*pi/180
    return rad

def r2d(rad):
    deg = rad/pi*180
    return deg

def rpm2radps(rpm):
    radps = rpm*(2*pi/60)
    return radps

def radps2rpm(rpm):
    rpm = radps/(2*pi/60)
    return rpm

def Nm2lbft(Nm):

    lbft = Nm / 1.356
    return lbft

def Nm2lbin(Nm):

    lbin = Nm / 0.1130
    return lbin

def lbft2Nm(lbft):

    Nm  = lbft * 1.356
    return lbft

def lbin2Nm(lbin):

    Nm = lbin * 0.1130
    return lbin




