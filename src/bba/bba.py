import os
import sys
import numpy as np
from sympy import *

def RotMatBLUH(ph, om, ka):
    rotm = np.zeros((3,3))
    _rotMatBLUH(ph, om, ka, rotm)

    return rotm

def GetBLUHAngles(rm):

    ka = atan2(rm[0, 1], rm[1, 1])
    #print ka * 180 / pi
    ph = atan2(rm[2, 0], rm[2, 2])
    om = atan2(-rm[2, 1], sqrt(rm[1, 1] * rm[1, 1] + rm[0, 1] * rm[0, 1]))

    return np.array([ph, om, ka])

def _rotMatBLUH(ph, om, ka, rotm):

    so = sin(om)
    co = cos(om)
    sp = sin(ph)
    cp = cos(ph)
    sk = sin(ka)
    ck = cos(ka)

    rotm[0][0] = ck * cp + sk * so * sp
    rotm[1][0] = - sk * cp + ck * so * sp
    rotm[2][0] = co * sp

    rotm[0][1] = sk * co
    rotm[1][1] = ck * co
    rotm[2][1] = - so
    
    rotm[0][2] =  - ck * sp + sk * so * cp  
    rotm[1][2] =  sk * sp + ck * so * cp
    rotm[2][2] = co * cp

def GetSymbolicRotMatBLUH():

    ph = Symbol('ph')
    om = Symbol('om')
    ka = Symbol('ka')

    rotm = [[0,0,0], [0,0,0], [0,0,0]] 

    rm = _rotMatBLUH(ph, om, ka, rotm)

    return rotm

def CreateSymbolicColEquations():

    f = Symbol('f')
    x0 = Symbol('x0')
    y0 = Symbol('y0')
    
    xo_0 = Symbol('xo_0')
    yo_0 = Symbol('yo_0')
    zo_0 = Symbol('zo_0')
    
    xo_pt = Symbol('xo_pt')
    yo_pt = Symbol('yo_pt')
    zo_pt = Symbol('zo_pt')

    f = Symbol('f')
    f = Symbol('f')

    r = GetSymbolicRotMatBLUH()

    xc = xp_col_eq(x0, f, xo_0, yo_0, zo_0, xo_pt, yo_pt, zo_pt, r)
    yc = yp_col_eq(x0, f, xo_0, yo_0, zo_0, xo_pt, yo_pt, zo_pt, r)

    print diff(xc, 'om', 1)

def xp_col_eq(x0, f, xo_0, yo_0, zo_0 ,xo_pt, yo_pt, zo_pt, r):
    return x0 - f * (r[0][0]*(xo_pt-xo_0) + r[0][1]*(yo_pt-yo_0) + r[0][2]*(zo_pt-zo_0) ) / \
                    (r[2][0]*(xo_pt-xo_0) + r[2][1]*(yo_pt-yo_0) + r[2][2]*(zo_pt-zo_0) ) 

def yp_col_eq(y0, f, xo_0, yo_0, zo_0, xo_pt, yo_pt, zo_pt, r):
    return y0 - f * (r[1][0]*(xo_pt-xo_0) + r[1][1]*(yo_pt-yo_0) + r[1][2]*(zo_pt-zo_0) ) / \
                    (r[2][0]*(xo_pt-xo_0) + r[2][1]*(yo_pt-yo_0) + r[2][2]*(zo_pt-zo_0) )

def RunBundleBlockAdjutment(in_dt):
    
    CreateSymbolicColEquations()
    
    #print rotm[0][0]
    #print  diff(rotm[0][0], om, 1)

    print "ahoj"

    return


