#!/usr/bin/env python
import sys
import math
from sympy import *
from bba import xp_col_eq, yp_col_eq, _rotMatBLUH


# functions of this file use SimPy to create python source code of collinearity derivatives (see file bba_derivs.py")

def GetSymbolicRotMatBLUH():

    ph = Symbol('ph')
    om = Symbol('om')
    ka = Symbol('ka')

    rotm = [[0,0,0], [0,0,0], [0,0,0]] 

    rm = _rotMatBLUH(ph, om, ka, rotm, sin, cos)

    return rotm

def CreateSymbolicColEquations():

    f = Symbol('f')
    xi_c = Symbol('xi_c')
    yi_c = Symbol('yi_c')
    
    xo_c = Symbol('xo_c')
    yo_c = Symbol('yo_c')
    zo_c = Symbol('zo_c')
    
    xo_pt = Symbol('xo_pt')
    yo_pt = Symbol('yo_pt')
    zo_pt = Symbol('zo_pt')

    r1_dis = Symbol('r1_dis')

    r = GetSymbolicRotMatBLUH()

    xc = xp_col_eq(xi_c, yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r, r1_dis)
    yc = yp_col_eq(xi_c, yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r, r1_dis)

    return xc, yc

def GenerateDistSource(d, name, fd):

    fd.write("%s={\n" % name)
    for k, v in d.iteritems():
        fd.write("        '%s' : %s,\n" % (k, v))
    fd.write("}\n")

def Derivate(path):
    x_derivs = { "ph" : "PhiDerX",
                 "om" : "OmegaDerX",
                 "ka" : "KappaDerX",
                 "xo_c" : "XEODerX",
                 "yo_c" : "YEODerX",
                 "zo_c" : "ZEODerX",
                 "xo_pt" : "XPtDerX",
                 "yo_pt" : "YPtDerX",
                 "zo_pt" : "ZPtDerX",
                 "xi_c" : "XCenterDerX",
                 "yi_c" : "YCenterDerX",
                 "f" : "FocalDerX",
                 #"r1_dis" : "Rad1DistorX"
               }

    y_derivs = { "ph" : "PhiDerY",
                 "om" : "OmegaDerY",
                 "ka" : "KappaDerY",
                 "xo_c" : "XEODerY",
                 "yo_c" : "YEODerY",
                 "zo_c" : "ZEODerY",
                 "xo_pt" : "XPtDerY",
                 "yo_pt" : "YPtDerY",
                 "zo_pt" : "ZPtDerY",
                 "xi_c" : "XCenterDerY",
                 "yi_c" : "YCenterDerY",
                 "f" : "FocalDerY",
                 #"r1_dis" : "Rad1DistorY"
               }
    params =  "(" + ', '.join(x_derivs.keys()) + ", **kwargs):"

    xc, yc = CreateSymbolicColEquations()

    derivs = {}

    der_fd = open(path, "w")

    der_fd.write("#automaticaly generated by derivate.py\n\n")
    der_fd.write("from math import sin, cos\n\n")
    for un, func_name_x in x_derivs.iteritems():
        func_name_y = y_derivs[un]
        der_fd.write("def " + func_name_x + params + "\n")
        der_fd.write("    return " + str(diff(xc, un, 1)) + "\n")

        der_fd.write("def " + func_name_y + params + "\n")
        der_fd.write("    return " + str(diff(yc, un, 1)) + "\n\n")

    GenerateDistSource(x_derivs, "x_derivs",  der_fd)

    GenerateDistSource(y_derivs, "y_derivs",  der_fd)

    der_fd.close()
    return
if __name__ == '__main__':

    path = str(sys.argv)[0]
    Derivate(path)