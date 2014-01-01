#!/usr/bin/env python

import math
from sympy import *
from bba import xp_col_eq, yp_col_eq, _rotMatBLUH


def GetSymbolicRotMatBLUH():

    ph = Symbol('ph')
    om = Symbol('om')
    ka = Symbol('ka')

    rotm = [[0,0,0], [0,0,0], [0,0,0]] 

    rm = _rotMatBLUH(ph, om, ka, rotm)

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

    r = GetSymbolicRotMatBLUH()

    xc = xp_col_eq(xi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r)
    yc = yp_col_eq(yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r)

    return xc, yc

def GenerateDistSource(d, name, fd):

    fd.write("%s={\n" % name)
    for k, v in d.iteritems():
        fd.write("        '%s' : %s,\n" % (k, v))
    fd.write("}\n")

def Derivate():
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
               }

    params =  "(" + ', '.join(x_derivs.keys()) + ", **kwargs):"

    xc, yc = CreateSymbolicColEquations()

    derivs = {}

    der_fd = open("bba_derivs.py", "w")

    der_fd.write("#automaticaly generated by derivate.py\n\n")
    der_fd.write("from math import sin, cos\n\n")
    for un, func_name_x in x_derivs.iteritems():
        func_name_y = y_derivs[un]
        der_fd.write("def " + func_name_x + params + "\n")
        der_fd.write("    return " + str(diff(xc, un, 1)) + "\n")

        der_fd.write("def " + func_name_y + params + "\n")
        der_fd.write("    return " + str(diff(yc, un, 1)) + "\n\n")

        print un
        print diff(xc, un, 1)
        print diff(yc, un, 1)

    GenerateDistSource(x_derivs, "x_derivs",  der_fd)

    GenerateDistSource(y_derivs, "y_derivs",  der_fd)

    der_fd.close()

    return
if __name__ == '__main__':
    Derivate()