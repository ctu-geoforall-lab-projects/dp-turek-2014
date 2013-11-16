#!/usr/bin/env python

#%module
#% description: Bundle block djustment.
#%end

#%option G_OPT_F_INPUT
#% key: pointsfile
#% required: yes
#%end

#%option G_OPT_F_INTPUT
#% key: camfile
#% required: yes
#%end

import os
import sys
import numpy as np
from copy import deepcopy

from collections import namedtuple
from readpar import BingoParser, InputData
from helmert import HelmertTransform, Test
from relative_or import RelativeOrientation
from affine import Affine_Fit
import grass.script as grass
from math import atan2, sqrt, pi, sin, cos
from bba import RunBundleBlockAdjutment, RotMatBLUH, GetBLUHAngles

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def _readRelPhotoCombnationsFromBingoRelax(relax_file):
    fd = open(relax_file, "r")

    dt = []
    while 1:
        line =  fd.readline()
        if not line:
            break

        l = line.split()
        if len(l) < 10:
            continue

        if l[0]  != 'Block' and l[1] != "2":
            continue
        dt.append([int(l[3]), int(l[5])])

    fd.close()

    return dt

def plotPts(pts_list, lines_list):
    fig = plt.figure()

    i_col = 1
    cols = len(pts_list)
    for i, pts in  enumerate(pts_list):
        if pts is not None:
            ax = fig.add_subplot(1, cols, i_col, projection='3d')
            ax.scatter(pts[:,0], pts[:,1], pts[:,2], c='g', marker='^')
            i_col += 1

        for lns in lines_list:
            if lns is not None:
                for ln in lns:
                    ln = np.reshape(ln, (-1, 3))
                    ax.plot(ln[:,0], ln[:,1], ln[:,2])

    plt.show()

def ComputeDiffs(pts, phs):

    gcp_pts = []
    ro_pts = []

    pt_errs = []
    no_coords = []
    for pt in pts.itervalues():
        if pt.GetCoords() is not None:
            pt_errs.append([pt.GetId(), 
                         abs(np.linalg.norm(pt.GetCoords() - pt.GetResultCoords()))])
        else:
            no_coords.append(pt.GetId())

    pt_errs = np.array(pt_errs)
    print pt_errs
    print np.average(pt_errs[:,1])

    skip_pts = []
    for err in pt_errs:
        if err[1] > 0.3:
            skip_pts.append(int(err[0]))


    ph_errs = []
    for ph in phs.itervalues():
        if ph.GetEO() is not None:
            eo = ph.GetEO()
            res_eo = ph.GetResultExtOr()

            angles = GetBLUHAngles(eo[:,:3])
            #TODO

            res_coords = res_eo[:,3]
            res_angles = np.array(res_eo[1]) / pi * 180
            #print "angles %d" % ph.GetId()
            #print angles
            #print res_angles

            ph_errs.append([ph.GetId(), 
                            abs(np.linalg.norm(eo[:,3] - res_coords))] )
        else:
            no_coords.append(ph.GetId())

    ph_errs = np.array(ph_errs)

    skip_phs = []
    for err in ph_errs:
        if err[1] > 0.5:
            skip_phs.append(int(err[0]))

    print skip_phs
    print skip_pts
    print ph_errs
    print np.average(ph_errs[:,1])

    #return skip_phs, skip_pts
    return ph_errs, pt_errs
    #print np.average(ph_errs[:,2:], axis=1)

def PlotRelOrs(pts, phs, ro=None):

    res_pts = []
    plot_pts = []
    for pt in pts.itervalues():
        if pt.GetResultCoords() is not None and pt.GetCoords() is not None:
            res_pts.append(pt.GetResultCoords())
            plot_pts.append(pt.GetCoords())


    lines = []
    for ph in phs.itervalues():
        if ph.GetResultExtOr() is not None and ph.GetEO() is not None:
            res_eo = ph.GetResultExtOr()
            res_pts.append(res_eo[:, 3])
            plot_pts.append(ph.GetEO()[:,3])
            lines.append(GetCameraLines(ph.GetEO()))

            r = res_eo[:, :3]
            eo_res = np.hstack((r, np.array([ph.GetEO()[:,3]]).T))
            lines.append(GetCameraLines(eo_res, inv=True))
            
            lines.append(GetCameraLines(ph.GetEO()))

    plot_pts = np.array(plot_pts)
    res_pts = np.array(res_pts)

    plot_list = [plot_pts]
    #if ro:
    #    plot_list.append(ro)

    plotPts(plot_list, lines)

def GetCameraLines(eo, inv=False):

    size = 1
    apts = np.array([
                  [size, 0, 0],
                  [0, size, 0],
                  [0, 0, -size]])
    origin_mat = np.tile(eo[:,3], (apts.shape[0] , 1))
    axis_pts = (eo[:,:3].T.dot(apts.T)).T + origin_mat
    

    if not inv:
        apts = np.array([
                      [size, 0, 0],
                      [0, size, 0],
                      [0, 0, size]])
    
        origin_mat = np.tile(eo[:,3], (apts.shape[0] , 1))
        axis_pts = (eo[:,:3].T.dot(apts.T)).T + origin_mat
        
    """
    if not inv:
        dz = axis_pts[:, 2] - eo[2,3]
        axis_pts[:, 2]  = eo[2,3] - dz
    """

    return  np.hstack((origin_mat, axis_pts))

def main():
    

    path = "/home/ostepok/Dev/GRASS/diplomka/Bingo_project/withEO_3GCP"
    in_dt = InputData(BingoParser(path))

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()
    
    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()
    rel_or_phs = _readRelPhotoCombnationsFromBingoRelax("/home/ostepok/Dev/GRASS/diplomka/Bingo_project/withEO_3GCP/relax.lis")
    
    ros = RelativeOrientation(in_dt, cam_m, distor, rel_or_phs)
    np.set_printoptions(suppress=True)

    HelmertTransform(pts, phs)
    Test(ros, in_dt, cam_m, distor, rel_or_phs)
    PlotRelOrs(pts, phs)

    skip_phs, skip_pts = ComputeDiffs(pts, phs)
    
    ph_errs, pt_errs = ComputeDiffs(pts, phs)

    skip_pts = []
    skip_phs = []
    for i in range(3):
        RunBundleBlockAdjutment(in_dt, skip_phs, skip_pts)
        it_ph_errs, it_pt_errs = ComputeDiffs(pts, phs)
        print it_ph_errs[:,1]
        it_ph_errs = np.expand_dims(it_ph_errs[:,1], axis=0)
        it_pt_errs = np.expand_dims(it_pt_errs[:,1], axis=0)

        print it_ph_errs.shape
        print ph_errs.shape

        ph_errs = np.hstack((ph_errs, it_ph_errs.T)) 
        pt_errs = np.hstack((pt_errs, it_pt_errs.T)) 

    print ph_errs
    print pt_errs

    np.savetxt('/home/ostepok/Desktop/ph_res.csv', ph_errs, fmt='%.4f', delimiter=',')
    np.savetxt('/home/ostepok/Desktop/pt_res.csv', pt_errs, fmt='%.4f', delimiter=',')


    return 

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
