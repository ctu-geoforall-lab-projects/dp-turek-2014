#!/usr/bin/env python

#%module
#% description: Bundle block djustment.
#%end

#%option  G_OPT_M_DIR
#% key: bingo_data_dir
#% description: Path to the file of Bingo-F protocol.
#% required: yes
#%end

#%option G_OPT_M_DIR
#% key: protocol_dir
#% description: Directory where the protocols will be generated.
#% required: yes
#%end

#%flag
#% key: f
#% description: Adjust as free network.
#%end


import os
import sys
import numpy as np
import timeit

import grass.script as grass


from plots import PlotScene
from readpar import InputData
from parsers import BingoParser, FabioDataParser
from helmert import HelmertTransform, Test
from relative_or import RelativeOrientation
from affine import Affine_Fit
from bba import RunBundleBlockAdjutment, RotMatBLUH, GetBLUHAngles

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

def solvePnP(in_dt):

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()
    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()
    import cv2


    for ph in phs.itervalues():
        img_pts = []
        obj_pts = []
        pts = ph.GetPoints()
        i = 0
        for pt in pts.itervalues():
            gcp = pt.GetGcp()[0]
            phpts = pt.GetPhotoPoints()
            phpt = phpts[ph]
            img_pts.append([phpt[0], phpt[1]])
            obj_pts.append([gcp[0], gcp[1], gcp[2]])
            i += 1
            #if i == 4:
            #    break
        img_pts = np.array([img_pts])
        obj_pts = np.array([obj_pts])
    
        ret, r, t = cv2.solvePnP(obj_pts, img_pts, cam_m, distor, flags=cv2.EPNP)
        ret, r, t = cv2.solvePnP(obj_pts, img_pts, cam_m, distor, rvec=r, tvec=t, useExtrinsicGuess=True)

        r = cv2.Rodrigues(r)[0]
        t = t[:,0]
        c = - r.T.dot(t)
        eo = ph.GetEO()
        c_eo = eo[:,3]

        p = np.zeros((3, 4))
        p[:, :3] = r
        p[:, 3] = t
        #print "photo"
        #print c_eo - c
        ph.SetP(p)

    return

    #cv2.solvePnP(objectPoints, imagePoints, cameraMatrix, distCoeffs[, rvec[, tvec[, useExtrinsicGuess[, flags]]]])  retval, rvec, tvec

def main():
    in_dt = InputData(BingoParser(options["bingo_data_dir"]))

    # "/home/ostepok/Dev/GRASS/diplomka/data/Bingo_project/withEO_3GCP"
    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()
    
    cam_m, distor = in_dt.GetCameras().values()[0].GetParams()
    rel_or_phs = _readRelPhotoCombnationsFromBingoRelax("/home/ostepok/Dev/GRASS/diplomka/data/Bingo_project/withEO_3GCP/relax.lis")
    #[[1, 2], [1, 3],[2, 3]]
    #rel_or_phs = [[1, 2], [1, 3],[2, 3]]
    
    # computes relative orientations of pairs, merges them into common relative system
    ros = RelativeOrientation(in_dt, rel_or_phs, '/home/ostepok/Desktop/rel_protocol.txt')

    np.set_printoptions(precision=3, suppress=True)

    HelmertTransform(pts, phs)
    #print "ahoj"
    PlotScene(pts, phs)

    #solvePnP(in_dt)
    #PlotRelOrs(pts, phs)
    
    #Test(ros, in_dt, cam_m, distor, rel_or_phs)
    
    skip_pts = []
    skip_phs = []

    ph_errs, pt_errs = RunBundleBlockAdjutment(in_dt, skip_phs, skip_pts, '/home/ostepok/Desktop/protocol.txt',
        '/home/ostepok/Desktop/protocol_measurements.txt')

   
    PlotScene(pts, phs)
    return 

if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
