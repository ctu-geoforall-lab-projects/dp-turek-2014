import os
import sys
import numpy as np
import scipy as sc


from plots import PlotRelOrs
import grass.script as grass
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from math import sqrt
import cv2
from copy import deepcopy

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


class RelOr:
    def __init__(self, id, ph1_id, ph2_id, pm, pt_ids, sps_pts, tie_pts):
        self.id = id
        self.ph1_id = ph1_id
        self.ph2_id = ph2_id
        self.pm = pm
        self.pt_ids = pt_ids
        self.sps_pts = sps_pts
        self.tie_pts = tie_pts
        self.scale = None

def AddRelOrToIndex(ph1_id, ph2_id, ro_id, ph_ids_to_ro_id):
    if ph1_id not in ph_ids_to_ro_id:
        ph_ids_to_ro_id[ph1_id] = {ph2_id : ro_id}
    else:
        ph_ids_to_ro_id[ph1_id][ph2_id] = ro_id

    if ph2_id  not in ph_ids_to_ro_id:
        ph_ids_to_ro_id[ph2_id] = {ph1_id : ro_id}
    else:
        ph_ids_to_ro_id[ph1_id][ph2_id] = ro_id


def RelativeOrientation(in_dt, cam_m, distor, rel_or_phs, protocol_path):
    ph_ids_to_ro_id = {}
    ros = []

    #rel_or_phs = [(598, 576), (597, 598)]  
    #rel_or_phs = [(576, 598), (597, 598)]  

    #rel_or_phs = [(598, 576), (598, 597)]
    #rel_or_phs = [(576, 598), (598, 597), (597, 583)]  


    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    for i_ro_id, ro_phs in enumerate(rel_or_phs):
        tie_point_ids, tie_pts = in_dt.GetTiePoints(ro_phs)
        sps_pts, P, p1, p2 = computeRelativeOrientation(tie_pts, 
                                                        cam_m, 
                                                        distor)

        tie_pts = np.hstack((p1[0, ...] , p2[0, ...]))

        ro = RelOr(i_ro_id, ro_phs[0], ro_phs[1], P, tie_point_ids, sps_pts, tie_pts)
        ros.append(ro)

        if i_ro_id == 0:
            pm1 = np.eye(3, 4)
            phs[ro_phs[0]].SetP(pm1)
            phs[ro_phs[1]].SetP(ro.pm)
            out_pts = ro.pt_ids
        
            ro.scale = 1

            for i_pt, pt_id in enumerate(ro.pt_ids):
                pts[pt_id].AddRelOrCoords(i_ro_id, sps_pts[i_pt][:3])
        else:
            if ph_ids_to_ro_id.has_key(ro_phs[0]):
                conn_ph = phs[ro_phs[0]]
            elif ph_ids_to_ro_id.has_key(ro_phs[1]):
                conn_ph = phs[ro_phs[1]]
    
            if ph_ids_to_ro_id.has_key(ro.ph1_id):
                ro_conn = ros[ph_ids_to_ro_id[ro.ph1_id].values()[0]]
            else:
                ro_conn = ros[ph_ids_to_ro_id[ro.ph2_id].values()[0]]

            eo_ph, eo_trans = ConnectRelOr(phs, ro, ro_conn)
            out_pts = compute3Dpoints(eo_ph, eo_trans, tie_pts[:,:2], tie_pts[:,2:])

            for i_pt, pt_add_id in enumerate(ro.pt_ids):
                pts[pt_add_id].AddRelOrCoords(i_ro_id, out_pts[i_pt])

        apts =  {}
        for pid in ro.pt_ids:
            apts[pts[pid]] = pts[pid]

        aphs =  {}
        for phid in ro_phs:
            aphs[phs[phid]] = phs[phid]


        if ro_phs == [553, 591]:
        #if ro_phs == [597, 585]: 
        # or ro_phs == [576, 598]:
            PlotRelOrs(apts, aphs, i_ro_id)
        
        AddRelOrToIndex(ro_phs[0], ro_phs[1], i_ro_id, ph_ids_to_ro_id)

    createProtocol(protocol_path, ros, cam_m, distor, in_dt)
        
    computeMissingPoints(in_dt)
    return ros, ph_ids_to_ro_id

def computeMissingPoints(in_dt):

    phs = in_dt.GetPhotos()
    pts = in_dt.GetPoints()

    for pt in pts.itervalues():
        if not pt.GetRelOrsCoords():

            ph_pts = pt.GetPhotoPoints()
            ph_pts_list = ph_pts.keys()
            found = False
            for i, ph_pt1 in enumerate(ph_pts_list[:-1]):
                if ph_pt1.GetP() is not None:
                    for ph_pt2 in ph_pts_list[i + 1:]:
                        if ph_pt2.GetP() is not None:

                            p1 = ph_pt1.GetP()
                            p2 = ph_pt2.GetP()

                            tie_pts1 = np.array([ph_pts[ph_pt1]])
                            tie_pts2 = np.array([ph_pts[ph_pt2]])

                            cam_m1, distor1 = ph_pt1.GetCamera().GetParams()
                            cam_m2, distor2 = ph_pt2.GetCamera().GetParams()

                            tp1_u = undistortPoints(tie_pts1, cam_m1, distor1)
                            tp2_u = undistortPoints(tie_pts2, cam_m1, distor1)

                            out_pts = compute3Dpoints(p1, p2, tp1_u, tp2_u)
                            print out_pts
                            pt.AddRelOrCoords(-1, out_pts[0])
                        found = True
                        break
                if found:
                    break

def ConnectRelOr(phs, add_rel_or, conn_rel_or):
    def GetPivotPhoto(rel_or1, rel_or2):
        if rel_or1.ph1_id == rel_or2.ph1_id:
            pivot_ph_id = rel_or1.ph1_id
        elif rel_or1.ph1_id == rel_or2.ph2_id:
            pivot_ph_id = rel_or1.ph1_id
        else:
            pivot_ph_id = rel_or1.ph2_id

        return pivot_ph_id

    pivot_ph_id = GetPivotPhoto(add_rel_or, conn_rel_or)
    pivot_ph = phs[pivot_ph_id]

    pm_ro = add_rel_or.pm
    ro_pt_ids = add_rel_or.pt_ids
    t_ro = pm_ro[:, 3]
    R_ro = pm_ro[:, :3]

    #ro_pts = add_rel_or.sps_pts[:,:3]
    if pivot_ph == add_rel_or.ph2_id:
        #ro_pts = R_ro.dot(ro_pts) + t_ro
        phpts1 = add_rel_or.tie_pts[:,:2]
        phpts2 = add_rel_or.tie_pts[:,2:]
        add_ph = phs[add_rel_or.ph1_id]
    else:
        phpts2 = add_rel_or.tie_pts[:,2:]
        phpts1 = add_rel_or.tie_pts[:,:2]
        add_ph = phs[add_rel_or.ph2_id]

    def computeScale(add_rel_or, conn_rel_or, phs):

        pivot_ph_id = GetPivotPhoto(add_rel_or, conn_rel_or)

        if add_rel_or.ph1_id == pivot_ph_id:
            eo_add_ro = np.eye(3, 4)
        else:
            eo_add_ro = add_rel_or.pm

        if conn_rel_or.ph1_id == pivot_ph_id:
            eo_conn_ro = np.eye(3, 4)
        else:
            eo_conn_ro = conn_rel_or.pm
        scales = []

        t_conn = eo_conn_ro[:, 3]
        R_conn = eo_conn_ro[:, :3] 

        t_add = eo_add_ro[:, 3]
        R_add = eo_add_ro[:, :3] 

        for i_pt_add, pt_add_id in enumerate(add_rel_or.pt_ids):
            if pt_add_id in conn_rel_or.pt_ids:
                i_pt_conn = conn_rel_or.pt_ids.index(pt_add_id)
                pt_conn = conn_rel_or.sps_pts[i_pt_conn, :3]
                pt_conn = R_conn.dot(pt_conn) + t_conn

                pt_add = add_rel_or.sps_pts[i_pt_add, :3]
                pt_add = R_add.dot(pt_add) + t_add

                scales.append(np.linalg.norm(pt_conn) / np.linalg.norm(pt_add))

        scales = np.array(scales)
        #if inv:
        #    print "invers"
        #    scales = 1/scales

        med = np.median(scales)
        print "scale"
        print med
        print np.std(scales)

        return med

    eo_ph = pivot_ph.GetP()
    t_ph = eo_ph[:, 3]
    R_ph = eo_ph[:, :3]

    scale = computeScale(add_rel_or, conn_rel_or, phs)

    add_rel_or.scale = scale * conn_rel_or.scale

    t_trans = R_ro.dot(t_ph) + add_rel_or.scale * t_ro
    R_trans = R_ro.dot(R_ph)
    
    eo_trans = np.hstack((R_trans, np.expand_dims(t_trans, axis=1)))
    add_ph.SetP(eo_trans)

    return eo_ph, eo_trans 

def undistortPoints(pts, cam_m, distor):
    tp = np.expand_dims(pts, axis=0)
    tp_u = cv2.undistortPoints(tp, cam_m, distor)
    return tp_u

def createProtocol(protocol_path, ros, cam_m, distor, in_dt):
    prot_fd = open(protocol_path, 'w')


    total_err = 0
    for ro in ros:
        tie_point_ids, tie_pts = in_dt.GetTiePoints([ro.ph1_id, ro.ph2_id])
        err = RelOrProtocol(tie_pts, ro.sps_pts, ro.pm, cam_m, distor)

        err_avrg =  sum(err) / len(err)
        prot_fd.write("%5d, %5d, %15.8f\n" % (ro.ph1_id, ro.ph2_id, err_avrg))

        total_err +=  err_avrg

    prot_fd.write("Total error: %15.8f\n" % (total_err/len(ros)))
    prot_fd.close()


def RelOrProtocol(tie_pts, pts_sps, P, cam_m, distor):
    
    tp1, tp2, tp1_u, tp2_u = UndistorTiePoints(tie_pts, cam_m, distor)

    err1 = computeError(tp2, pts_sps, P[:,:3], P[:,3], cam_m, distor)

    pm1 = np.eye(3, 4)
    err2 = computeError(tp1, pts_sps, pm1[:,0:3], pm1[:,3], cam_m, distor)
    err = err1 + err2


    return err


def UndistorTiePoints(tie_pts, cam_m, distor):

    tp1_u = undistortPoints(tie_pts[:, :2], cam_m, distor)
    tp2_u = undistortPoints(tie_pts[:, 2:4], cam_m, distor)

    shape = (1, len(tie_pts), 2)
    tp1 = np.empty(shape)
    tp2 = np.empty(shape)
    tp1[0,:,0] = cam_m[0, 0] * tp1_u[0,:,0] + cam_m[0, 2]
    tp1[0,:,1] = cam_m[1, 1] * tp1_u[0,:,1] + cam_m[1, 2]

    tp2[0,:,0] = cam_m[0, 0] * tp2_u[0,:,0] + cam_m[0, 2]
    tp2[0,:,1] = cam_m[1, 1] * tp2_u[0,:,1] + cam_m[1, 2]

    return tp1, tp2, tp1_u, tp2_u

def computeRelativeOrientation(tie_pts, cam_m, distor):

    tp1, tp2, tp1_u, tp2_u  = UndistorTiePoints(tie_pts, cam_m, distor)


    #idxs = np.nonzero(tie_pts[:,4])

    #F, mask = cv2.findFundamentalMat(tp1[0, :], tp2[0, :], param1=0.1, param2=0.95, method = cv2.FM_RANSAC)
    #em = cam_m.T.dot(F).dot(cam_m)
    

    em, mask = cv2.findEssentialMat(tp1[0, :], tp2[0, :], 
                                    threshold=0.001, 
                                    prob=0.95, 
                                    focal=cam_m[0, 0], 
                                    pp=(cam_m[0, 2], cam_m[1, 2]))

    F = np.linalg.inv(cam_m.T).dot(em).dot(np.linalg.inv(cam_m))

    p1, p2 = cv2.correctMatches(em, tp1_u, tp2_u)
    #p1 = tp1_u
    #p2 = tp2_u
    """
    print "chyby"
    GetPtsErrors(p1, tp1_u)
    GetPtsErrors(p2, tp2_u)
    """

    pts, R, t, mask = cv2.recoverPose(em, p1, p2)
    P = np.hstack((R, t))
    pm1 = np.eye(3, 4)

    pts_sps = compute3Dpoints(pm1, P, p1, p2)


    return pts_sps, P, p1, p2

    """
    #coords_f = '/home/ostepok/Desktop/coords.txt'
    #pts_fd = open(coords_f, mode='w')

        #pts_fd.write("%f|%f|%f\n" % tuple(pt))
    #pts_fd.close()

    p = grass.pipe_command("v.in.ascii", cat="0",
                                        input=coords_f,
                                        quiet=True,
                                        overwrite=True,
                                        output="rel_or_%i" % i_em,
                                        columns='x double precision, y double precision, z double precision',
                                        flags="z",
                                        z="3")
    p.wait()
    print p.stdout
    """
    #OutputArray _t, double focal, Point2d pp, InputOutputArray _mask)
    #for pt in tie_pts:
        #print "points"
        #pt_h = TriangulatePoint(pt[:2], pt[2:4], pm1, pms2[1])
        #print pt_h
        #pt_h = cv2.triangulatePoints(pm1, pms2[1], pt[:2], pt[2:4])
        #print pt_h
    #print PointDepth(pms2[1], pt_h)

def xp(f, X0, Y0, Z0 ,X,Y,Z,a):
    return f * (a[0,0]*(X-X0) + a[0,1]*(Y-Y0) + a[0,2]*(Z-Z0) ) / (a[2,0]*(X-X0) + a[2,1]*(Y-Y0) + a[2,2]*(Z-Z0) ) 

def yp(f, X0, Y0, Z0, X,Y,Z ,a):
    return f * (a[1,0]*(X-X0) + a[1,1]*(Y-Y0) + a[1,2]*(Z-Z0) ) / (a[2,0]*(X-X0) + a[2,1]*(Y-Y0) + a[2,2]*(Z-Z0) )

def _create5ptInFile(tie_pts, pts5_fd):

    i = 0
    for i_pt in tie_pts:
        print "%f %f %f %f\n" % tuple(i_pt)        
        i += 1
        if i == 5:
            break

        pts5_fd.write("%f %f %f %f\n" % tuple(i_pt))

    pts5_fd.flush()

def computeError(phs_pts, sps_pts, R, t, cam_m, distor):

    rvec, J = cv2.Rodrigues(R)
    #cam_m = np.eye(3)

    phs_pts_rep, J = cv2.projectPoints(sps_pts, rvec, t, cam_m, distor)
    sh = phs_pts_rep.shape
    phs_pts_rep = phs_pts_rep.reshape((sh[1], sh[0], sh[2]))

    s = phs_pts[0,...] - phs_pts_rep[0,...]
    errs = []
    for i in s:
        errs.append(np.linalg.norm(i))

    return np.array(errs)

def _create5ptOutFile(em_fd):

    ems = []
    i_line = 0;
    em_row = 0;

    em = np.zeros((3, 3))
    while 1:
        line = em_fd.readline()
        if not line:
            break

        l = line.split()
        if not l:
            em_row = 0;
            ems.append(em)
            em = np.zeros((3, 3))
            continue

        em[em_row, :3] = map(float, l)
        em_row += 1
    return ems

def GetPtsErrors(corr_pts, pts):
    errs = 1.0 - abs(corr_pts / pts)
    errs *= 100
    print "average"
    print np.average(errs, axis=1)
    print "max"
    print np.amax(errs, axis=1)
    print np.argmax(errs, axis=1)


def compute3Dpoints(P1, P2, npts1, npts2):

    ptsh3d = cv2.triangulatePoints(P1, P2, npts1.T, npts2.T).T

    pts_sps = deepcopy(ptsh3d[:,:3])
    for i, pt in enumerate(ptsh3d):
        pt = (pt / pt[3])[:3]
        pts_sps[i, :] = pt[:3]

    return pts_sps

def CameraMotionFromEssentialMatrix(em):
    U, W, V = np.linalg.svd(em, full_matrices=True)
    W = np.zeros((3, 3))
    W[0,1] = -1.0
    W[1,0] = 1.0
    W[2,2] = 1.0

    R1 = U.dot(W).dot(V.T)
    R2 =  U.dot(W.T).dot(V.T)

    t1 = np.expand_dims(U[2, :], axis=0)
    t2 = -t1

    pms = []
    for r in [R1, R2]:
        for t in [t1, t1]:
            pms.append(np.concatenate((r, t.T), axis=1))

    return pms

def TriangulatePoint(pt1, pt2, P1, P2):

    A = np.zeros((4, 4))

    A[0,:] = pt1[0] * P1[2,:] - P1[0,:]
    A[1,:] = pt1[1] * P1[2,:] - P1[1,:]
    
    A[2, :] = pt2[0] * P2[2,:] - P2[0,:]
    A[3, :] = pt2[1] * P2[2,:] - P2[1,:]

    pt_3d = np.zeros(4);
    U, W, V = np.linalg.svd(A, full_matrices=True)
    VT = V.T

    pt_3d[0] = VT[3,0];
    pt_3d[1] = VT[3,1];
    pt_3d[2] = VT[3,2];
    pt_3d[3] = VT[3,3];

    return pt_3d;

def PointDepth(P, pth):

    pth_phs = P.dot(pth)
    
    det = np.linalg.det(P[:,:3])
    if det > 0:
        sign = 1.0
    else:
        sign = -1.0

    w = pth_phs[2]
    W = P[0, 2];

    a = P[2, 0]
    b = P[2, 1]
    c = P[2, 2]

    m3_l = sqrt(a*a + b*b + c*c)

    return (w/W)*(sign/m3_l);

#sdef ExportCameras()
    """
    cmd = ["./5pt/5point_upravy/bin/Release/5point", pts5_fd.name,  em_fd.name]
    p = Popen(cmd, stdout=PIPE)
    p.wait()
    out, err = p.communicate()
    print out
    print err

    print "ok"
    ems = _create5ptOutFile(em_fd);
    """
