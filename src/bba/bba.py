import os
import sys
import numpy as np
import math
from bba_derivs import x_derivs, y_derivs
from sympy import *

import scipy.misc as misc
import scipy.sparse as sparse
from  scipy.sparse.linalg import splu

from collections import namedtuple
from sympy.matrices import *

def RotMatBLUH(ph, om, ka):
    rotm = np.zeros((3,3))
    _rotMatBLUH(ph, om, ka, rotm)

    return rotm

def GetBLUHAngles(rm):

    ka = math.atan2(rm[0, 1], rm[1, 1])
    #print ka * 180 / pi
    ph = math.atan2(rm[2, 0], rm[2, 2])
    om = math.atan2(-rm[2, 1], sqrt(rm[1, 1] * rm[1, 1] + rm[0, 1] * rm[0, 1]))

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

def xp_col_eq(xi_c, f, xo_c, yo_c, zo_c ,xo_pt, yo_pt, zo_pt, r, **kwargs):
    return xi_c - f * (r[0][0]*(xo_pt-xo_c) + r[0][1]*(yo_pt-yo_c) + r[0][2]*(zo_pt-zo_c)) / \
                      (r[2][0]*(xo_pt-xo_c) + r[2][1]*(yo_pt-yo_c) + r[2][2]*(zo_pt-zo_c)) 

def yp_col_eq(yi_c, f, xo_c, yo_c, zo_c, xo_pt, yo_pt, zo_pt, r, **kwargs):
    return yi_c - f * (r[1][0]*(xo_pt-xo_c) + r[1][1]*(yo_pt-yo_c) + r[1][2]*(zo_pt-zo_c)) / \
                      (r[2][0]*(xo_pt-xo_c) + r[2][1]*(yo_pt-yo_c) + r[2][2]*(zo_pt-zo_c))

def GetParamsVals(ph, pt):

    params = {}

    cam = ph.GetCamera()

    cam_m, distor = cam.GetParams()

    params['f'] = cam_m[0, 0]
    params['xi_c'] = cam_m[0, 2]
    params['yi_c'] = cam_m[1, 2]

    ph_eo = ph.GetEO()

    o_c = ph_eo[:,3]
    eo_r = ph_eo[:,:3]

    # reflection caused by different handness of Opencv and BLUH systems 
    a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
    angles = GetBLUHAngles(a.dot(eo_r))
    for i, ang in enumerate(['ph', 'om', 'ka']):
        params[ang] = angles[i]

    for i, ax in enumerate(['xo_c', 'yo_c', 'zo_c']):
        params[ax] = o_c[i]

    o_pt, control = pt.GetGcp()
    if o_pt is not None and not control:
        params["gcp"] = True 
    else:
        params["gcp"] = False
        o_pt = pt.GetCoords()

    if o_pt is None:
        return None

    for i, ax in enumerate(['xo_pt', 'yo_pt', 'zo_pt']):
        params[ax] = o_pt[i]

    phpts = pt.GetPhotoPoints()

    if len(phpts) < 2:
        return None

    i_pt = phpts[ph]
    params['xi_pt'] = i_pt[0]
    params['yi_pt'] = i_pt[1]

    return params
    #.evalf(<args>).

BbaIndexes = namedtuple('BbaIndexes', 'Lrows_idxs,' 
                                      'pts_idxs,  idxs_pts, '
                                      'phs_idxs, idxs_phs, '
                                      'cams_idxs, idxs_cams, '
                                      'cams_n, photos_n, tie_pts_n, '
                                      'image_points_n')

def GetL(in_dt, skip_phs, skip_pts):

    cams = in_dt.GetCameras()
        
    L = []
    L_XO = []


    phs = in_dt.GetPhotos()

    #TODO for testing
    use_photos = phs.keys()
    #use_photos = [598, 597]

    cams_idxs = {}
    idxs_cams = []

    phs_idxs = {}
    idxs_phs = []
    
    pts_idxs = {}
    idxs_pts = []
    
    Lrows_idxs = []

    photos_n = 0
    image_points_n = 0
    tie_pts_n = 0
    cams_n = 0

    for cam in cams.itervalues():
        phs = cam.GetPhotos()

        for ph in phs.itervalues():
            if ph in skip_phs:
                continue

            ph_pts = ph.GetPoints()

            photo_points_n = 0
            if ph.GetId() not in use_photos: 
                continue

            #print "PHOTO"
            #print ph.GetId()

            for pt in ph_pts.itervalues():
                if pt in skip_pts:
                    continue

                #if not i use_photos: 
                pt_phs = pt.GetPhotoPoints().keys()
                if len((set(pt_phs) & set(use_photos))) < 2:
                    continue   

                v = GetParamsVals(ph, pt)
                if v is not None:
                    if cam not in idxs_cams:
                        idxs_cams.append(cam) 
                        cams_idxs[cam] = cams_n
                        cams_n += 1

                    if ph not in idxs_phs:
                        idxs_phs.append(ph) 
                        phs_idxs[ph] = photos_n
                        photos_n += 1

                    if pt not in pts_idxs and not v['gcp']:
                        pts_idxs[pt] = tie_pts_n
                        idxs_pts.append(pt) 
                        tie_pts_n +=  1

                    Lrows_idxs.append((cam, ph, pt, v))

                    r = RotMatBLUH(v['ph'], v['om'], v['ka'])

                    L_XO.append(xp_col_eq(r=r, **v))
                    L_XO.append(yp_col_eq(r=r, **v))

                    L.append(v['xi_pt'])
                    L.append(v['yi_pt'])

                    photo_points_n += 1
                else:
                    pass
                    #print "skip %d" % pt.GetId()

        image_points_n += photo_points_n
        bba_idxs = BbaIndexes(Lrows_idxs, 
                              pts_idxs, idxs_pts, phs_idxs, idxs_phs, cams_idxs, idxs_cams, 
                              cams_n, photos_n, tie_pts_n, image_points_n)
    return L, L_XO, bba_idxs

def AdjustData(Xa, bba_idxs, unk, apo):
    bi = bba_idxs

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    for cam_idx in range(bi.cams_n):
        cam = bi.idxs_cams[cam_idx]
        cam_x_col = apo.cam_0col + cam_idx * cam_step

        cam_io = {}
        for i, un in enumerate(unk.cam):
            cam_io[un] = Xa[cam_x_col + i]

        cam.SetIO(cam_io)

    for pt_idx in range(bi.tie_pts_n):

        pt = bi.idxs_pts[pt_idx]
        pt_x_col = apo.pt_0col + pt_idx * pt_step

        pt_coords = np.zeros((3))

        for i in range(pt_step):
            pt_coords[i] = Xa[pt_x_col + i]
        pt.SetCoords(pt_coords)

    for ph_idx in range(bi.photos_n):
        ph = bi.idxs_phs[ph_idx]
        ph_x_col = apo.ph_0col + ph_idx * ph_step

        ph_eo = np.zeros((ph_step))
        for i in range(ph_step):
            ph_eo[i] = Xa[ph_x_col + i]

        # reflection back TODO error prone, merge it
        a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
        eo_r = a.T.dot(RotMatBLUH(*ph_eo[:3]))
        eo = np.hstack((eo_r, np.array([ph_eo[3:6]]).T))
        ph.SetEO(eo)


Unknowns = namedtuple('Unknowns', 'cam, ph, pt')

def CreateIterationProtocol(i_num, cov, Xa, dx, e, bi, apo, prot_fd, unk, gcps):

    prot_fd.write(_('\n\n\n\n*********************************************************************'))
    prot_fd.write(_('\n\nIteration: %d\n\n' % (i_num + 1)))

    prot_fd.write(_('Ground control point differences:\n'))

    sigmas = np.sqrt(np.diagonal(cov))
    prot_fd.write("%10s, %10s\n" % ("gcp id", "diff"))
    for gcp in gcps.itervalues():
        if gcp.GetCoords() is None or not gcp.GetGcp()[1]:
            continue

        dist = np.linalg.norm(gcp.GetGcp()[0] - gcp.GetCoords())
        bingo_dist = np.linalg.norm(gcp.GetGcp()[0] - gcp.GetResultCoords())

        prot_fd.write("%10d, %10.4f, %10.4f" % (gcp.GetId(), dist, bingo_dist))
        prot_fd.write("\n")

    prot_fd.write("\n\n")

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    def _printParamsLine(x, unk_l, x_col):

        for i, un in enumerate(unk_l):
            val = x[x_col + i]
            if un in ['ph', 'om', 'ka']:
                val = val * 180 / pi

            prot_fd.write("%15.4f" % (val))

    def _printParamCaptions(unk_l):
            for i, un in enumerate(unk_l):
                prot_fd.write("%15s" % (un))
            prot_fd.write("\n\n")

    def _printFeature(Xa, dx, sigmas, unk_l, x_col):
        _printParamsLine(Xa, unk_l, x_col)
        prot_fd.write("\n")
        _printParamsLine(dx, unk_l, x_col)
        prot_fd.write("\n")
        _printParamsLine(sigmas, unk_l, x_col)
        prot_fd.write("\n\n")

    prot_fd.write(_('Interior orientaitons adjustment results\n'))

    for cam_idx in range(bi.cams_n):
        if cam_idx == 0:
            _printParamCaptions(unk.cam)

        cam = bi.idxs_phs[cam_idx]
        prot_fd.write(_("Camera %d\n" % cam.GetId()))
        cam_x_col = apo.cam_0col + cam_idx * cam_step

        _printFeature(Xa, dx, sigmas, unk.cam, cam_x_col)

    prot_fd.write(_('Exterior orientaitons adjustment results\n'))

    for ph_idx in range(bi.photos_n):
        if ph_idx == 0:
            _printParamCaptions(unk.ph)

        ph = bi.idxs_phs[ph_idx]
        prot_fd.write(_("Photo %d\n" % ph.GetId()))
        ph_x_col = apo.ph_0col + ph_idx * ph_step

        _printFeature(Xa, dx, sigmas, unk.ph, ph_x_col)

    prot_fd.write(_('Object points orientaitons adjustment results\n'))

    for pt_idx in range(bi.tie_pts_n):
        if pt_idx == 0:
            _printParamCaptions(unk.pt)

        pt = bi.idxs_pts[pt_idx]
        prot_fd.write(_("Point %d\n" % pt.GetId()))
        pt_x_col = apo.pt_0col + pt_idx * pt_step

        _printFeature(Xa, dx, sigmas, unk.pt, pt_x_col)

def RunBundleBlockAdjutment(in_dt,  skip_phs, skip_pts, protocol_path):

    prot_fd = open(protocol_path, 'w')

    cam_unk = []
    
    ph_unk = ['ph', 'om', 'ka', 'xo_c', 'yo_c', 'zo_c']
    pt_unk = ['xo_pt', 'yo_pt', 'zo_pt']
    unk =  Unknowns(cam_unk, ph_unk, pt_unk)

    for i in range(15):
    
        if i == 8:
            cam_unk = ['f', 'xi_c', 'yi_c']
            unk =  Unknowns(cam_unk, ph_unk, pt_unk)
        """
        if i == 8:
            ph_unk = []
            pt_unk = []
            cam_unk = ['f', 'xi_c', 'yi_c']
            unk =  Unknowns(cam_unk, ph_unk, pt_unk)
        """
        cov, Xa, dx, e, bi, apo = RunBBAIteration(in_dt, skip_phs, skip_pts, unk)
        CreateIterationProtocol(i, cov, Xa, dx, e, bi, apo, prot_fd, unk, in_dt.GetGcps())

    prot_fd.close()

AdjustedParametersOrder = namedtuple('AdjustedParameters', 'cam_0col, ph_0col, pt_0col, cols_n')
def GetAdjParamsOrder(unk, bi):
    cam_zero_col = 0
    ph_zero_col = bi.cams_n * len(unk.cam)
    pt_zero_col = bi.photos_n * len(unk.ph) + ph_zero_col

    col_n = ph_zero_col + pt_zero_col + bi.tie_pts_n * len(unk.pt)

    return cam_zero_col, ph_zero_col, pt_zero_col, col_n

def RunBBAIteration(in_dt,  skip_phs, skip_pts, unk):

    L, L_XO, bi = GetL(in_dt, skip_phs, skip_pts)
    
    L_XO = np.array(L_XO)
    L = np.array(L)

    cam_step = len(unk.cam)
    ph_step = len(unk.ph)
    pt_step = len(unk.pt)

    apo = AdjustedParametersOrder(*GetAdjParamsOrder(unk, bi))

    A = np.zeros((len(L), apo.cols_n))
    X0 = np.zeros((apo.cols_n))

    prev_cam = None
    prev_ph = None
    for i, d in enumerate(bi.Lrows_idxs):
        cam, ph, pt, v = d
        if prev_ph is None:
            i_ph_idx = 0
            prev_ph = ph
        elif prev_ph != ph:
            prev_ph = ph
            i_ph_idx += 1

        if prev_cam is None:
            i_cam_idx = 0
            prev_cam = cam
        elif prev_ph != ph:
            prev_cam = cam
            i_cam_idx += 1

        row_x = i * 2
        row_y = i * 2 + 1


        cam_idx = apo.cam_0col + bi.cams_idxs[cam] * cam_step
        for i_un, un in enumerate(unk.cam):
            c = cam_idx + i_un
            A[row_x, c] = x_derivs[un](**v)
            A[row_y, c] = y_derivs[un](**v)
            X0[c] = v[un]

        ph_idx = apo.ph_0col + bi.phs_idxs[ph] * ph_step
        for i_un, un in enumerate(unk.ph):
            c = ph_idx + i_un
            A[row_x, c] = x_derivs[un](**v)
            A[row_y, c] = y_derivs[un](**v)
            X0[c] = v[un]

        if not v['gcp']:
            pt_idx = apo.pt_0col + bi.pts_idxs[pt] * pt_step
            for i_un, un in enumerate(unk.pt):
                c = pt_idx + i_un
                A[row_x, c] = x_derivs[un](**v)
                A[row_y, c] = y_derivs[un](**v)
                X0[c] = v[un]

    #print np.linalg.cholesky(N)
    
    #n_sps = sparse.csc_matrix(N)
    #lu_obj = splu(n_sps)
    #Ni = lu_obj.solve(np.eye(N.shape[0])).T
    
    """
    print N.dot(Ni)
    Ns = Matrix(N)
    Ni = Ns.inv()
    print Ni.dot(N)
    """

    l = L - L_XO
    P = np.eye(A.shape[0])

    N = A.T.dot(P).dot(A)
    n = A.T.dot(P).dot(l)
    #Ni = np.linalg.inv(N)
    Ni = np.linalg.pinv(N)
    dx = Ni.dot(n)

    Xa = X0 + dx
    e = A.dot(dx) - l

    AdjustData(Xa, bi, unk, apo)

    L1, L_XO1, bi = GetL(in_dt,  skip_phs, skip_pts)
    e1 = L_XO1 - L

    print "e"
    print e - e1

    #for in photos_n * ph_step:

    """
    xx = np.where(A != 0.0)
    A[xx] = 200
    misc.imsave('/home/ostepok/Desktop/A.png', A)
    np.savetxt('/home/ostepok/Desktop/A.csv', A, fmt='%.18e', delimiter=',')
    np.savetxt('/home/ostepok/Desktop/res.csv', A, fmt='%.3e', delimiter=',')
    """
    r = A.shape[0] - A.shape[1]
    sigma0_2 = e.T.dot(P).dot(e) / r
    cov =  sigma0_2 * Ni

    return cov, Xa, dx, e, bi, apo

