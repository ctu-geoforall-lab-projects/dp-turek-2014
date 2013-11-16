import os
import sys
import numpy as np
import math
from bba_derivs import x_derivs, y_derivs
from sympy import *

import scipy.misc as misc
import scipy.sparse as sparse
from  scipy.sparse.linalg import splu

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

    #eo_r[1, :] = eo_r[1, :] * -1
    a = np.array([[1,0,0],[0,1,0],[0,0,-1]])
    angles = GetBLUHAngles(a.dot(eo_r))
    for i, ang in enumerate(['ph', 'om', 'ka']):
        params[ang] = angles[i]

    for i, ax in enumerate(['xo_c', 'yo_c', 'zo_c']):
        params[ax] = o_c[i]

    o_pt = pt.GetGcpCoords()
    if o_pt is not None:
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


def GetL(in_dt, skip_phs, skip_pts):

    phs = in_dt.GetPhotos()
    
    params = {}
    
    L = []
    L_XO = []

    photos_n = 0
    image_points_n = 0
    rows_order = []
    pts_idxs = {}
    idxs_pts = []

    tie_pts_n = 0

    use_photos = phs.keys()
    #use_photos = [598, 597]
    idxs_phs = []

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

            v = GetParamsVals(ph, pt)
            #if not i use_photos: 
            pt_phs = pt.GetPhotoPoints().keys()
            if len((set(pt_phs) & set(use_photos))) < 2:
                continue   

            if v is not None:
                if not params.has_key(ph):
                    params[ph] = {}
                    photos_n += 1
                    idxs_phs.append(ph) 

                if not pts_idxs.has_key(pt) and not v['gcp']:
                    pts_idxs[pt] = tie_pts_n
                    idxs_pts.append(pt) 
                    tie_pts_n +=  1

                params[ph][pt] = v
                rows_order.append((ph, pt, v))

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

    return L, L_XO, rows_order, pts_idxs, idxs_pts, idxs_phs, tie_pts_n, image_points_n, photos_n

def RunBundleBlockAdjutment(in_dt,  skip_phs, skip_pts):

    unknowns = ['ph', 'om', 'ka', 'xo_c', 'yo_c', 'zo_c', 'xo_pt', 'yo_pt', 'zo_pt', 'yi_c', 'xi_c', 'f']
    """
    xc, yc = CreateSymbolicColEquations()

    derivs = {}
    for un in unknowns:
        derivs[un] = (diff(xc, un, 1), diff(yc, un, 1))
        print un
        print diff(xc, un, 1)
        print diff(yc, un, 1)

    return
    """

    L, L_XO, rows_order, pts_idxs, idxs_pts, idxs_phs, tie_pts_n, image_points_n, photos_n = GetL(in_dt, skip_phs, skip_pts)

    
    L_XO = np.array(L_XO)
    L = np.array(L)


    ph_unknowns = ['ph', 'om', 'ka', 'xo_c', 'yo_c', 'zo_c']
    pt_unknowns = ['xo_pt', 'yo_pt', 'zo_pt']

    ph_step = len(ph_unknowns)
    pt_step = len(pt_unknowns)

    params_n = photos_n * ph_step + tie_pts_n * pt_step
    A = np.zeros((len(L), params_n))

    i_ph = 0
    i_phpt = 0

    pt_zero_idx = photos_n * ph_step

    prev_ph = None

    X0 = np.zeros((params_n))
    for i, d in enumerate(rows_order):
        ph, pt, v = d
        if prev_ph is None:
            i_ph_idx = 0
            prev_ph = ph
        elif prev_ph != ph:
            prev_ph = ph
            i_ph_idx += 1

        row_x = i * 2
        row_y = i * 2 + 1

        for i_un, un in enumerate(ph_unknowns):
            c = i_ph_idx * ph_step + i_un
            A[row_x, c] = x_derivs[un](**v)
            A[row_y, c] = y_derivs[un](**v)
            X0[c] = v[un]

        if not v['gcp']:
            pt_idx = pt_zero_idx + pts_idxs[pt] * pt_step
            for i_un, un in enumerate(pt_unknowns):
                c = pt_idx + i_un
                A[row_x, c] = x_derivs[un](**v)
                A[row_y, c] = y_derivs[un](**v)
                X0[c] = v[un]

        i_phpt = 0
        i_ph = 0



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
    Ni = np.linalg.inv(N)
    
    dx = Ni.dot(n)

    Xa = X0 + dx
    e = A.dot(dx) - l

    for pt_idx in range(tie_pts_n):
        pt = idxs_pts[pt_idx]
        pt_x_idx = pt_zero_idx + pt_idx * pt_step

        pt_coords = np.zeros((3))
        pt_coords_dx = np.zeros((3))

        for i in range(pt_step):
            pt_coords[i] = Xa[pt_x_idx + i]
            pt_coords_dx[i] = dx[pt_x_idx + i]

        pt.SetCoords(pt_coords)

    for ph_idx in range(photos_n):
        ph = idxs_phs[ph_idx]
        ph_x_idx = ph_idx * ph_step

        ph_eo = np.zeros((ph_step))
        ph_eo_dx = np.zeros((ph_step))
        for i in range(ph_step):
            ph_eo[i] = Xa[ph_x_idx + i]
            ph_eo_dx[i] = dx[ph_x_idx + i]

        a = np.array([[1,0,0],[0,1,0],[0,0,-1]])

        eo_r = a.T.dot(RotMatBLUH(*ph_eo[:3]))

        eo = np.hstack((eo_r, np.array([ph_eo[3:]]).T))

        old_eo = ph.GetEO()

        ph.SetEO(eo)


    L1, L_XO1, rows_order, pts_idxs, idxs_pts, idxs_phs, tie_pts_n, image_points_n, photos_n = GetL(in_dt,  skip_phs, skip_pts)
    e1 = L_XO1 - L

    print "e"
    print e - e1

    #for in photos_n * ph_step:

    """
    xx = np.where(A != 0.0)
    A[xx] = 200
    misc.imsave('/home/ostepok/Desktop/A.png', A)
    """
    print independent_columns(A, tol = 1e-05)
    np.savetxt('/home/ostepok/Desktop/A.csv', A, fmt='%.18e', delimiter=',')


    np.savetxt('/home/ostepok/Desktop/res.csv', A, fmt='%.3e', delimiter=',')
    return


def independent_columns(A, tol = 1e-05):
    """
    Return an array composed of independent columns of A.

    Note the answer may not be unique; this function returns one of many
    possible answers.

    http://stackoverflow.com/q/13312498/190597 (user1812712)
    http://math.stackexchange.com/a/199132/1140 (Gerry Myerson)
    http://mail.scipy.org/pipermail/numpy-discussion/2008-November/038705.html
        (Anne Archibald)
    http://stackoverflow.com/questions/13312498/how-to-find-degenerate-rows-columns-in-a-covariance-matrix
    >>> A = np.array([(2,4,1,3),(-1,-2,1,0),(0,0,2,2),(3,6,2,5)])
    >>> independent_columns(A)
    np.array([[1, 4],
              [2, 5],
              [3, 6]])
    """
    Q, R = np.linalg.qr(A)
    independent = np.where(np.abs(R.diagonal()) < tol)[0]
    return independent


