import os
from math import pi
from collections import Counter
import numpy as np

from bba import RotMatBLUH
#AtSettings = namedtuple('AtSettings', 'scale, focusl, iter, atiter, lm')

#TODO include standard deviations
class Point:
    def __init__(self, pt_id):
        self.ph_pts = {}
        self.gcp = None
        self.pt_id = pt_id
        self.ro_c = {}
        self.coords = None

        self.res_c = None
        self.control_gcp = False
    def GetPhotoPoints(self):
        return self.ph_pts

    def AddPhotoPoint(self, ph, coords):
        self.ph_pts[ph] = coords
        ph.AddPoint(self)

    def SetGcp(self, coords, control):
        self.gcp = coords
        self.control_gcp = control

    def GetGcp(self):
        return self.gcp, self.control_gcp

    def IsControlGcp():
        return self.control_gcp

    def AddRelOrCoords(self, rel_or_id, coords):
        self.ro_c[rel_or_id] = coords

    def GetRelOrsCoords(self):
        return self.ro_c

    def SetCoords(self, coords):
        self.coords = coords

    def GetCoords(self):
        return self.coords

    def GetId(self):
        return self.pt_id

    def SetResultCoords(self, coords = None):
        self.res_c = coords

    def GetResultCoords(self):
        return self.res_c

    def __hash__(self):
        return hash(self.GetId())

    def __eq__(self, other):
        if isinstance(other, Point):
            return self.GetId() == other.GetId()
        return self.GetId() == other

class Photo:
    def __init__(self, ph_id, cam):
        self.ph_id = ph_id 
        self.ph_pts = {}
        self.cam = cam
        cam.AddPhoto(self)

        self.res_eo = None

        self.eo = None
        self.c = None

        self.neigh_phs = None
    def HelmertTransformEO(self, ht_mat, ht_scale):
        self.eo[:,:3] = self.eo[:,:3].dot(ht_mat[:,:3].T)
        self.eo[:,3] = ht_mat[:,3] + ht_scale * ht_mat[:,:3].dot(self.eo[:,3])

    def SetEO(self, eo):
        self.eo = eo

    def SetP(self, P):
        if self.eo is None:
            self.eo = np.empty((3,4))
        self.eo[:,:3] = P[:,:3]
        self.eo[:,3] = - self.eo[:,:3].T.dot(P[:,3])

    def GetEO(self):
        return self.eo

    def GetP(self):
        if self.eo is None:
            return None


        P = np.empty((3,4))
        P[:,:3] = self.eo[:,:3]
        P[:,3] = - self.eo[:,:3].dot(self.eo[:,3])

        return P

    def SetPMat(self, eo):
        self.eo = eo

    def GetEO(self):
        return self.eo

    def SetCoords(self, c):
        self.c = c

    def GetCoords(self):
        return self.c

    def AddPoint(self, pt):
        self.ph_pts[pt] = pt

    def GetId(self):
        return self.ph_id

    def GetPoints(self):
        return self.ph_pts

    def SetResultExtOr(self, eo):
        self.res_eo =  eo

    def GetResultExtOr(self):
        return self.res_eo

    def GetCamera(self):
        return self.cam

    def GetNeighboroughPhotos(self):
        if self.neigh_phs is None:
            phs = [ph for pt in self.ph_pts.itervalues() 
                      for ph in pt.GetPhotoPoints().iterkeys() 
                      if ph != self]

            self.neigh_phs = Counter(phs)
    
        return self.neigh_phs

    def __hash__(self):
        return hash(self.GetId())

    def __eq__(self, other):
        if isinstance(other, Photo):
            return self.GetId() == other.GetId()

        return self.GetId() == other

class Camera:
    def __init__(self, cam_id, cam_m, distor):
        self.phs = {}
        self.cam_id = cam_id

        self.cam_m = cam_m
        self.distor = distor

    def GetChipSize(self):
        return self.chip_size

    def SetChipSize(self, chip_size):
        self.chip_size = chip_size

    def SetIO(self, io):

        if 'f' in io: 
            self.cam_m[0, 0] = io['f']
            self.cam_m[1, 1] = io['f']

        if 'xi_c' in io: 
            self.cam_m[0, 2] = io['xi_c']

        if 'yi_c' in io:
            self.cam_m[1, 2] = io['yi_c']

        if 'r1_dis' in io:
            self.distor[0] = io['r1_dis']

    def GetParams(self):
        return self.cam_m, self.distor

    def AddPhoto(self, ph):
        self.phs[ph] = ph

    def GetPhotos(self):
        return self.phs

    def GetId(self):
        return self.cam_id

    def __hash__(self):
        return hash(self.GetId())

    def __eq__(self, other):
        if isinstance(other, Camera):
            return self.GetId() == other.GetId()

        return self.GetId() == other
"""
def _GPtId(pt_id): 
    return ('pt', pt_id)

def _GPhId(ph_id):
    return ('ph', ph_id)

def _PtGId(gpt_id): 
    return gpt_id[1]

def _PhGId(gph_id):
    return gph_id[1]
"""

class InputData:
    def __init__(self, parser):
        self.parser = parser

        self.parser.Parse()
        self.cams, self.phs, self.pts, self.gcps = self.parser.GetData()

    def GetPhotosConnectivity(self):

        neigh_phs = dict( ((ph, ph2), num) 
                            for ph in sorted(self.phs.values())
                            for ph2, num in ph.GetNeighboroughPhotos().iteritems()
                            if ph2 > ph)
            

        return sorted(neigh_phs, key=neigh_phs.get, reverse=True), neigh_phs

    def GetCameras(self):
        return self.cams

    def GetTiePoints(self, ph_ids):

        tie_pts_phs = map(lambda ph_id : self.phs[ph_id], ph_ids)

        phs_pts = map(lambda ph : ph.GetPoints(), tie_pts_phs)

            
        phs_pts = map(lambda ph : set(ph), phs_pts)

        tie_pts = set.intersection(*phs_pts)
        tie_pt_ids = map(lambda pt : pt.GetId(), tie_pts)

        tie_pt_phpts = map(lambda pt : pt.GetPhotoPoints(), tie_pts)

        tie_pts_arr = np.array([phpts[ph_id] for phpts in tie_pt_phpts 
                                             for ph_id in ph_ids])

        sh = tie_pts_arr.shape
        tie_pts_arr = tie_pts_arr.reshape(sh[0] / len(ph_ids), ph_ids)

        return tie_pt_ids, tie_pts_arr

    def GetResults(self, pt_ids):
        pts = [self.pts[pt_id] for pt_id in pt_ids]

        return np.array([pt.GetResultCoords() for pt in pts])

    def GetPhotos(self):
        return self.phs

    def GetPoints(self):
        return self.pts

    def GetGcps(self):
        return self.gcps
