import os
from math import pi
from collections import namedtuple
from collections import Counter
import numpy as np
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

    def GetPhotoPoints(self):
        return self.ph_pts

    def AddPhotoPoint(self, ph, coords):
        self.ph_pts[ph] = coords
        ph.AddPoint(self)

    def SetGcp(self, coords):
        self.gcp = coords

    def GetGcpCoords(self):
        return self.gcp

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

        self.res_c = None
        self.res_rot = None

        self.eo = None
        self.c = None

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

    def SetResultExtOr(self, coords = None, rot = None):
        self.res_c = coords
        self.res_rot = rot

    def GetResultExtOr(self):
        return self.res_c, self.res_rot

    def GetNeighboroughPhotos(self):

        phs = [ph for pt in self.ph_pts.itervalues() 
                  for ph in pt.GetPhotoPoints().iterkeys() 
                  if ph != self]

        phs_count = Counter(phs)
        return phs_count

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

    def GetParams(self):
        return self.cam_m, self.distor

    def AddPhoto(self, ph):
        self.phs[ph] = ph

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
        self.cams, self.phs, self.pts = self.parser.GetData()

        self.phs[577].GetNeighboroughPhotos()
        
    def GetPhotosConnectivity(self):

        neigh_phs = dict( ((ph, ph2), num) 
                            for ph in sorted(self.phs.values())
                            for ph2, num in ph.GetNeighboroughPhotos().iteritems()
                            if ph2 > ph)
            

        return sorted(neigh_phs, key=neigh_phs.get, reverse=True)

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

class BingoParser:
    def __init__(self, path):   
        self.path = path

    def Parse(self):

        self.cams, self.pts = self._readCameraFile()
        self.phs = self._readPhotosFile()
        self._readPhotoPointsFile(self.pts, self.phs)

        self.res_eo = self._readResult(self.phs, self.pts)

    def GetData(self):
        return self.cams, self.phs, self.pts

    def GetResult(self):
        return self.res_eo

    def _readResult(self, phs, pts):
        f = open(self.path + '/old-itera.dat')

        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "*" or not line.strip():
                continue

            lin=line.rstrip("")
            l = lin.split()
            if not l:
                continue
            if "ORIA" in l[0]:

                ph_id = int(l[1])

                eo = map(float, l[2:-1])
                eo[3:] = map(lambda d : d / 200 * pi, eo[3:])

                phs[ph_id].SetResultExtOr(coords=eo[:3], rot=eo[3:])

            elif "CORD" in l[0]:
                pt_id = int(l[1])
                c = map(float, l[2:5])
                pts[pt_id].SetResultCoords(c)

        return phs

    def _readCameraFile(self):
        f = open(self.path + '/geoin.dat')

        cams = {}
        pts = {}

        while 1:
            line = f.readline()
            if not line:
                break
           
            lin=line.strip()
            
            l = lin.split()
            if not l:
                continue
           
            if l[0] == '*' or l[0] == 'C':
                continue
            #TODO implemented for just one camera 
            if 'CAPA' in l[0]:
                fl = float(l[3])
                x0 = float(l[4])
                y0 = float(l[5])

            if 'ECCA' in l[0]:
                lge = float(l[2])
                lgn = float(l[3])
                lgh = float(l[4])

                distor = np.array([0, 0, 0, 0], dtype=float)
                cam_m = np.array([[fl, 0, x0],
                                  [0, fl, y0],
                                  [0, 0, 1]])

                self.cam = cams[l[1]] = Camera(l[0], cam_m, distor)

            if 'CONT' in l[0] or 'CHCK' in l[0]:
                gcp_id = int(l[1])
                e = float(l[2])
                n = float(l[3])
                h = float(l[4])

                c = np.array([e, n, h])

                pts[gcp_id] = Point(gcp_id)
                pts[gcp_id].SetGcp(c)
        
        return cams, pts

    def _readPhotosFile(self):
        f = open(self.path + '/gps.dat')

        phs = {}
        count = 0

        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue

            lin=line.rstrip("\r\n")
            l = lin.split()
            if not l:
                continue
            if "LINE" in l[0]:
                continue

            ph_id = int(l[0])
            #TOOD just one cam
            cam = self.cam

            eo = map(float, l)
            phs[ph_id] = Photo(ph_id, self.cam) 
                               #np.array(eo), TODO
                               #None
                     
        return phs

    def _readPhotoPointsFile(self, pts, phs):
        f = open(self.path + '/image.dat')

        count = 0

        new_ph = True

        #TODO works with just one camere
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue

            lin=line.strip()
            l = lin.split()
            if not l:
                continue

            if '-99' == l[0]:
                new_ph = True
                continue

            if new_ph:
                #if len(l) > 1:
                #TODO camera
                ph_id = int(l[0])
                new_ph = False
                continue

            try:
                pt_id = int(l[0])
            except:
                pt_id = int(l[0][1:])

            #TODO 
            c = np.array([float(l[1]), float(l[2])])

            if not pts.has_key(pt_id):
                pts[pt_id] = Point(pt_id)

            pts[pt_id].AddPhotoPoint(phs[ph_id], c)

class GRASSFilesParser:
    # Authors: 
    # original author Man Chao, Malaysia 
    # pythonized by Yann Chemin
    # further development:Stepan Turek

    def __init__(self, FN):   
        self.FN = FN

    def Parse(self):
        pass

    def Read_Photo_file(self):
        f = open(self.FN["PHOTOINIT"])

        m = {}
        count = 0
        while 1:
            line = f.readline()
            if not line:
                break

            if line[0] == "#":
                continue

            lin=line.rstrip()
            if count == 0:
                l=lin.split()

                att_setts = AtSettings(*map(int, l))
            else:
                l = lin.split("\t")
                p_id =  map(int, l[0:2])
                p_coords = map(float, l[2:])

                p = Photo(*(p_id + p_coords))
                m[p_id[0]] = p

            count=count+1

        photos = Photos(att_setts, m)
        return photos

    Gcp = namedtuple('Gcp', 'id, e, n, elev')
    Gcps = namedtuple('Gcps', 'melev, m')

    #Read GCP File
    def Read_GCP_file(self):
        f = open(self.FN["GCP"])
        
        m = {}
        count = 0    
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue        
            else:
                lin=line.rstrip("\r\n")
                l = lin.split()
                g_id = int(l[0])
                g_coords = map(float, l[1:])
                
                g = Gcp(g_id, np.array(g_coords))
                m [g_id] = g
                count=count+1

                melev = g.elev

        #TODO MElev!!!!
        gcps = Gcps(melev, m)
        return gcps
        
    #Read CAMERA File
    def Read_Cam_file(self):
        f = open(self.FN["CAMERA"])

        cams = {}
        count = 0
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                pass        
            else:
                lin=line.rstrip("\r\n")
                l = lin.split()

                cam_id = int(l[0])
                #TODO magic number (more general scaling)?
                c = Camera(cam_id, float(l[1]) / (-1000.0))
                cams[cam_id] = c
                count=count+1
        
        return cams

    #Read measure Points
    def ReadMeasurePoints(self, photo_id, adj_row):

        conPhoto = os.path.join(self.FN['GROUP'], str(photo_id), "CONTROL_POINTS")
        f = open(conPhoto)

        m_pts = []
        count = 0
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue
            else:
                lin=line.rstrip("\r\n")
                l = lin.split()
                if l[6] == "1":
                    l = map(float, l)
                    p = MeasurePt(l[5], photo_id, l[0], adj_row - l[1], l[2] / (1000.0))
                    m_pts.append(p)
                    count=count+1

        return m_pts


    """
    RefPoints = namedtuple('RefPoints', 'e, n, x, y, status')

    def ReadRefPointsFile(self, id):
        ref_pts_path = os.path.join(self.FN['GROUP'], str(id), "REF_POINTS")
        f = open(ref_pts_path)

        ref_points = []
        count = 0
        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                pass        
            else:
                lin=line.rstrip("\r\n")
                l = lin.split()
                #TODO magic number (more general scaling)?
                l = map(float, l)
                l[4] = int(l[4]) 

                r = RefPoints(*l)
                ref_points.append(r)
                count=count+1
        f.close()
        return ref_points

    def ReadCellHeader(self, id):
        ref_pts_path = os.path.join(self.FN['CELLHD'], str(id))
        f = open(ref_pts_path)

        ref_points = []
        i = 0
        while 1:
            line = f.readline()
            if not line:
                break
            if i == 7:
                rows = int(line.split()[1])
                break
            i=i+1

        f.close()

        return rows
    """