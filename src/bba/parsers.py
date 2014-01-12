import os
from math import pi
import numpy as np

from bba import RotMatBLUH
from readpar import Point, Camera, Photo

from relative_or_zal import undistortPoints

def undistortPoint(pt, cam_m, distor):

    x_c = cam_m[0, 1]
    y_c = cam_m[0, 2]

    x = (pt[0] - x_c)    
    y = (pt[1] - y_c)

    k1 = distor[0]
    k2 = distor[1]

    p1 = distor[2]
    p2 = distor[3]
    
    r_sq = x ** 2 + y ** 2

    x_u = x * (1 + k1 * r_sq + k2 * r_sq ** 2) + p1 * (r_sq + 2 * x ** 2) + 2 * p2 * x * y
    y_u = y * (1 + k1 * r_sq + k2 * r_sq ** 2) + p2 * (r_sq + 2 * y ** 2) + 2 * p1 * x * y
    return x_u, y_u

class BingoParser:
    def __init__(self, path):   
        self.path = path
        self.r =  (313400.00, 5154800, 400)
        self.r =  (0, 0, 0)
    
    def Parse(self):

        self.cams, self.gcps = self._readCameraFile()
        self.phs = self._readPhotosFile()
        self.pts = self.gcps.copy()
        self._readPhotoPointsFile(self.pts, self.phs)

        self.res_eo = self._readResult(self.phs, self.pts)

        self.cams.values()[0].SetChipSize((11,8))

    def GetData(self):
        return self.cams, self.phs, self.pts, self.gcps

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

                r =  RotMatBLUH(*eo[3:])
                t = np.array([eo[:3]])
                t[0,0] = t[0,0] - self.r[0]
                t[0,1] = t[0,1] - self.r[1]
                t[0,2] = t[0,2] - self.r[2]

                eo = np.hstack((r, t.T))

                phs[ph_id].SetResultExtOr(eo=eo)

            elif "CORD" in l[0]:
                pt_id = int(l[1])
                c = map(float, l[2:5])
                c[0] = c[0] - self.r[0]
                c[1] = c[1] - self.r[1]
                c[2] = c[2] - self.r[2]

                pts[pt_id].SetResultCoords(c)

        return phs

    def _readCameraFile(self):
        f = open(self.path + '/geoin.dat')

        cams = {}
        gcps = {}

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
                control = False
                if 'CHCK' in l[0]:
                    control = True

                gcp_id = int(l[1])
                e = float(l[2]) - self.r[0]
                n = float(l[3]) - self.r[1]
                h = float(l[4]) - self.r[2]

                c = np.array([e, n, h])

                gcps[gcp_id] = Point(gcp_id)
                gcps[gcp_id].SetGcp(c, control)
                

        return cams, gcps

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

class FabioDataParser:
    def __init__(self, path):   
        self.path = path

    def Parse(self):

        self.cams = self._readCameraFile()
        self.gcps = self._readGcpsObjCoords()
        self.pts = self.gcps.copy()
        self.phs = self._readPhotosFile(self.pts)
    
    def GetData(self):
        return self.cams, self.phs, self.pts, self.gcps

    def GetResult(self):
        return self.res_eo

    def _readCameraFile(self):
        f = open(self.path + '/io.txt')

        cams = {}
        while 1:
            line = f.readline()
            if not line:
                break
           
            lin=line.strip()
            
            l = lin.split()
            if not l:
                continue
           
            if l[0] == '#':
                continue

            if 'f' in l[0]:
                fl = float(l[1])

            if 'pix_s' in l[0]:
                self.size = float(l[1]) * 100000.0
            #if 'size' in l[0]:
            #    self.size = float(l[1]) / 10
            if 'pp' in l[0]:
                x0 = float(l[1])
                y0 = float(l[2])

                self.cam_m = np.array([[fl * self.size, 0, x0 * self.size],
                                      [0, fl * self.size, y0 * self.size],
                                      [0, 0, 1]])

                cam_m = np.array([[fl, 0, 0],
                                      [0, fl, 0],
                                      [0, 0, 1.0]])

            if 'dist' in l[0]:
                self.distor = np.array(map(float, l[1:]), dtype=float)
                distor = np.zeros(4)

        self.cam = cams[l[1]] = Camera("cam1", cam_m, distor)

        return cams

    def _readGcpsObjCoords(self):
        f = open(self.path + '/obj_coord.txt')
        gcps = {}

        while 1:
            line = f.readline()
            if not line:
                break

            l = line.strip().split()
            if not l:
                continue      
            
            if l[0] == '#':
                continue

            control = False
            if 'CHCK' in l[-1]:
                control = True

            gcp_id = int(l[0])
            e = float(l[2]) 
            n = float(l[3]) 
            h = float(l[4]) 

            c = np.array([e, n, h])

            gcps[gcp_id] = Point(gcp_id)
            gcps[gcp_id].SetGcp(c, control)

        return gcps
    
    def _readPhotosFile(self, pts):
        f = open(self.path + '/img_coord.txt')

        phs = {}
        count = 0

        while 1:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue

            l=line.strip().split()
            if not l:
                continue

            ph_id = int(l[1])
            #TOOD just one cam
            cam = self.cam

            eo = map(float, l)

            if not phs.has_key(ph_id):
                phs[ph_id] = Photo(ph_id, self.cam) 
                               #np.array(eo), TODO
                               #None
            pt_id = int(l[0])

            #TODO
            c1 = np.array([float(l[2]), float(l[3]) ])
            c = undistortPoint(c1, self.cam_m, self.distor)

            c2 = undistortPoints(np.array([c1]), self.cam_m, self.distor)
            c2 = c2 * self.cam_m[0,0] / self.size

            #print c2[0][0][:]
            if not pts.has_key(pt_id):
                pts[pt_id] = Point(pt_id)
            pts[pt_id].AddPhotoPoint(phs[ph_id], c2[0][0][:] )
        return phs

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
