import numpy as np
from random import random
import struct


class Polygon:

  def __init__(self, point1, point2, point3):
    self.points = [point1, point2, point3]
    self.plane_eq = self.equation_plane(*self.points)
    self.S_polygon = self.calculate_square(*self.points)
    self.normal = self.plane_eq[:3]
  
  def equation_plane(self, point1, point2, point3):
    x1, y1, z1 = point1
    x2, y2, z2 = point2
    x3, y3, z3 = point3
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)
    return a, b, c, d
  
  def calculate_square(self, point1, point2, point3):
    A = self.rasst(point1, point2)
    B = self.rasst(point2, point3)
    C = self.rasst(point1, point3)
    p = (A + B + C)/2
    return (p * (p - A) * (p - B) * (p - C)) ** 0.5
  
  def rasst(self, point1, point2):
    return ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2) ** 0.5
  
  def is_in_polygon(self, point):
    S1, S2, S3 = self.calculate_square(point, self.points[0], self.points[1]), self.calculate_square(point, self.points[1], self.points[2]), self.calculate_square(point, self.points[0], self.points[2])
    if abs(S1 + S2 + S3 - self.S_polygon) <= 0.000001:
      return True
    return False

  def __str__(self):
    return "{0}x + {1}y + {2}z + {3}=0".format(*self.plane_eq)


class Checker:

  def __init__(self, filename):
    normal, vertex1, vertex2, vertex3 = self.BinarySTL(filename)
    self.model = [Polygon(vertex1[i].tolist(), vertex2[i].tolist(), vertex3[i].tolist())  for i in range(len(vertex1))]
    self.eqs = np.array([poly.plane_eq for poly in self.model])
    self.Ss = np.array([poly.S_polygon for poly in self.model])
    self.point1 = np.transpose(np.array([poly.points[0] for poly in self.model]))
    self.point2 = np.transpose(np.array([poly.points[1] for poly in self.model]))
    self.point3 = np.transpose(np.array([poly.points[2] for poly in self.model]))
    self.r12 = np.sum((self.point2 - self.point1) ** 2, axis=0) ** 0.5
    self.r13 = np.sum((self.point3 - self.point1) ** 2, axis=0) ** 0.5
    self.r23 = np.sum((self.point2 - self.point3) ** 2, axis=0) ** 0.5
  
  def calculate_square(self, a, b, c):
    p = (a + b + c)/2
    pa = p - a
    pb = p - b
    pc = p - c
    return (p * pa * pb * pc)**0.5
  
  def BinarySTL(self, fname):
    fp = open(fname, 'rb')
    Header = fp.read(80)
    nn = fp.read(4)
    Numtri = struct.unpack('i', nn)[0]
    record_dtype = np.dtype([
                   ('normals', np.float32,(3,)),  
                   ('Vertex1', np.float32,(3,)),
                   ('Vertex2', np.float32,(3,)),
                   ('Vertex3', np.float32,(3,)) ,              
                   ('atttr', '<i2',(1,) )
    ])
    data = np.fromfile(fp , dtype = record_dtype , count =Numtri)
    fp.close()
    Normals = data['normals']
    Vertex1= data['Vertex1']
    Vertex2= data['Vertex2']
    Vertex3= data['Vertex3']
    return Normals, Vertex1, Vertex2, Vertex3
  
  def check(self, x, y, z):
    all = 0
    mxon = 0
    coords = np.array([x, y, z, 1]).reshape(4, 1)
    for i in range(1000):
      upperz, lowerz = 0, 0
      m, n, p = (random() - 0.5) * 20, (random() - 0.5) * 20, (random() - 0.5) * 20
      napr = np.array([m, n, p, 0]).reshape(4, 1)
      t = -1 * (np.matmul(self.eqs, coords))/(np.matmul(self.eqs, napr))
      t = t.reshape(-1)
      is_on_plane = t[abs(t) == 0]
      mxon = max(mxon, len(is_on_plane))
      coords0 = np.zeros((11, t.shape[0]))
      coords0[0] = x + m * t
      coords0[1] = y + n * t
      coords0[2] = z + p * t
      coords0[3] = np.sum((coords0[:3] - self.point1) ** 2, axis=0) ** 0.5
      coords0[4] = np.sum((coords0[:3] - self.point2) ** 2, axis=0) ** 0.5
      coords0[5] = np.sum((coords0[:3] - self.point3) ** 2, axis=0) ** 0.5
      coords0[6] = self.calculate_square(coords0[3], coords0[4], self.r12)
      coords0[7] = self.calculate_square(coords0[3], coords0[5], self.r13)
      coords0[8] = self.calculate_square(coords0[4], coords0[5], self.r23)
      coords0[9] = (coords0[6] + coords0[7] + coords0[8]) - self.Ss
      coords0[10] = coords0[9] <= 0.000001
      maskup = coords0[2] >= z
      masklow = coords0[2] <= z
      coords0 = np.transpose(coords0)
      ups = coords0[maskup]
      lows = coords0[masklow]
      upperz = np.sum(ups[:, 10])
      lowerz = np.sum(lows[:, 10])
      all += 1 - 2 * ((upperz % 2) == 0)
      all += 1 - 2 * ((lowerz % 2) == 0)
    if abs(all) <= 6 and mxon > 0:
      return True
    if all >= 0:
      return True
    return False
