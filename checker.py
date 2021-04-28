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
        p = (A + B + C) / 2
        return (p * (p - A) * (p - B) * (p - C)) ** 0.5

    def rasst(self, point1, point2):
        return ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2) ** 0.5

    def is_in_polygon(self, point):
        S1, S2, S3 = self.calculate_square(point, self.points[0], self.points[1]), self.calculate_square(point,
                                                                                                         self.points[1],
                                                                                                         self.points[
                                                                                                             2]), self.calculate_square(
            point, self.points[0], self.points[2])
        if abs(S1 + S2 + S3 - self.S_polygon) <= 0.000001:
            return True
        return False

    def is_parallel(self, napr):
        scalar = 0
        for i in range(3):
            scalar += napr[i] * self.normal[i]
        if scalar <= 0.00001:
            return True
        return False

    def __str__(self):
        return "{0}x + {1}y + {2}z + {3}=0".format(*self.plane_eq)


class Checker:

    def __init__(self, filename):
        normal, vertex1, vertex2, vertex3 = self.BinarySTL(filename)
        self.model = [Polygon(vertex1[i].tolist(), vertex2[i].tolist(), vertex3[i].tolist()) for i in
                      range(len(vertex1))]

    def BinarySTL(self, fname):
        fp = open(fname, 'rb')
        Header = fp.read(80)
        nn = fp.read(4)
        Numtri = struct.unpack('i', nn)[0]
        record_dtype = np.dtype([
            ('normals', np.float32, (3,)),
            ('Vertex1', np.float32, (3,)),
            ('Vertex2', np.float32, (3,)),
            ('Vertex3', np.float32, (3,)),
            ('atttr', '<i2', (1,))
        ])
        data = np.fromfile(fp, dtype=record_dtype, count=Numtri)
        fp.close()
        Normals = data['normals']
        Vertex1 = data['Vertex1']
        Vertex2 = data['Vertex2']
        Vertex3 = data['Vertex3']
        return Normals, Vertex1, Vertex2, Vertex3

    def check(self, x, y, z):
        all = 0
        for i in range(100):
            upperz, lowerz = 0, 0
            m, n, p = (random() - 0.5) * 20, (random() - 0.5) * 20, (random() - 0.5) * 20
            for poly in self.model:
                eq = poly.plane_eq
                t = -1 * (x * eq[0] + y * eq[1] + z * eq[2] + eq[3]) / (m * eq[0] + n * eq[1] + p * eq[2])
                x0, y0, z0 = x + m * t, y + n * t, z + p * t
                if poly.is_in_polygon([x0, y0, z0]):
                    if z0 >= z:
                        upperz += 1
                    else:
                        lowerz += 1
            all += 1 - 2 * ((upperz % 2) == 0)
            all += 1 - 2 * ((lowerz % 2) == 0)
        if all >= 0:
            return True
        return False
