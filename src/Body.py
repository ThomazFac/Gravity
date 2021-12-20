import numpy as np


class Body:
    def __init__(self, r=(0, 0, 0), dr=(0, 0, 0), ddr=(0, 0, 0), M=0):
        self.r = r  # Position [m]
        self.dr = dr  # Velocity [m/s]
        self.ddr = ddr  # Acceleration [m/s^2]
        self.m = M  # Mass [kg]

    def __str__(self):
        return '({self.r[0]}, {self.r[1]}, {self.r[2]})'.format(self=self)

    def vel(self):
        v = (self.dr[0] ** 2 + self.dr[1] ** 2 + self.dr[2] ** 2) ** 0.5
        return v

    def acc(self):
        a = (self.ddr[0] ** 2 + self.ddr[1] ** 2 + self.ddr[2] ** 2) ** 0.5
        return a


def mod(self):
    v = (self.r[0] ** 2 + self.r[1] ** 2 + self.r[2] ** 2) ** 0.5
    return v


def ang(self, other):
    # Angle between the position of two bodies (2D for now)
    x0, x1 = self[0], other[0]
    y0, y1 = self[1], other[1]

    t = np.arctan2((y1 - y0), (x1 - x0))

    return t


def dist(self, other):
    # Distance between two bodies
    return ((self.r[0] - other.r[0]) ** 2 + (self.r[1] - other.r[1]) ** 2 + (self.r[2] - other.r[2]) ** 2) ** 0.5

