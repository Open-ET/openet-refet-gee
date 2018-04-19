import math


def _deg2rad(deg):
    return deg * math.pi / 180.0

def _rad2deg(rad):
    return rad * 180.0 / math.pi

def _c2f(c):
    return c * (9.0 / 5) + 32

def _f2c(f):
    return (f - 32) * (5.0 / 9)
