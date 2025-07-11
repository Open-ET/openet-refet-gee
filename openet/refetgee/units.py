import math


def deg2rad(deg):
    return deg * math.pi / 180.0


def rad2deg(rad):
    return rad * 180.0 / math.pi


def c2f(c):
    return c * (9.0 / 5) + 32


def f2c(f):
    return (f - 32) * (5.0 / 9)
