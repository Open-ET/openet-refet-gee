import math

import pytest

import openet.refetgee.units as units


def test_deg2rad(d=30, r=(math.pi / 6)):
    assert units.deg2rad(d) == pytest.approx(r)


def test_rad2deg(r=(math.pi / 6), d=30):
    assert units.rad2deg(r) == pytest.approx(d)


def test_c2f(c=20, f=68):
    assert units.c2f(c) == f


def test_f2c(f=68, c=20):
    assert units.f2c(f) == c
