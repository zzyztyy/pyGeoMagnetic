import numpy as np

from pyIGRF.calculate import igrf12syn

FACT = 180./np.pi


def igrf12(lat, lon, alt=0., year=2005.):
    """
    :return
         D is declination (+ve east)
         I is inclination (+ve down)
         H is horizontal intensity
         X is north component
         Y is east component
         Z is vertical component (+ve down)
         F is total intensity
    """
    X, Y, Z, F = igrf12syn(0, year, 1, alt, lat, lon)
    D = FACT * np.arctan2(Y, X)
    H = np.sqrt(X * X + Y * Y)
    I = FACT * np.arctan2(Z, H)
    return D, I, H, X, Y, Z, F


def igrf12sv(lat, lon, alt=0, year=2005):
    """
         Annual variation
         D is declination (+ve east)
         I is inclination (+ve down)
         H is horizontal intensity
         X is north component
         Y is east component
         Z is vertical component (+ve down)
         F is total intensity
    """
    X, Y, Z, F = igrf12syn(0, year, 1, alt, lat, lon)
    H = np.sqrt(X * X + Y * Y)
    DX, DY, DZ, DF = igrf12syn(1, year, 1, alt, lat, lon)
    DD = (60.0 * FACT * (X * DY - Y * DX)) / (H * H)
    DH = (X * DX + Y * DY) / H
    DS = (60.0 * FACT * (H * DZ - Z * DH)) / (F * F)
    DF = (H * DH + Z * DZ) / F
    return DD, DS, DH, DX, DY, DZ, DF


if __name__ == '__main__':
    DATE = 2005
    ITYPE = 1
    ALT = 300
    CLT = 40
    XLN = 116
    print(igrf12(CLT, XLN, ALT, DATE))
