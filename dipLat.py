import numpy as np

import igrf

FACT = 180./np.pi


def dipLat(lat, lon, alt, year=2005.):
    mag = igrf.igrf12syn(0, year, 1, alt, lat, lon)
    print(mag)
    mlat = np.tan(mag[1]/2/FACT)*FACT
    return mlat


if __name__ == '__main__':
    lat = 7.7
    lon = 116
    alt = 0
    year = 2005
    print(dipLat(lat, lon, alt, year))
