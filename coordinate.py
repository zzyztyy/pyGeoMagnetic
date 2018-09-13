import numpy as np


def geodetic2geocentric(theta, alt):
    # conversion from geodetic to geocentric coordinates
    # (using the WGS84 spheroid)
    ct = np.cos(theta)
    st = np.sin(theta)
    a2 = 40680631.6
    b2 = 40408296.0
    one = a2 * st * st
    two = b2 * ct * ct
    three = one + two
    rho = np.sqrt(three)
    r = np.sqrt(alt * (alt + 2.0 * rho) + (a2 * one + b2 * two) / three)
    cd = (alt + rho) / r
    sd = (a2 - b2) / rho * ct * st / r
    one = ct
    ct = ct * cd - st * sd
    st = st * cd + one * sd
    gclat = np.arctan2(st, ct)
    gclon = np.arctan2(sd, cd)
    return gclat, gclon, r


def geocentric2cartesian(lat, lon, rho):
    x = rho*np.cos(lat)*np.cos(lon)
    y = rho*np.cos(lat)*np.sin(lon)
    z = rho*np.sin(lat)
    return x, y, z


def cartesian2geocentric(x, y, z):
    rho = np.sqrt(x*x+y*y+z*z)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, np.sqrt(x*x+y*y))
    return lat, lon, rho


def span(vec, alpha, beta):
    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)
    spanMat = np.array([[-sb*ca, -sb*sa, -cb],
                       [-sa,     ca,     0],
                       [cb*ca,   cb*sa,  -sb]])
    return np.dot(spanMat, vec)
