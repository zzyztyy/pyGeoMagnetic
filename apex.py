import numpy as np

from coordinate import cartesian2geocentric, geocentric2cartesian,\
    geodetic2geocentric, rotateVector, geocentric2geodetic
from pyIGRF.calculate import igrf12syn
from pyIGRF.loadCoeffs import getCoeffs


FACT = 180./np.pi
R = 6371.2


def itrace(Y, YOLD, YAPX, YP, BX, BY, BZ, BB, SGN, DS, NSTP):
    """
    Follow a geomagnetic field line until passing its apex
    Cartesian component magnetic field (partial) derivitives steer the trace

    This uses the 4-point Adams formula after initialization.
    First 7 iterations advance point by 3 steps.

    INPUTS:
     YAPX  = Matrix of cartesian coordinates of the three points
             about the apex.  Set in itrace.

    FLDCOMD has geomagnetic field at current trace point:
     BX    = X component (nT)
     BY    = Y component (nT)
     BZ    = Z component (nT)
     BB    = Magnitude   (nT)

    itrace has field line tracing variables determined in findApex:
     NSTP  = Step count.
     Y     = Array containing current tracing point cartesian coordinates.
     YOLD  = Array containing previous tracing point cartesian coordinates.
     SGN   = Determines direction of trace.
     DS    = Step size (arc length in km).

    REFERENCES:
     Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
     Greenbelt, Maryland
    """
    YP[0][3] = SGN * BX / BB
    YP[1][3] = SGN * BY / BB
    YP[2][3] = SGN * BZ / BB
    D2 = DS / 2.
    D6 = DS / 6.
    D12 = DS / 12.
    D24 = DS / 24.
    if NSTP <= 7:
        for I in range(0, 3):
            if NSTP == 1:
                YP[I][0] = YP[I][3]
                YOLD[I] = Y[I]
                YAPX[I][0] = Y[I]
                Y[I] = YOLD[I] + DS * YP[I][0]
            elif NSTP == 2:
                YP[I][1] = YP[I][3]
                Y[I] = YOLD[I] + D2 * (YP[I][1] + YP[I][0])
            elif NSTP == 3:
                Y[I] = YOLD[I] + D6 * (2. * YP[I][3] + YP[I][1] + 3. * YP[I][0])
            elif NSTP == 4:
                YP[I][1] = YP[I][3]
                YAPX[I][1] = Y[I]
                YOLD[I] = Y[I]
                Y[I] = YOLD[I] + D2 * (3.0 * YP[I][1] - YP[I][0])
            elif NSTP == 5:
                Y[I] = YOLD[I] + D12 * (5. * YP[I][3] + 8. * YP[I][1] - YP[I][0])
            elif NSTP == 6:
                YP[I][2] = YP[I][3]
                YOLD[I] = Y[I]
                YAPX[I][2] = Y[I]
                Y[I] = YOLD[I] + D12 * (23. * YP[I][2] - 16. * YP[I][1] + 5. * YP[I][0])
            elif NSTP == 7:
                YAPX[I][0] = YAPX[I][1]
                YAPX[I][1] = YAPX[I][2]
                Y[I] = YOLD[I] + D24 * (9. * YP[I][3] + 19. * YP[I][2] - 5. * YP[I][1] + YP[I][0])
                YAPX[I][2] = Y[I]

    else:
        for I in range(0, 3):
            YAPX[I][0] = YAPX[I][1]
            YAPX[I][1] = Y[I]
            YOLD[I] = Y[I]
            Y[I] = YOLD[I] + D24 * (55. * YP[I][3] - 59. * YP[I][2] + 37. * YP[I][1] - 9. * YP[I][0])
            YAPX[I][2] = Y[I]
            for J in range(0, 3):
                YP[I][J] = YP[I][J + 1]


def northPole(date):
    """
    Calculate the position of north pole at the first order approximation.
    :param date: years(float year)
    :return: lat(float deg) lon(float deg)
    """
    g, h = getCoeffs(date)
    colat = np.arccos(g[1][0]/np.sqrt(g[1][0]**2 + g[1][1]**2 + h[1][1]**2))
    elong = np.arctan2(h[1][1], g[1][1])
    return colat*FACT-90, elong*FACT-180


def traceToApex(lat, lon, alt, date=2005.):
    """
    Follow a geomagnetic field line until passing its apex.
    :param lat: latitude of start position(float deg)
    :param lon: longitude of start position(float deg)
    :param alt: altitude of start position(float deg)
    :param date: years(year)
    :return: trace: dots of the field line(list([float, float, float]) km)
    """
    nlat, nlon = northPole(date)
    nlat, nlon = nlat/FACT, nlon/FACT
    ctp = np.cos(np.pi/2-nlat)
    stp = np.sin(np.pi/2-nlat)

    gccolat, plon, gcrho = geodetic2geocentric(np.pi/2-lat/FACT, alt)
    gclat, gclon = np.pi/2-gccolat, lon/FACT
    x0, y0, z0 = geocentric2cartesian(gclat, gclon, gcrho)
    bx, by, bz, bb = igrf12syn(0, date, 2, gcrho, gclat * FACT, gclon * FACT)

    trace = []
    step = 0
    sgn = -np.sign(bz)
    arrive = False
    Y = [x0, y0, z0]
    YOLD = [0, 0, 0]
    YAPX= [[0., 0., 0.],
           [0., 0., 0.],
           [0., 0., 0.]]
    YP = [[0., 0., 0., 0.],
          [0., 0., 0., 0.],
          [0., 0., 0., 0.]]

    while not arrive and step < 100:
        stngml = ctp*np.sin(gclat)+stp*np.cos(gclat)*np.cos(gclon-nlon)
        cgml2 = max(1-stngml**2, 0.25)
        DS = gcrho*0.06/cgml2-370

        bx, by, bz, bb = igrf12syn(0, date, 2, gcrho, gclat * FACT, gclon * FACT)
        # lstBdown = [lstBdown[1], lstBdown[2], bz]
        Bx, By, Bz = rotateVector([bx, by, bz], gclon, gclat)

        itrace(Y, YOLD, YAPX, YP, Bx, By, Bz, bb, sgn, DS, step)
        step += 1
        gclat, gclon, gcrho = cartesian2geocentric(Y[0], Y[1], Y[2])
        trace.append(Y[:])

        if step >= 7:
            RC = np.sqrt(YAPX[0][2] ** 2 + YAPX[1][2] ** 2 + YAPX[2][2] ** 2)
            RP = np.sqrt(YAPX[0][1] ** 2 + YAPX[1][1] ** 2 + YAPX[2][1] ** 2)
            arrive = RC < RP

    # interpolate to where Bdown/B < 0.00002 to find cartesian coordinates at dip equator
    dot = [trace[-3], trace[-2], trace[-1]]
    lstBdown = [0., 0., 0.]
    step = 0
    while abs(bz/bb) > 0.000002:
        if step >= 3:
            dot1, dot2, dot3 = dot[-3], dot[-2], dot[-1]
            Ax = secDegInterpolate(lstBdown[0], lstBdown[1], lstBdown[2], dot1[0], dot2[0], dot3[0], 0.)
            Ay = secDegInterpolate(lstBdown[0], lstBdown[1], lstBdown[2], dot1[1], dot2[1], dot3[1], 0.)
            Az = secDegInterpolate(lstBdown[0], lstBdown[1], lstBdown[2], dot1[2], dot2[2], dot3[2], 0.)
            dot.append([Ax, Ay, Az])
            gclat, gclon, gcrho = cartesian2geocentric(Ax, Ay, Az)
            bx, by, bz, bb = igrf12syn(0, date, 2, gcrho, gclat * FACT, gclon * FACT)
            lstBdown = [lstBdown[1], lstBdown[2], bz]
        else:
            Ax, Ay, Az = dot[step]
            gclat, gclon, gcrho = cartesian2geocentric(Ax, Ay, Az)
            bx, by, bz, bb = igrf12syn(0, date, 2, gcrho, gclat * FACT, gclon * FACT)
            lstBdown = [lstBdown[1], lstBdown[2], bz]
        step += 1
    trace[-1] = dot[-1]
    return trace, sgn


def secDegInterpolate(x1, x2, x3, y1, y2, y3, xfit):
    """second degree interpolation used by findApex"""
    x12 = x1-x2
    x13 = x1-x3
    x23 = x2-x3
    xf1 = xfit-x1
    xf2 = xfit-x2
    xf3 = xfit-x3
    yfit = (y1*x23*xf2*xf3-y2*x13*xf1*xf3+y3*x12*xf1*xf2)/(x12*x13*x23)
    return yfit


def qdCoordination(start, apex, northPoleLat, northPoleLon, sgn):
    """
    Calculate the magnetic coordination according to Apex position.
    :param start: position of start [x, y, z] (list(float), km)
    :param apex: position of Apex [x, y, z] (list(float), km)
    :param northPoleLat: latitude of north pole (float, rad)
    :param northPoleLon: longitude of north pole (float, rad)
    :param sgn: minus for south, positive for north (-1 or 1)
    :return mlat: magnetic latitude (float, deg)
            mlon: magnetic longitude (float, deg)
    """
    sgclat, sgclon, sr = cartesian2geocentric(start[0], start[1], start[2])
    slat, sh = geocentric2geodetic(sgclat, sr)
    agclat, agclon, ar = cartesian2geocentric(apex[0], apex[1], apex[2])
    alat, ah = geocentric2geodetic(agclat, ar)
    mlat = np.arccos(min(np.sqrt((R+sh)/(R+ah)), 1))*sgn*FACT
    # print(alat*FACT, agclon*FACT, ah)

    xlon = np.arctan2(apex[1], apex[0])
    elon = northPoleLon
    ang = xlon-elon
    cang = np.cos(ang)
    sang = np.sin(ang)
    r = ar
    cte = apex[2]/r
    ste = np.sqrt(1.-cte*cte)
    ctp = np.cos(np.pi/2-northPoleLat)
    stp = np.sin(np.pi/2-northPoleLat)
    stfcpa = ste*ctp*cang-cte*stp
    stfspa = sang*ste
    mlon = np.arctan2(stfspa, stfcpa)*FACT % 360.
    return mlat, mlon


def gd2qd(lat, lon, alt=0., date=2005.):
    """
    transform from geodetic to quasi-dipole coordination
    :param lat: latitude (float, deg)
    :param lon: longitude (float, deg)
    :param alt: altitude (float, km)
    :param date: time (float, year)
    :return mlat: magnetic latitude (float, deg)
            mlon: magnetic longitude (float, deg)
    """
    trace, sgn = traceToApex(lat, lon, alt, date)
    nlat, nlon = northPole(date)
    mlat, mlon = qdCoordination(trace[0], trace[-1], nlat/FACT, nlon/FACT, -sgn)
    return mlat, mlon


def lineToApex(mlat, mlon, alt, date=2005.):
    nlat, nlon = northPole(date)
    nlat, nlon = nlat / FACT, nlon / FACT
    cb = np.cos(np.pi / 2 - nlat)
    sb = np.sin(np.pi / 2 - nlat)

    A = np.pi - mlon
    sA = np.sin(A)
    cA = np.cos(A)

    ha = (R + alt)/(np.cos(mlat)**2) - R

    arrive = False
    north_a = np.pi/2 - nlat
    south_a = np.pi/2 + nlat
    alat, alon = 0., 0.
    # 控制a，二分迭代找apex
    while not arrive:
        a = (north_a+south_a)/2
        bz, bb, alat, alon = tempB(a, ha, cb, sb, cA, sA, nlon, date)
        # print(tempB(north_a, ha, cb, sb, cA, sA, nlon, date)[0], tempB(south_a, ha, cb, sb, cA, sA, nlon, date)[0])
        # print(a*FACT, north_a*FACT, south_a*FACT, alon*FACT, bz/bb)
        # print('.')
        # print(bz, bb, bz/bb)
        if abs(bz/bb) < 0.000002:
            arrive = True
        else:
            if bz > 0:
                north_a = a
            else:
                south_a = a
    return alat, alon, ha


def tempB(a, h, cb, sb, cA, sA, nlon, date):
    sgn = np.sign(sA)
    sa = np.sin(a)
    ca = np.cos(a)
    sB = sb*sA/sa
    cB = np.sqrt(1-sB*sB)
    cC = (sA*sB*ca*cb-cA*cB)/(1-sb*sb*sA*sA)
    C = np.arccos(cC)*FACT

    lon = (sgn*np.arccos(cC) + nlon) % (2*np.pi)
    lat = np.pi/2 - a

    bx, by, bz, bb = igrf12syn(0, date, 1, h, lat*FACT, lon*FACT)
    return bz, bb, lat, lon


# def itrace_r(Y, YOLD, YAPX, YP, BX, BY, BZ, BB, SGN, DS, NSTP):
#     YP[0][3] = SGN * BX / BB
#     YP[1][3] = SGN * BY / BB
#     YP[2][3] = SGN * BZ / BB
#     D2 = DS / 2.
#     D6 = DS / 6.
#     D12 = DS / 12.
#     D24 = DS / 24.
#     if NSTP <= 7:
#         for I in range(0, 3):
#             if NSTP == 1:
#                 YP[I][0] = YP[I][3]
#                 YOLD[I] = Y[I]
#                 YAPX[I][0] = Y[I]
#                 Y[I] = YOLD[I] + DS * YP[I][0]
#             elif NSTP == 2:
#                 YP[I][1] = YP[I][3]
#                 Y[I] = YOLD[I] + D2 * (YP[I][1] + YP[I][0])
#             elif NSTP == 3:
#                 Y[I] = YOLD[I] + D6 * (2. * YP[I][3] + YP[I][1] + 3. * YP[I][0])
#             elif NSTP == 4:
#                 YP[I][1] = YP[I][3]
#                 YAPX[I][1] = Y[I]
#                 YOLD[I] = Y[I]
#                 Y[I] = YOLD[I] + D2 * (3.0 * YP[I][1] - YP[I][0])
#             elif NSTP == 5:
#                 Y[I] = YOLD[I] + D12 * (5. * YP[I][3] + 8. * YP[I][1] - YP[I][0])
#             elif NSTP == 6:
#                 YP[I][2] = YP[I][3]
#                 YOLD[I] = Y[I]
#                 YAPX[I][2] = Y[I]
#                 Y[I] = YOLD[I] + D12 * (23. * YP[I][2] - 16. * YP[I][1] + 5. * YP[I][0])
#             elif NSTP == 7:
#                 YAPX[I][0] = YAPX[I][1]
#                 YAPX[I][1] = YAPX[I][2]
#                 Y[I] = YOLD[I] + D24 * (9. * YP[I][3] + 19. * YP[I][2] - 5. * YP[I][1] + YP[I][0])
#                 YAPX[I][2] = Y[I]
#
#         # if NSTP == 6 or NSTP == 7:  # signal if apex passed
#         #     RC = np.sqrt(YAPX[0][2]**2 + YAPX[1][2]**2 + YAPX[2][2]**2)
#         #     RP = np.sqrt(YAPX[0][1]**2 + YAPX[1][1]**2 + YAPX[2][1]**2)
#         #     if RC < RP:
#         #         return True
#     else:  # NSTP > 7
#         for I in range(0, 3):
#             YAPX[I][0] = YAPX[I][1]
#             YAPX[I][1] = Y[I]
#             YOLD[I] = Y[I]
#             Y[I] = YOLD[I] + D24 * (55. * YP[I][3] - 59. * YP[I][2] + 37. * YP[I][1] - 9. * YP[I][0])
#             YAPX[I][2] = Y[I]
#             for J in range(0, 3):
#                 YP[I][J] = YP[I][J + 1]
#         # RC = np.sqrt(Y[0]**2 + Y[1]**2 + Y[2]**2)
#         # RP = np.sqrt(YOLD[0]**2 + YOLD[1]**2 + YOLD[2]**2)
#         # if RC < RP:
#         #     return True
#     # return 0


def traceToStart(alat, alon, ha, date, sgn, hs):
    nlat, nlon = northPole(date)
    nlat, nlon = nlat / FACT, nlon / FACT
    ctp = np.cos(np.pi / 2 - nlat)
    stp = np.sin(np.pi / 2 - nlat)

    gccolat, plon, gcrho = geodetic2geocentric(np.pi / 2 - alat, ha)
    gclat, gclon = np.pi / 2 - gccolat, alon
    x0, y0, z0 = geocentric2cartesian(gclat, gclon, gcrho)

    alt = 0.
    trace = []
    step = 0
    arrive = False
    Y = [x0, y0, z0]
    YOLD = [0, 0, 0]
    YAPX = [[0., 0., 0.],
            [0., 0., 0.],
            [0., 0., 0.]]
    YP = [[0., 0., 0., 0.],
          [0., 0., 0., 0.],
          [0., 0., 0., 0.]]

    while not arrive and step < 100:
        stngml = ctp * np.sin(gclat) + stp * np.cos(gclat) * np.cos(gclon - nlon)
        cgml2 = max(1 - stngml ** 2, 0.25)
        DS = gcrho * 0.06 / cgml2 - 370

        bx, by, bz, bb = igrf12syn(0, date, 2, gcrho, gclat * FACT, gclon * FACT)
        Bx, By, Bz = rotateVector([bx, by, bz], gclon, gclat)

        itrace(Y, YOLD, YAPX, YP, Bx, By, Bz, bb, sgn, DS, step)
        step += 1
        gclat, gclon, gcrho = cartesian2geocentric(Y[0], Y[1], Y[2])
        trace.append(Y[:])

        lat, alt = geocentric2geodetic(gclat, gcrho)
        arrive = alt < hs

    # interpolate to where h = hs to find cartesian coordinates at start point
    dot = [trace[-3], trace[-2], trace[-1]]
    h = [0., 0., 0.]
    step = 0
    while abs(alt-hs) > 0.1:
        if step >= 3:
            dot1, dot2, dot3 = dot[-3], dot[-2], dot[-1]
            Ax = secDegInterpolate(h[0], h[1], h[2], dot1[0], dot2[0], dot3[0], hs)
            Ay = secDegInterpolate(h[0], h[1], h[2], dot1[1], dot2[1], dot3[1], hs)
            Az = secDegInterpolate(h[0], h[1], h[2], dot1[2], dot2[2], dot3[2], hs)
            dot.append([Ax, Ay, Az])
            gclat, gclon, gcrho = cartesian2geocentric(Ax, Ay, Az)
            lat, alt = geocentric2geodetic(gclat, gcrho)
            h = [h[1], h[2], alt]
        else:
            Ax, Ay, Az = dot[step]
            gclat, gclon, gcrho = cartesian2geocentric(Ax, Ay, Az)
            lat, alt = geocentric2geodetic(gclat, gcrho)
            h = [h[1], h[2], alt]
        step += 1
    trace[-1] = dot[-1]
    return trace


def qd2gd(mlat, mlon, alt=0., date=2005.):
    sgn = np.sign(mlat)
    alat, alon, ha = lineToApex(mlat/FACT, mlon/FACT, alt, date)
    # print(alat*FACT, alon*FACT, ha)
    trace = traceToStart(alat, alon, ha, date, sgn, alt)
    start = trace[-1]
    gclat, gclon, gcr = cartesian2geocentric(start[0], start[1], start[2])
    lat, alt = geocentric2geodetic(gclat, gcr)
    lon = gclon
    return lat*FACT, lon*FACT, alt
