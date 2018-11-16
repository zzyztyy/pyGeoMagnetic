import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

from pyGeoMagApex import coordinate, apex, dipLat

FACT = 180./np.pi
R = 6371.2
textPath = "TextFile/"


def sphere(center, radius, ax):
    # data
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
    y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
    z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

    # surface plot
    ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')


def draw(trace):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # center and radius
    center = [0, 0, 0]
    radius = 6371.2
    sphere(center, radius, ax)
    ax.plot([x[0] for x in trace], [x[1] for x in trace], [x[2] for x in trace])
    ax.set_zlabel('Z')
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.show()


def testApex():
    latS = -90
    latE = 90
    lonS = -180
    lonE = 180
    alt = 100.
    step = 1
    with open('qd2gd_gd2qd_bias2.txt', 'a') as fout:
        for lat in range(latS, latE, step):
            for lon in range(lonS, lonE, step):
                print(lat+step/2, lon+step/2)
                glat, glon = apex.qd2gd(lat + step / 2, lon + step / 2, alt)
                mlat, mlon = apex.gd2qd(glat, glon, alt)
                b = [lat + step / 2, lon + step / 2, mlat, mlon, glat, glon]
                b = [str(round(x, 2)) for x in b]
                fout.write(' '.join(b) + '\n')
                # print(mlat, mlon)
                # tlat, tlon, talt, trace = apex.qd2gd(mlat, mlon, alt)
                # print(tlat, tlon)


def testDipLat():
    lat = 40
    lon = 116
    alt = 0
    year = 2005
    print(dipLat.dipLat(lat, lon, alt, year))


def testCoordinate():
    lat = 40
    lon = 116
    alt = 0
    # for alt in range(0, 2000, 200):
    #     gclat_list = []
    #     gcr_list = []
    #     plon_list = []
    #     for lat in range(0, 90):
    #         ct = np.cos(np.pi / 2 - lat/FACT)
    #         st = np.sin(np.pi / 2 - lat/FACT)
    #         a2 = 40680631.6
    #         b2 = 40408296.0
    #         rho2 = a2 * st * st + b2 * ct * ct
    #         rho = np.sqrt(rho2)
    #         gclat, plon, gcr = coordinate.geodetic2geocentric(np.pi/2-lat/FACT, alt)
    #         gclat_list.append(gclat)
    #         gcr_list.append(gcr-alt-rho)
    #         plon_list.append(plon)
    #     plt.plot(range(0, 90), gcr_list)
    # plt.show()
    gclat, plon, gcr = coordinate.geodetic2geocentric(np.pi / 2 - lat / FACT, alt)
    print(90-gclat*FACT, gcr)
    lat, alt = coordinate.geocentric2geodetic(np.pi / 2 - gclat, gcr)
    print(lat*FACT, alt)


def testTrace():
    trace, sgn = apex.traceToApex(40, 116, 100, 2005, 77 / FACT, -77 / FACT, 1)
    return trace


def testDS():
    new_time, old_time = [], []
    dlat, dlon = [], []
    for i in range(500):
        print(i)
        lat = np.random.random() * 90
        lon = np.random.random() * 180
        atime = time.time()
        mlat, mlon = apex.gd2qd(lat, lon, ds=1)
        btime = time.time()
        mlat2, mlon2 = apex.gd2qd(lat, lon, ds=0)
        ctime = time.time()
        new_time.append(btime - atime)
        old_time.append(ctime - btime)
        dlat.append(abs(mlat - mlat2))
        dlon.append(abs(mlon - mlon2))
    plt.plot(new_time)
    plt.plot(old_time)
    plt.show()
    print(np.mean(dlat), np.median(dlat), np.max(dlat))
    print(np.mean(dlon), np.median(dlon), np.max(dlon))


def testBias():
    latBin1, lonBin1 = [], []
    latBin2, lonBin2 = [], []
    with open(textPath + 'qd2geo.txt', 'r') as f1:
        text = f1.readlines()
        for a in text:
            b = a.split('#')
            for i in range(360):
                c = b[i]
                lat, lon = c.split()
                latBin1.append(float(lat))
                lonBin1.append(float(lon) % 360)
    # latBin1 = np.array(latBin1).reshape(180, 360)
    # lonBin1 = np.array(lonBin1).reshape(180, 360)

    with open(textPath + 'qd2geoNew.txt', 'r') as f2:
        text = f2.readlines()
        for a in text:
            b = a.split('#')
            for i in range(360):
                c = b[i]
                lat, lon = c.split()
                latBin2.append(float(lat))
                lonBin2.append(float(lon))
    # latBin2 = np.array(latBin2).reshape(180, 360)
    # lonBin2 = np.array(lonBin2).reshape(180, 360)

    latBia = []
    lonBia = []
    for i in range(180*360):
        latBia.append(latBin1[i]-latBin2[i])
        delta = lonBin1[i] - lonBin2[i]
        if delta > 180:
            delta -= 360
        elif delta < -180:
            delta += 360
        lonBia.append(delta)
    latBia = np.array(latBia).reshape(180, 360)
    lonBia = np.array(lonBia).reshape(180, 360)

    plt.figure(figsize=[8, 8])
    plt.subplots_adjust(top=0.98, bottom=0.03, left=0.1, right=1.0, hspace=0.1)
    cmap = 'bwr'
    latMax = 0.5
    lonMax = 0.5
    fontsize = 14
    plt.subplot(2, 1, 1)
    # m = Basemap()
    # m.drawcoastlines()
    plt.title('LatBia', fontsize=fontsize)
    # plt.imshow(latBia, cmap=cmap, vmax=latMax, vmin=-latMax)
    plt.pcolor(range(-180, 181), range(-90, 90), latBia, cmap=cmap, vmax=latMax, vmin=-latMax)
    plt.colorbar()
    plt.xticks(range(-180, 181, 30), fontsize=fontsize)
    plt.yticks(range(-90, 91, 30), fontsize=fontsize)
    plt.xlabel('MagLon', fontsize=fontsize)
    plt.ylabel('MagLat', fontsize=fontsize)

    plt.subplot(2, 1, 2)
    plt.title('LonBia', fontsize=fontsize)
    # m = Basemap()
    # m.drawcoastlines()
    plt.title('LonBia', fontsize=fontsize)
    plt.pcolor(range(-180, 181), range(-90, 90), lonBia, cmap=cmap, vmax=lonMax, vmin=-lonMax)
    plt.colorbar()
    plt.xticks(range(-180, 181, 30), fontsize=fontsize)
    plt.yticks(range(-90, 91, 30), fontsize=fontsize)
    plt.xlabel('MagLon', fontsize=fontsize)
    plt.ylabel('MagLat', fontsize=fontsize)
    plt.show()


def testTime():
    lat, lon, alt = [], [], []
    x, y, z = [], [], []
    t1, t2 = [], []
    rate = []
    with open(textPath + "qd2gd_time.txt", 'r') as f:
        text = f.readlines()
        for a in text:
            b = a.split()
            b = [float(x) for x in b]
            lat.append(b[0])
            lon.append(b[1])
            alt.append(b[2])
            gclat, d, rho = coordinate.geodetic2geocentric(np.pi/2-b[0], b[2])
            xt, yt, zt = coordinate.geocentric2cartesian(gclat, b[1], rho)
            x.append(xt)
            y.append(yt)
            z.append(zt)
            t1.append(b[3])
            t2.append(b[4])
            rate.append(b[5])
    plt.figure(figsize=[10, 5])
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=1.0, hspace=0.1)
    cmap = 'jet'
    plt.scatter(lon, lat, c=rate, cmap=cmap, alpha=0.5, vmax=0.08)
    plt.colorbar()
    plt.xlabel("MagLon")
    plt.ylabel("MagLat")
    title = "Rate_qd2gd"
    plt.title(title)
    plt.axis([-180, 180, -90, 90])
    # plt.savefig(title+'.png')
    plt.show()

    # plt.scatter(rate, alt)
    # plt.show()


def testConsist():
    tlat = []
    tlon = []
    mlat = []
    mlon = []
    with open(textPath + "qd2gd_gd2qd_bias2.txt", 'r') as f:
        text = f.readlines()
        for a in text:
            b = a.split()
            b = [float(x) for x in b]
            mlat.append(b[0])
            mlon.append(b[1])
            tlat.append(b[2]-b[0])
            dlon = b[3] - b[1]
            if dlon > 180:
                dlon -= 360
            elif dlon < -180:
                dlon += 360
            tlon.append(dlon)
    latBia = np.array(tlat).reshape(180, 360)
    lonBia = np.array(tlon).reshape(180, 360)

    plt.figure(figsize=[8, 8])
    plt.subplots_adjust(top=0.98, bottom=0.03, left=0.1, right=1.0, hspace=0.1)
    cmap = 'bwr'
    latMax = 0.1
    lonMax = 0.1
    fontsize = 14
    plt.subplot(2, 1, 1)
    # m = Basemap()
    # m.drawcoastlines()
    plt.title('LatConsist', fontsize=fontsize)
    # plt.imshow(latBia, cmap=cmap, vmax=latMax, vmin=-latMax)
    plt.pcolor(range(-180, 181), range(-90, 90), latBia, cmap=cmap, vmax=latMax, vmin=-latMax)
    plt.colorbar()
    plt.xticks(range(-180, 181, 30), fontsize=fontsize)
    plt.yticks(range(-90, 91, 30), fontsize=fontsize)
    plt.xlabel('MagLon', fontsize=fontsize)
    plt.ylabel('MagLat', fontsize=fontsize)

    plt.subplot(2, 1, 2)
    plt.title('LonBia', fontsize=fontsize)
    # m = Basemap()
    # m.drawcoastlines()
    plt.title('LonConsist', fontsize=fontsize)
    plt.pcolor(range(-180, 181), range(-90, 90), lonBia, cmap=cmap, vmax=lonMax, vmin=-lonMax)
    plt.colorbar()
    plt.xticks(range(-180, 181, 30), fontsize=fontsize)
    plt.yticks(range(-90, 91, 30), fontsize=fontsize)
    plt.xlabel('MagLon', fontsize=fontsize)
    plt.ylabel('MagLat', fontsize=fontsize)
    plt.show()


if __name__ == '__main__':
    print('start')
    testConsist()
