#!/usr/bin/env python
import sys
from datetime import datetime
from math import radians, degrees, sin, cos, sqrt, atan2

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from obspy.geodetics.base import gps2dist_azimuth

from .decorators import retry

# Python 2
try:
    from urllib2 import urlopen, HTTPError
# Python 3
except ImportError:
    from urllib.request import urlopen, HTTPError


def progress_bar(count, total, status=''):
    """Show progress bar for the desired task."""
    barlen = 30
    filledlen = int(round(barlen * count / float(total)))

    percents = round(100. * count / float(total), 1)
    pbar = '#' * filledlen + '-' * (barlen - filledlen)

    if count < total:
        sys.stdout.write('[%s] %s%%  %s\r' % (pbar, percents, status))
        sys.stdout.flush()
    else:
        sys.stdout.write('[%s] %s%% %s\n' % (pbar, percents,
            'Done downloading data' + ' '*(len(status)-21)))
        sys.stdout.flush()


@retry(HTTPError, tries=5, delay=0.5, backoff=1)
def access_url(url):
    """Return list of lines in url; if HTTPError, repeat."""
    for line in urlopen(url).readlines():
        yield line.rstrip()


def format_time(ogtime):
    """Increase readability of given catalog times."""
    return pd.Timestamp(' '.join(ogtime.split('T'))[:-1])


def to_epoch(ogtime):
    """Convert formatted time to Unix/epoch time."""
    fstr = '%Y-%m-%dT%H:%M:%S.%fZ'
    epoch = datetime(1970, 1, 1)
    finaltime = (datetime.strptime(ogtime, fstr) - epoch).total_seconds()

    return finaltime


def get_azs_and_dists(cat1, cat2, cat1mids, cat2mids):
    """Calculate azimuths for all matches between two catalogs."""
    azimuths = []
    dists = []

    for idx, eid in enumerate(cat1mids):
        mask1 = cat1['id'] == eid
        mask2 = cat2['id'] == cat2mids[idx]

        lat1 = cat1[mask1].iloc[0]['latitude']
        lon1 = cat1[mask1].iloc[0]['longitude']
        lat2 = cat2[mask2].iloc[0]['latitude']
        lon2 = cat2[mask2].iloc[0]['longitude']

        dist, azi, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
        azimuths.append(azi)
        dists.append(dist/1000.)

    return azimuths, dists


def round2bin(number, binsize, direction):
    """Round number to nearest histogram bin edge (either 'up' or 'down')."""
    if direction == 'down':
        return number - (number % binsize)
    elif direction == 'up':
        return number - (number % binsize) + binsize


def get_map_bounds(cat1, cat2=None):
    """Generate map bounds and gridsize."""
    if cat2 is not None:
        minlat1, maxlat1 = cat1['latitude'].min(), cat1['latitude'].max()
        minlon1, maxlon1 = cat1['longitude'].min(), cat1['longitude'].max()
        minlat2, maxlat2 = cat2['latitude'].min(), cat2['latitude'].max()
        minlon2, maxlon2 = cat2['longitude'].min(), cat2['longitude'].max()

        minlat, maxlat = min(minlat1, minlat2), max(maxlat1, maxlat2)
        minlon, maxlon = min(minlon1, minlon2), max(maxlon1, maxlon2)
    else:
        minlat, maxlat = cat1['latitude'].min(), cat1['latitude'].max()
        minlon, maxlon = cat1['longitude'].min(), cat1['longitude'].max()

    latdiff, londiff = (maxlat-minlat) / 5., (maxlon-minlon) / 5.
    lllat, lllon = max(minlat-latdiff, -90), max(minlon-londiff, -180)
    urlat, urlon = min(maxlat+latdiff, 90), min(maxlon+londiff, 180)

    if (lllon < 175) and (urlon > 175) and (len(cat1[
            cat1['longitude'].between(-100, 100)]) == 0):
        lllon = cat1[cat1['longitude'] > 0].min()['longitude']
        urlon = 360 + cat1[cat1['longitude'] < 0].max()['longitude']
        clon = 180
    else:
        clon = 0

    gridsize = max(urlat-lllat, urlon-lllon) / 30.
    hgridsize, tgridsize = gridsize / 2., gridsize / 10.

    return lllat, lllon, urlat, urlon, gridsize, hgridsize, tgridsize, clon


def round2center(num, gridsize):
    """Round number to nearest grid-square center."""
    hgridsize = gridsize / 2.

    return num - (num % gridsize) + hgridsize


def round2lon(num):
    """Round number to nearest timezone longitude."""
    return 15 * round(num / 15.)


def add_centers(catalog, gridsize):
    """Add corresponding centers to catalog."""
    zippedlatlon = list(zip(round2center(catalog['latitude'], gridsize),
                            round2center(catalog['longitude'], gridsize)))
    newcat = catalog.reset_index()
    newcat.loc[:, 'center'] = pd.Series(zippedlatlon)
    newcat = newcat.set_index('index')
    newcat.index.names = ['']

    return newcat


def group_lat_lons(catalog, minmag=-5):
    """Group detections by nearest grid-square center and return min/max of
    counts."""
    if not catalog['mag'].isnull().all():
        magmask = catalog['mag'] >= minmag
        groupedlatlons = catalog[magmask].groupby('center')
        groupedlatlons = groupedlatlons.count().sort_index()
    elif catalog['mag'].isnull().all() and (minmag != -5):
        groupedlatlons = catalog.groupby('center').count().sort_index()
        print("No magnitude data in catalog - plotting all events")
    else:
        groupedlatlons = catalog.groupby('center').count().sort_index()

    groupedlatlons = groupedlatlons[['id']]
    groupedlatlons.columns = ['count']
    cmin = min(list(groupedlatlons['count']) or [0])
    cmax = max(list(groupedlatlons['count']) or [0])

    return groupedlatlons, cmin, cmax


def range2rgb(rmin, rmax, numcolors):
    """Create a list of red RGB values using colmin and colmax with numcolors
    number of colors."""
    colors = np.linspace(rmax/255., rmin/255., numcolors)
    colors = [(min(1, x), max(0, x-1), max(0, x-1)) for x in colors]

    return colors


def draw_grid(lats, lons, col, alpha=1):
    """Draw rectangle with vertices given in degrees."""
    latlons = list(zip(lons, lats))
    poly = Polygon(latlons, facecolor=col, alpha=alpha, edgecolor='k',
                   zorder=11, transform=ccrs.PlateCarree())
    plt.gca().add_patch(poly)


def trim_times(cat1, cat2, otwindow=16):
    """Trim catalogs so they span the same time window."""
    mintime = max(cat1['time'].min(), cat2['time'].min())
    maxtime = min(cat1['time'].max(), cat2['time'].max())
    adjmin = mintime - otwindow
    adjmax = maxtime + otwindow

    cat1trim = cat1[cat1['time'].between(adjmin, adjmax, inclusive=True)
                   ].copy()
    cat2trim = cat2[cat2['time'].between(adjmin, adjmax, inclusive=True)
                   ].copy()
    cat1trim = cat1trim.reset_index(drop=True)
    cat2trim = cat2trim.reset_index(drop=True)

    return cat1trim, cat2trim


def WW2000(mcval, mags, binsize):
    """Wiemer and Wyss (2000) method for determining a and b values."""
    mags = mags[~np.isnan(mags)]
    mags = np.around(mags, 1)
    mc_vec = np.arange(mcval-1.5, mcval+1.5+binsize/2., binsize)
    max_mag = max(mags)
    corr = binsize / 2.
    bvalue = np.zeros(len(mc_vec))
    std_dev = np.zeros(len(mc_vec))
    avalue = np.zeros(len(mc_vec))
    rval = np.eros(len(mc_vec))

    for idx in range(len(mc_vec)):
        mval = mags[mags >= mc_vec[idx]-0.001]
        mag_bins_edges = np.arange(mc_vec[idx]-binsize/2., max_mag+binsize,
                                   binsize)
        mag_bins_centers = np.arange(mc_vec[idx], max_mag+binsize/2., binsize)

        cdf = np.zeros(len(mag_bins_centers))

        for jdx in range(len(cdf)):
            cdf[jdx] = np.count_nonzero(~np.isnan(mags[
                mags >= mag_bins_centers[jdx]-0.001]))

        bvalue[idx] = np.log10(np.exp(1))/(np.average(mval)
            - (mc_vec[idx]-corr))
        std_dev[idx] = bvalue[idx] / sqrt(cdf[0])

        avalue[idx] = np.log10(len(mval)) + bvalue[idx] * mc_vec[idx]
        log_l = avalue[idx] - bvalue[idx] * mag_bins_centers
        lval = 10.**log_l

        bval, _ = np.histogram(mval, mag_bins_edges)
        sval = abs(np.diff(lval))
        rval[idx] = (sum(abs(bval[:-1] - sval)) / len(mval)) * 100

    ind = np.where(rval <= 10)[0]

    if len(ind) != 0:
        idx = ind[0]
    else:
        idx = list(rval).index(min(rval))

    mcval = mc_vec[idx]
    bvalue = bvalue[idx]
    avalue = avalue[idx]
    std_dev = std_dev[idx]
    mag_bins = np.arange(o, max_mag+binsize/2., binsize)
    lval = 10.**(avalue - bvalue * mag_bins)

    return mcval, bvalue, avalue, lval, mag_bins, std_dev


def get_data(catalog, dirname, startyear=2000, endyear=2000, minmag=-5,
             maxmag=12):
    """Download catalog data from earthquake.usgs.gov"""
    year = startyear
    catalog = catalog.lower()
    alldata = []
    catname = catalog if catalog else 'preferred'
    fname = '%s%s-%s.csv' % (catname, startyear, endyear)

    catstring = '&catalog={0}'.format(catname) if catname != 'preferred'\
                 else ''

    bartotal = 12 * (endyear - startyear + 1)
    barcount = 1

    while year <= endyear:
        month = 1
        yeardata = []

        while month <= 12:
            if month in [4, 6, 9, 11]:
                endday = 30
            elif month == 2:
                checkly = (year % 4)

                if checkly == 0:
                    endday = 29
                else:
                    endday = 28
            else:
                endday = 31

            startd = '-'.join([str(year), str(month)])
            endd = '-'.join([str(year), str(month), str(endday)])

            url = ('https://earthquake.usgs.gov/fdsnws/event/1/query.csv'
                   '?starttime={0}-1%2000:00:00&endtime={1}%2023:59:59'
                   '&orderby=time-asc{2}&minmagnitude={3}'
                   '&maxmagnitude={4}').format(startd, endd, catstring,
                   str(minmag), str(maxmag))
            monthdata = list(access_url(url))

            if (month != 1) or (year != startyear):
                del monthdata[0]

            yeardata.append(monthdata)

            progress_bar(barcount, bartotal, 'Downloading data ...')
            barcount += 1
            month += 1

        alldata.append(yeardata)
        year += 1

    # Flatten list
    alldata = [item for sublist in alldata for item in sublist]
    alldata = [item for sublist in alldata for item in sublist]

    with open('%s/%s' % (dirname, fname), 'w') as openfile:
        for event in alldata:
            openfile.write('%s\n' % event.decode())
    alldatadf = pd.read_csv('%s/%s' % (dirname, fname))

    return alldatadf

