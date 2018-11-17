#!/usr/bin/env python
"""Code for creating figures comparing two catalogs spanning the same time
frame. Run `QCmulti.py -h` for command line options.
"""
import os
import sys
import errno
import argparse
import time
import shutil
from datetime import datetime
from math import sqrt, degrees, radians, sin, cos, atan2, pi, ceil

import markdown
from scipy import stats
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth

# Python 2
try:
    from urllib2 import urlopen, HTTPError
# Python 3
except ImportError:
    from urllib.request import urlopen, HTTPError

import QCutils as qcu
from decorators import retry, printstatus


###############################################################################
###############################################################################
###############################################################################


@printstatus('Creating basic catalog summary')
def basic_cat_sum(catalog, catname, dirname):
    """Gather basic catalog summary statistics."""
    lines = []

    lines.append('Catalog name: %s\n\n' % catname.upper())

    lines.append('First date in catalog: %s\n' % catalog['time'].min())
    lines.append('Last date in catalog: %s\n\n' % catalog['time'].max())

    lines.append('Total number of events: %s\n\n' % len(catalog))

    lines.append('Minimum latitude: %s\n' % catalog['latitude'].min())
    lines.append('Maximum latitude: %s\n' % catalog['latitude'].max())
    lines.append('Minimum longitude: %s\n' % catalog['longitude'].min())
    lines.append('Maximum longitude: %s\n\n' % catalog['longitude'].max())

    lines.append('Minimum depth: %s\n' % catalog['depth'].min())
    lines.append('Maximum depth: %s\n' % catalog['depth'].max())
    lines.append('Number of 0 km depth events: %s\n'
                 % len(catalog[catalog['depth'] == 0]))
    lines.append('Number of NaN depth events: %s\n\n'
                 % len(catalog[pd.isnull(catalog['depth'])]))

    lines.append('Minimum magnitude: %s\n' % catalog['mag'].min())
    lines.append('Maximum magnitude: %s\n' % catalog['mag'].max())
    lines.append('Number of 0 magnitude events: %s\n'
                 % len(catalog[catalog['mag'] == 0]))
    lines.append('Number of NaN magnitude events: %s'
                 % len(catalog[pd.isnull(catalog['mag'])]))

    with open('%s_summary.txt' % catname, 'w') as sumfile:
        for line in lines:
            sumfile.write(line)


@printstatus('Creating summary of comparison criteria and statistics')
def comp_criteria(cat1, cat1name, cat1mids, cat2, cat2name, cat2mids, dirname,
                  otwindow=16, distwindow=100):
    """Trim catalogs and summarize comparison criteria/statistics."""
    lines = []

    nummatches = len(cat1mids)
    unq1 = len(cat1) - len(cat1mids)
    unq2 = len(cat2) - len(cat2mids)

    if nummatches > 0:
        newcat1 = cat1[cat1['id'].isin(cat1mids)].reset_index(drop=True)
        newcat2 = cat2[cat2['id'].isin(cat2mids)].reset_index(drop=True)

        mintime = min(newcat1['time'].min(), newcat2['time'].min())
        mintime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(mintime))
        maxtime = max(newcat1['time'].max(), newcat2['time'].max())
        maxtime = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(maxtime))

        lines.append('Overlapping time period: %s to %s\n\n'
                     % (mintime, maxtime))

        lines.append('-- Matching Criteria --\n')
        lines.append('Time window: %s s\n' % otwindow)
        lines.append('Distance window: %s km\n\n' % distwindow)

        lines.append('-- Matching Results --\n')
        lines.append('Number of associated events: %s\n' % nummatches)
        lines.append('Number of unassociated %s events: %s\n'
                     % (cat1name, unq1))
        lines.append('Number of unassociated %s events: %s\n\n'
                     % (cat2name, unq2))

        lines.append('Minimum matched latitude: %s\n' %
                     min(newcat1['latitude'].min(), newcat2['latitude'].min()))
        lines.append('Maximum matched latitude: %s\n' %
                     max(newcat1['latitude'].max(), newcat2['latitude'].max()))
        lines.append('Minimum matched longitude: %s\n' %
                     min(newcat1['longitude'].min(),
                         newcat2['longitude'].max()))
        lines.append('Maximum matched longitude: %s\n\n' %
                     max(newcat1['longitude'].max(),
                         newcat2['longitude'].max()))

        lines.append('Minimum matched depth: %s\n' %
                     min(newcat1['depth'].min(), newcat2['depth'].min()))
        lines.append('Maximum matched depth: %s\n\n' %
                     max(newcat1['depth'].max(), newcat2['depth'].max()))

        lines.append('Minimum matched magnitude: %s\n' %
                     min(newcat1['mag'].min(), newcat2['mag'].min()))
        lines.append('Maximum matched magnitude: %s' %
                     max(newcat1['mag'].max(), newcat2['mag'].max()))

    else:
        lines.append('-- Matching Criteria --\n')
        lines.append('Time window: %s s\n' % otwindow)
        lines.append('Distance window: %s km\n\n' % distwindow)

        lines.append('-- Matching Results --\n')
        lines.append('NO MATCHING EVENTS FOUND')

    with open('%s_comparisoncriteria.txt' % dirname, 'w') as compfile:
        for line in lines:
            compfile.write(line)


@printstatus('Matching events')
def match_events(cat1, cat2, dirname, otwindow=16, distwindow=100):
    """Match events within two catalogs."""
    cat1ids, cat2ids = [], []
    matchlines = [('Matching events using %ss time threshold and %skm '
                   'distance threshold\n') % (otwindow, distwindow),
                  '***********************\n']
    pcolumns = ['convtime', 'id', 'latitude', 'longitude', 'depth', 'mag']
    sep = '-----------------------\n'

    for i in range(len(cat1)):
        cat2ix = cat2[cat2['time'].between(cat1.ix[i]['time'] - otwindow,
            cat1.ix[i]['time'] + otwindow)].index.values

        if len(cat2ix) != 0:
            dists = np.array([gps2dist_azimuth(cat1.ix[i]['latitude'],
                cat1.ix[i]['longitude'], cat2.ix[x]['latitude'],
                cat2.ix[x]['longitude'])[0] / 1000. for x in cat2ix])
            dtimes = np.array([abs(cat1.ix[i]['time'] - cat2.ix[x]['time'])
                              for x in cat2ix])
            carr = dists + 5*dtimes

            ind = np.argmin(carr)

            if (dists[ind] < distwindow) and (dtimes[ind] < otwindow):
                cat1event = cat1.ix[i][pcolumns]
                cat2event = cat2.ix[cat2ix[ind]][pcolumns]
                dmag = cat1event['mag'] - cat2event['mag']
                diffs = map('{:.2f}'.format, [dists[ind], dtimes[ind], dmag])

                mline1 = ' '.join([str(x) for x in cat1event[:]]) + ' ' +\
                         ' '.join(diffs) + '\n'
                mline2 = ' '.join([str(x) for x in cat2event[:]]) + '\n'
                matchlines.extend((sep, mline1, mline2))

                cat1ids.append(cat1event['id'])
                cat2ids.append(cat2event['id'])

    cat1matched = cat1[cat1['id'].isin(cat1ids)].reset_index(drop=True)
    cat2matched = cat2[cat2['id'].isin(cat2ids)].reset_index(drop=True)

    with open('%s_matches.txt' % dirname, 'w') as matchfile:
        for mline in matchlines:
            matchfile.write(mline)

    return cat1ids, cat2ids, cat1matched, cat2matched


@printstatus('Finding closest unassociated events')
def find_closest(cat1, cat1name, cat1mids, cat2, dirname):
    """Find closest event for unassociated events."""
    cat1un = cat1[~cat1['id'].isin(cat1mids)].reset_index(drop=True)

    clines = ['Closest unassociated events for %s events\n' % cat1name,
               '***********************\n'
               'date time id latitude longitude depth magnitude '
               '(distance) (Δ time) (Δ magnitude)\n']
    pcolumns = ['convtime', 'id', 'latitude', 'longitude', 'depth', 'mag']
    sep = '-----------------------\n'

    for i in range(len(cat1un)):
        cat2ix = cat2[cat2['time'].between(cat1un.ix[i]['time'],
            cat1.ix[i]['time'] + 300)].index.values
        x = 600
        while len(cat2ix) == 0:
            cat2ix = cat2[cat2['time'].between(cat1un.ix[i]['time'],
                cat1.ix[i]['time'] + x)].index.values
            x += 6000

        dists = np.array([gps2dist_azimuth(cat1un.ix[i]['latitude'],
            cat1un.ix[i]['longitude'], cat2.ix[x]['latitude'],
            cat2.ix[x]['longitude'])[0] / 1000. for x in cat2ix])
        dtimes = np.array([abs(cat1un.ix[i]['time'] - cat2.ix[x]['time'])
                          for x in cat2ix])
        carr = dists + 5*dtimes

        ind = np.argmin(carr)

        cat1event = cat1un.ix[i][pcolumns]
        cat2event = cat2.ix[cat2ix[ind]][pcolumns]
        dmag = cat1event['mag'] - cat2event['mag']
        diffs = map('{:.2f}'.format, [dists[ind], dtimes[ind], dmag])

        cline1 = ' '.join([str(x) for x in cat1event[:]]) + '\n'
        cline2 = ' '.join([str(x) for x in cat2event[:]]) + ' ' +\
                 ' '.join(diffs) + '\n'
        clines.extend((sep, cline1, cline2))

    with open('%s_closestunassociated.txt' % dirname, 'w') as unfile:
        for cline in clines:
            unfile.write(cline)


@printstatus('Mapping events from both catalogs')
def map_events(cat1, cat1name, cat2, cat2name, cat1mids, cat2mids, dirname):
    """Map matching events between catalogs."""
    if len(cat1mids) == 0:
        return

    lllat, lllon, urlat, urlon, _, _, _, clon = qcu.get_map_bounds(cat1, cat2)

    cat1lons, cat1lats, cat2lons, cat2lats = [], [], [], []
    for i, mid in enumerate(cat1mids):
        cat1lons.append(cat1[cat1['id'] == mid]['longitude'].get_values()[0])
        cat1lats.append(cat1[cat1['id'] == mid]['latitude'].get_values()[0])
        cat2lons.append(cat2[cat2['id'] == cat2mids[i]]['longitude'
                            ].get_values()[0])
        cat2lats.append(cat2[cat2['id'] == cat2mids[i]]['latitude'
                            ].get_values()[0])

    plt.figure(figsize=(12, 7))
    mplmap = plt.axes(projection=ccrs.PlateCarree(central_longitude=clon))
    mplmap.set_extent([lllon, urlon, lllat, urlat], ccrs.PlateCarree())
    mplmap.coastlines('50m')
    mplmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                     linewidth=1, color='gray', alpha=0.5, linestyle='--')

    for i, lat in enumerate(cat1lats):
        mplmap.plot([cat1lons[i], cat2lons[i]], [lat, cat2lats[i]],
                    color='k', transform=ccrs.PlateCarree())

    mplmap.scatter(cat1lons, cat1lats, color='b', s=2, zorder=4,
                   transform=ccrs.PlateCarree(), label=cat1name)
    mplmap.scatter(cat2lons, cat2lats, color='r', s=2, zorder=4,
                   transform=ccrs.PlateCarree(), label=cat2name)
    mplmap.add_feature(cfeature.NaturalEarthFeature('cultural',
        'admin_1_states_provinces_lines', '50m', facecolor='none',
        edgecolor='k', zorder=9))
    mplmap.add_feature(cfeature.BORDERS)
    plt.legend()

    plt.savefig('%s_mapmatcheddetecs.png' % dirname, dpi=300)
    plt.close()


@printstatus('Mapping unassociated events')
def map_unique_events(cat, catname, mids):
    """Map unassociated events from a catalog."""
    if len(mids) == len(cat):
        return

    cat = cat[~cat['id'].isin(mids)].reset_index(drop=True)
    lllat, lllon, urlat, urlon, _, _, _, clon = qcu.get_map_bounds(cat)

    plt.figure(figsize=(12, 7))
    mplmap = plt.axes(projection=ccrs.PlateCarree(central_longitude=clon))
    mplmap.coastlines('50m')
    mplmap.scatter(cat['longitude'].tolist(), cat['latitude'].tolist(),
                   color='r', s=2, zorder=4, transform=ccrs.PlateCarree())
    mplmap.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                     linewidth=1, color='gray', alpha=0.5, linestyle='--')
    mplmap.add_feature(cfeature.NaturalEarthFeature('cultural',
        'admin_1_states_provinces_lines', '50m', facecolor='none',
        edgecolor='k', zorder=9))
    mplmap.add_feature(cfeature.BORDERS)
    plt.title('%s unassociated events' % catname, fontsize=20, y=1.08)

    plt.savefig('%s_uniquedetecs.png' % catname, dpi=300)
    plt.close()


@printstatus('Graphing polar histogram of azimuths and distances')
def make_az_dist(cat1, cat1name, cat2, cat2name, cat1mids, cat2mids, dirname,
                 distwindow=100, numbins=16):
    """Make polar scatter/histogram of azimuth vs. distance."""
    if len(cat1mids) == 0:
        return

    azimuths, distances = qcu.get_azs_and_dists(cat1, cat2, cat1mids, cat2mids)

    width = 2*pi / numbins
    razimuths = list(map(radians, azimuths))
    bins = np.linspace(0, 2*pi, numbins+1)
    azhist = np.histogram(razimuths, bins=bins)[0]
    hist = (float(distwindow)/max(azhist)) * azhist
    bins = (bins + width/2)[:-1]

    plt.figure(figsize=(6, 6))
    ax1 = plt.subplot(111, projection='polar')
    ax1.scatter(razimuths, distances, color='b', s=10)
    bars = ax1.bar(bins, hist, width=width)
    ax1.set_theta_zero_location('N')
    ax1.set_rmax(distwindow)
    ax1.set_theta_direction(-1)
    ax1.set_rlabel_position(112.5)
    ax1.set_title('%s location relative to %s' % (cat1name, cat2name),
                  fontsize=20)

    for _, hbar in list(zip(hist, bars)):
        hbar.set_facecolor('b')
        hbar.set_alpha(0.2)

    plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.11)

    plt.savefig('%s_polarazimuth.png' % dirname, dpi=300)
    plt.close()


@printstatus('Comparing parameters of matched events')
def compare_params(cat1, cat1name, cat2, cat2name, cat1mids, cat2mids, param,
                   dirname):
    """Compare parameters of matched events."""
    if len(cat1mids) == 0:
        return

    cat1params, cat2params = [], []

    for ix, eid in enumerate(cat1mids):
        param1 = float(cat1[cat1['id'] == eid][param])
        param2 = float(cat2[cat2['id'] == cat2mids[ix]][param].get_values()[0])

        cat1params.append(param1)
        cat2params.append(param2)

    minparam = min(min(cat1params), min(cat2params))
    maxparam = max(max(cat1params), max(cat2params))
    xes = range(int(minparam), ceil(maxparam))

    mval, bval, rval, _, _ = stats.linregress(cat1params, cat2params)
    linegraph = [mval*x + bval for x in xes]
    r2val = rval*rval

    aparam = param if param != 'mag' else 'magnitude'
    tparam = aparam.capitalize()

    plt.figure(figsize=(8, 8))
    plt.scatter(cat1params, cat2params, edgecolor='b', facecolor=None)
    plt.plot(xes, linegraph, c='r', linewidth=1, label='best fit')
    plt.plot(xes, xes, c='k', linewidth=1, label='m = 1')
    plt.legend(loc='upper left')
    plt.xlim(minparam, maxparam)
    plt.ylim(minparam, maxparam)
    plt.xlabel('%s %s' % (cat1name, aparam), fontsize=14)
    plt.ylabel('%s %s' % (cat2name, aparam), fontsize=14)
    plt.axes().set_aspect('equal', 'box')

    plt.title('%s correlation' % tparam, fontsize=20)

    plt.savefig('%s_compare%s.png' % (dirname, param), dpi=300)
    plt.close()


@printstatus('Graphing parameter differences between matched events')
def make_diff_hist(cat1, cat2, cat1mids, cat2mids, param, binsize, dirname,
                   title='', xlabel=''):
    """Make histogram of parameter differences between matched detections."""
    if len(cat1mids) == 0:
        return

    paramdiffs = []

    for idx, eid in enumerate(cat1mids):
        c1mask = cat1['id'] == eid
        c2mask = cat2['id'] == cat2mids[idx]

        if param == 'distance':
            cat1lat = cat1[c1mask]['latitude'].values[0]
            cat1lon = cat1[c1mask]['longitude'].values[0]
            cat2lat = cat2[c2mask]['latitude'].values[0]
            cat2lon = cat2[c2mask]['longitude'].values[0]
            pardiff = gps2dist_azimuth(cat1lat, cat1lon, cat2lat, cat2lon
                                      )[0] / 1000.
            paramdiffs.append(pardiff)
        else:
            cat1param = cat1[c1mask][param].values[0]
            cat2param = cat2[c2mask][param].values[0]
            if np.isnan(cat1param) or np.isnan(cat2param):
                continue
            pardiff = cat1param - cat2param
            paramdiffs.append(pardiff)

    minpardiff, maxpardiff = min(paramdiffs), max(paramdiffs)
    pardiffdown = qcu.round2bin(minpardiff, binsize, 'down')
    pardiffup = qcu.round2bin(maxpardiff, binsize, 'up')
    numbins = int((pardiffup-pardiffdown) / binsize)
    pardiffbins = np.linspace(pardiffdown, pardiffup+binsize,
                              numbins+2) - binsize/2.

    plt.figure(figsize=(12, 6))
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel('Count', fontsize=14)

    plt.hist(paramdiffs, pardiffbins, alpha=1, color='b', edgecolor='k')

    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.11)
    plt.tick_params(labelsize=12)
    plt.xlim(pardiffdown-binsize/2., pardiffup+binsize/2.)
    plt.ylim(0)

    plt.savefig('%s_%sdiffs.png' % (dirname, param), dpi=300)
    plt.close()


###############################################################################
###############################################################################
###############################################################################


def create_figures():
    """Generate and save all relevant figures and text files."""
    parser = argparse.ArgumentParser()

    parser.add_argument('catalog1', nargs='?', type=str,
                        help='pick first catalog to download data from; to \
                        download data from all catalogs, use "preferred"; if \
                        using -sf, give catalog name')
    parser.add_argument('catalog2', nargs='?', type=str,
                        help='pick second catalog to download data from; to \
                        download data from all catalogs, use "preferred"; if \
                        using -sf, give catalog name')
    parser.add_argument('startyear', nargs='?', type=int,
                        help='pick starting year; if using -sf, give first \
                        year in catalog')
    parser.add_argument('endyear', nargs='?', type=int,
                        help='pick end year (to get a single year of data, \
                        enter same year as startyear); if using -sf, give \
                        last year in catalog')

    parser.add_argument('-mr', '--magrange', nargs=2, type=float,
                        default=[-5, 12],
                        help='give the magnitude range for downloading data \
                        (default range is from -5 to 12)')
    parser.add_argument('-sf', '--specifyfiles', nargs=2, type=str,
                        help='specify two existing .csv files to use')
    parser.add_argument('-fd', '--forcedownload', action='store_true',
                        help='forces downloading of data even if .csv file \
                        exists')
    parser.add_argument('-nm', '--nomatches', action='store_false',
                        help='do not include list of matching events in HTML \
                        report')

    args = parser.parse_args()

    minmag, maxmag = args.magrange

    if args.specifyfiles is None:

        if not args.catalog1:
            sys.stdout.write('No first catalog specified. Exiting...\n')
            sys.exit()
        elif not args.catalog2:
            sys.stdout.write('No second catalog specified. Exiting...\n')
            sys.exit()
        elif not args.startyear:
            sys.stdout.write('No starting year specified. Exiting...\n')
            sys.exit()
        elif not args.endyear:
            sys.stdout.write('No ending year specified. Exiting...\n')
            sys.exit()

        cat1, cat2 = args.catalog1.lower(), args.catalog2.lower()
        startyear, endyear = map(int, [args.startyear, args.endyear])
        download = args.forcedownload

        dirname = '%s-%s%s-%s' % (cat1, cat2, startyear, endyear)

        if download:
            try:
                os.makedirs(dirname)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            datadf1 = qcu.get_data(cat1, dirname, startyear=startyear,
                endyear=endyear, minmag=minmag, maxmag=maxmag)
            datadf2 = qcu.get_data(cat2, dirname, startyear=startyear,
                endyear=endyear, minmag=minmag, maxmag=maxmag)
        else:
            # Python 2
            try:
                try:
                    datadf1 = pd.read_csv('%s/%s%s-%s.csv' %
                                          (dirname, cat1, startyear, endyear))
                    datadf2 = pd.read_csv('%s/%s%s-%s.csv' %
                                          (dirname, cat2, startyear, endyear))
                except IOError:
                    try:
                        os.makedirs(dirname)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                    datadf1 = qcu.get_data(cat1, dirname, startyear=startyear,
                        endyear=endyear, minmag=minmag, maxmag=maxmag)
                    datadf2 = qcu.get_data(cat2, dirname, startyear=startyear,
                        endyear=endyear, minmag=minmag, maxmag=maxmag)
            # Python 3
            except:
                try:
                    datadf1 = pd.read_csv('%s/%s%s-%s.csv' %
                                          (dirname, cat1, startyear, endyear))
                    datadf2 = pd.read_csv('%s/%s%s-%s.csv' %
                                          (dirname, cat2, startyear, endyear))
                except FileNotFoundError:
                    try:
                        os.makedirs(dirname)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                    datadf1 = qcu.get_data(cat1, dirname, startyear=startyear,
                        endyear=endyear, minmag=minmag, maxmag=maxmag)
                    datadf2 = qcu.get_data(cat2, dirname, startyear=startyear,
                        endyear=endyear, minmag=minmag, maxmag=maxmag)

    else:
        from shutil import copy2

        sfcat1, sfcat2 = args.specifyfiles
        cat1, cat2 = args.catalog1, args.catalog2
        dirname = '%s-%s%s-%s' % (cat1, cat2, args.startyear, args.endyear)

        try:
            os.makedirs(dirname)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        datadf1, datadf2 = pd.read_csv(sfcat1), pd.read_csv(sfcat2)
        try:
            copy2(sfcat1, dirname)
            copy2(sfcat2, dirname)
        except shutil.SameFileError:
            pass

    if len(datadf1) == 0:
        sys.stdout.write(('%s catalog has no data available for that time '
                          'period. Quitting...\n') % cat1.upper())
        sys.exit()

    if len(datadf2) == 0:
        sys.stdout.write(('%s catalog has no data available for that time '
                          'period. Quitting...\n') % cat2.upper())
        sys.exit()

    cat1, cat2 = cat1.upper(), cat2.upper()

    os.chdir(dirname)
    basic_cat_sum(datadf1, cat1, dirname)
    basic_cat_sum(datadf2, cat2, dirname)

    datadf1.loc[:, 'convtime'] = [' '.join(x.split('T')) for x in
                               datadf1['time'].tolist()]
    datadf1.loc[:, 'convtime'] = datadf1['convtime'].astype('datetime64[ns]')
    datadf1.loc[:, 'time'] = [qcu.to_epoch(x) for x in datadf1['time']]
    datadf2.loc[:, 'convtime'] = [' '.join(x.split('T')) for x in
                               datadf2['time'].tolist()]
    datadf2.loc[:, 'convtime'] = datadf2['convtime'].astype('datetime64[ns]')
    datadf2.loc[:, 'time'] = [qcu.to_epoch(x) for x in datadf2['time']]
    datadf1, datadf2 = qcu.trim_times(datadf1, datadf2)

    cat1ids, cat2ids, newcat1, newcat2 = match_events(datadf1, datadf2,
        dirname)

    if len(cat1ids) == 0:
        sys.stdout.write('*** No matching events found ***\n')

    comp_criteria(datadf1, cat1, cat1ids, datadf2, cat2, cat2ids, dirname)
    #find_closest(datadf1, cat1, cat1ids, datadf2, dirname)

    map_events(newcat1, cat1, newcat2, cat2, cat1ids, cat2ids, dirname)
    map_unique_events(datadf1, cat1, cat1ids)
    map_unique_events(datadf2, cat2, cat2ids)
    make_az_dist(newcat1, cat1, newcat2, cat2, cat1ids, cat2ids, dirname)
    compare_params(newcat1, cat1, newcat2, cat2, cat1ids, cat2ids, 'mag',
                   dirname)
    compare_params(newcat1, cat1, newcat2, cat2, cat1ids, cat2ids, 'depth',
                   dirname)
    make_diff_hist(newcat1, newcat2, cat1ids, cat2ids, 'time', 0.5, dirname,
                   xlabel='%s-%s time differences (sec)' % (cat1.upper(),
                   cat2.upper()))
    make_diff_hist(newcat1, newcat2, cat1ids, cat2ids, 'mag', 0.1, dirname,
                   xlabel='%s-%s magnitude differences' % (cat1.upper(),
                   cat2.upper()))
    make_diff_hist(newcat1, newcat2, cat1ids, cat2ids, 'depth', 2, dirname,
                   xlabel='%s-%s depth differences (km)' % (cat1.upper(),
                   cat2.upper()))
    make_diff_hist(newcat1, newcat2, cat1ids, cat2ids, 'distance', 2, dirname,
                   xlabel='%s-%s distances (km)' % (cat1.upper(),
                   cat2.upper()))

    return dirname, args.nomatches


def generate_html(dirname, matchbool):
    """Generate an HTML file containing all of the generated images and test
    files."""
    catalog1 = dirname.split('-')[0].upper()
    catalog2 = dirname.split('-')[1][:-4].upper()
    startyear = dirname.split('-')[1][-4:]
    endyear = dirname.split('-')[-1]

    with open('{0}_summary.txt'.format(catalog1)) as sum1file:
        cat1sum = '\t\t' + '\t\t'.join(sum1file.readlines())
    with open('{0}_summary.txt'.format(catalog2)) as sum2file:
        cat2sum = '\t\t' + '\t\t'.join(sum2file.readlines())
    with open('{0}_comparisoncriteria.txt'.format(dirname)) as compfile:
        compcrit = '\t\t' + '\t\t'.join(compfile.readlines())
    with open('{0}_matches.txt'.format(dirname)) as matchfile:
        matches = '\t\t' + '\t\t'.join(matchfile.readlines())

    if matchbool:
        tocm = '- [Summary of Matching Events](#matches)\n'
        strm = ('### Summary of Matching Events <a name="matches"></a>\n---\n'
                '{0}\n').format(matches)
    else:
        tocm, strm = '', ''

    toc = ('## Contents\n'
           '- [Basic Catalog Statistics](#catstats)\n'
           '    - [Catalog 1](#cat1stats)\n'
           '    - [Catalog 2](#cat2stats)\n'
           '- [Comparison Criteria](#compcrit)\n'
           '- [Matching Event Hypocenters](#matchh)\n'
           '- [Matching Event Magnitudes](#matchm)\n'
           '- [Unassociated Events](#missevs)\n'
           '{0}---\n').format(tocm)

    mdstring = ('# Report for {1} and {2} from {3} to {4}\n'
                '## Basic Catalog Statistics <a name="catstats"></a>\n'
                '### Catalog 1 <a name="cat1stats"></a>\n---\n'
                '{5}\n'
                '### Catalog 2 <a name="cat2stats"></a>\n---\n'
                '{6}\n'
                '### Comparison Criteria <a name="compcrit"></a>\n---\n'
                '{7}\n'
                '### Matching Event Hypocenters <a name="matchh"></a>\n---\n'
                '<img width="80%" src="{0}_mapmatcheddetecs.png">\n'
                '<img width="60%" src="{0}_polarazimuth.png">\n'
                '<img width="80%" src="{0}_distancediffs.png">\n'
                '<img width="80%" src="{0}_depthdiffs.png">\n'
                '<img width="80%" src="{0}_comparedepth.png">\n'
                '<img width="80%" src="{0}_timediffs.png">\n'
                '### Matching Event Magnitudes <a name="matchm"></a>\n---\n'
                '<img width="80%" src="{0}_magdiffs.png">\n'
                '<img width="80%" src="{0}_comparemag.png">\n'
                '### Unassociated Events <a name="missevs"></a>\n---\n'
                '<img width="80%" src="{1}_uniquedetecs.png">\n'
                '<img width="80%" src="{2}_uniquedetecs.png">\n'
                '{8}'
                ).format(dirname, catalog1, catalog2, startyear, endyear,
                         cat1sum, cat2sum, compcrit, strm)

    html = markdown.markdown(toc + mdstring)

    with open('{0}_report.html'.format(dirname), 'w') as htmlfile:
        htmlfile.write(html)


if __name__ == '__main__':

    try:
        dirname, matches = create_figures()
        generate_html(dirname, matches)
    except (KeyboardInterrupt, SystemError):
        sys.stdout.write('\nProgram canceled. Exiting...\n')
        sys.exit()
