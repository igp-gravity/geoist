#!/usr/bin/env python
"""Script for generating figures of catalog statistics. Run `QCreport.py -h`
for command line usage.
"""
import os
import sys
import errno
import argparse
from datetime import date, datetime
from math import sqrt, radians, cos
import shlex

import markdown
import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['Microsoft YaHei']  # 中文字体设置
matplotlib.rcParams['axes.unicode_minus'] = False
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from obspy.geodetics.base import gps2dist_azimuth

# Python 2
try:
    from urllib2 import urlopen, HTTPError
# Python 3
except ImportError:
    from urllib.request import urlopen, HTTPError

import geoist.catalog.QCutils as qcu
from geoist.catalog.decorators import retry, printstatus


###############################################################################
import pathlib 

china_border_file='//data//CutryL.shp'
#p2 = Path.cwd()+china_border_file
#china_border_file=os.path.normpath(china_border_file).replace(os.path.sep,"/")
china_border_file = str(pathlib.Path(__file__).parent/'data'/'CutryL.shp')
fault_file = str(pathlib.Path(__file__).parent/'data'/'gem_active_faults.shp')
##
network='us'
start_year='2010'
end_year='2012'
time_window=2.0
dist_window=15.0
min_mag = -5
max_mag = 12
#pathname='/home/zhangb/work/people/wangxiaoli/20180812/neic-catalog-qc-master/tmp'
#pathname=os.path.normpath(pathname).replace(os.path.sep,"/")
pathname = str(pathlib.Path(__file__).parent/'data')
catalog_file='us2010-2012.csv'
##
to_show={}
to_show['mulujichuxinxi']="true"
to_show['dizhenhuodong']="ture"
to_show['zhenyuanshendufenbu']="true"
to_show['zhenjitongji']="true"
to_show['dizhenpincitongji']="true"
to_show['chongfushijian']="true"
to_show['leijidizhenjushifang']="true"
to_show['shaixuanmulu']="true"
to_show['dizhenleixingtongji']="true"
for key in to_show:
    to_show[key]=to_show[key].lower()
###############################################################################

@printstatus('Creating basic catalog summary')
def basic_cat_sum(catalog, prefix, dup1, dup2, timewindow, distwindow):
    """Gather basic catalog summary statistics."""
    lines = []

    lines.append('目录: %s\n\n' % prefix[:-9].upper())

    lines.append('目录首次地震: %s\n' % catalog['time'].min())
    lines.append('目录结束地震: %s\n\n' % catalog['time'].max())

    lines.append('目录事件总数: %s\n\n' % len(catalog))

    lines.append('目录最小纬度: %s\n' % catalog['latitude'].min())
    lines.append('目录最大纬度: %s\n' % catalog['latitude'].max())
    lines.append('目录最小经度: %s\n' % catalog['longitude'].min())
    lines.append('目录最大经度: %s\n\n' % catalog['longitude'].max())

    lines.append('目录最小深度: %s\n' % catalog['depth'].min())
    lines.append('目录最大深度: %s\n' % catalog['depth'].max())
    lines.append('目录0km深度地震数目: %s\n'
                 % len(catalog[catalog['depth'] == 0]))
    lines.append('目录深度为空地震数目: %s\n\n'
                 % len(catalog[pd.isnull(catalog['depth'])]))

    lines.append('目录最小震级: %s\n' % catalog['mag'].min())
    lines.append('目录最大震级: %s\n' % catalog['mag'].max())
    lines.append('目录震级为0地震数据: %s\n'
                 % len(catalog[catalog['mag'] == 0]))
    lines.append('目录震级为空地震数据: %s\n\n'
                 % len(catalog[pd.isnull(catalog['mag'])]))

    lines.append('可能重复地震事件数目 (%ss 和 %skm 阈值): %d\n'
                 % (timewindow, distwindow, dup1))
    lines.append('可能重复地震事件数目 (16s 和 100km 阈值): %d'
                 % dup2)

    with open('%s_catalogsummary.txt' % prefix, 'w') as sumfile:
        for line in lines:
            sumfile.write(line)


def largest_ten(catalog, prefix):
    """Make a list of the 10 events with largest magnitude."""
    catalog = catalog.sort_values(by='mag', ascending=False)
    topten = catalog.head(n=10)
    topten = topten[['time', 'id', 'latitude', 'longitude', 'depth', 'mag']]

    with open('%s_largestten.txt' % prefix, 'w') as magfile:
        for event in topten.itertuples():
            line = ' '.join([str(x) for x in event[1:]]) + '\n'
            magfile.write(line)

@printstatus('Finding possible duplicates')
def my_list_duplicates(catalog, prefix, newcsv="newcat.csv",timewindow=2, distwindow=15,
                    magwindow=None, minmag=-5, locfilter=None):
    """Make a list of possible duplicate events."""
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T'))
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')
    catalog = catalog[catalog['mag'] >= minmag]
    if locfilter:
        catalog = catalog[catalog['place'].str.contains(locfilter, na=False)]
    cat = catalog[['time', 'convtime', 'id', 'latitude', 'longitude', 'depth',
                   'mag']].copy()
    cat.loc[:, 'time'] = [qcu.to_epoch(x) for x in cat['time']]
    n_cat=len(catalog.index)

    duplines1 = [('可能的重复地震（时间窗口：%ss，距离窗口：%skm ）\n') % (timewindow, distwindow),
                 '***********************\n'
                 '时间 纬度 经度 深度 震级'
                 ' (距离) (时间差) (震级差)\n']
    duplines2 = [('可能的重复地震（时间窗口：%ss，距离窗口：%skm ）\n') % (16, 100),
                 '***********************\n'
                 '时间 纬度 经度 深度 震级'
                 ' (距离) (时间差) (震级差)\n']
    sep = '-----------------------\n'
    # 将重复的地震标记为1丛，dupe_status记录该地震属于第几丛，第0丛表示
    # 没有重复地震
    dupe_status1=[0]*n_cat
    dupe_status2=[0]*n_cat
    # 下面两个变量记录当前记录到了第几丛
    thresh1dupes, thresh2dupes = 0, 0
    # 使用移动窗口判断重复地震，start_ind是窗口起始的index
    # 只需要判断新加进来的地震和窗口内的每一个地震是否重复
    # 窗口从第0个地震开始滑动
    start_ind=0
    for i in range(1,n_cat):
        # 只需要将新加进来的地震和（窗口内,新加地震之前）的地震进行比较
        for j in range(start_ind,i):
            dtime=(cat['convtime'][i]-cat['convtime'][j]).total_seconds()
            # 新加进来的地震和窗口内地震时间差超过16s，则新加地震之后的所有地震
            # 都会和该地震差超过16s，无需再将该地震和后续地震比较，因此将它移出
            # 窗口
            if(dtime>16):
                start_ind+=1
                continue
            dist=gps2dist_azimuth(cat['latitude'][i], cat['longitude'][i],
                            cat['latitude'][j], cat['longitude'][j])[0] / 1000.
            if (dist >= 100):
                continue
            # 新发现了一个地震丛
            if(dupe_status2[j]==0):
                thresh2dupes+=1
                dupe_status2[j]=thresh2dupes
            # 将新加入的地震和它的重复地震标为同一丛
            dupe_status2[i]=dupe_status2[j]
            dmag = cat['mag'][i] - cat['mag'][j]
#            diffs = map('{:.2f}'.format, [dist, dtime, dmag])
#            dupline1 = ' '.join([str(x) for x in cat.iloc[j,1:]]) + ' ' +\
#                ' '.join(diffs) + '\n'
#            dupline2 = ' '.join([str(x) for x in cat.iloc[i,1:]]) + '\n'
#            duplines2.extend((sep, dupline1, dupline2))
            if (dist < distwindow) and (abs(dtime) < timewindow):
                if(dupe_status1[j]==0):
                    thresh1dupes+=1
                    dupe_status1[j]=thresh1dupes
                dupe_status1[i]=dupe_status1[j]
#                duplines1.extend((sep, dupline1, dupline2))
    catalog['dup1']=dupe_status1
    catalog['dup2']=dupe_status2
    # generate duplines1
    mydup1=catalog[catalog['dup1']>0].sort_values(['dup1','time'])
    for i in range(1,thresh1dupes+1):
        tmpdup=mydup1[mydup1['dup1']==i]
        duplines1.extend(sep)
        for j,row in tmpdup.iterrows(): 
            dist=gps2dist_azimuth(row['latitude'], row['longitude'],
                    tmpdup['latitude'].iloc[0], tmpdup['longitude'].iloc[0])[0] / 1000.
            dtime=(row['convtime']-tmpdup['convtime'].iloc[0]).total_seconds()
            dmag = row['mag'] - tmpdup['mag'].iloc[0]
            diffs = map('{:.2f}'.format, [dist, dtime, dmag])
            #print(row)
            #print(type(row.iloc[4]),type(row.iloc[2]),type(row.iloc[1]),type(row.iloc[0]))
            catstr= '{} {:.2f} {:.2f} {:.2f} {:.2f}'.format(row.iloc[0],
                                                                row.iloc[1],
                                                                row.iloc[2],
                                                                row.iloc[3],
                                                                row.iloc[4])
            dupline1 = catstr + ' ' + ' '.join(diffs) + '\n'
            duplines1.extend(dupline1)
    # generate duplines2
    mydup2=catalog[catalog['dup2']>0].sort_values(['dup2','time'])
    #print(mydup2)
    for i in range(1,thresh2dupes+1):
        tmpdup=mydup2[mydup2['dup2']==i]
        duplines2.extend(sep)
        for j,row in tmpdup.iterrows(): 
            dist=gps2dist_azimuth(row['latitude'], row['longitude'],
                    tmpdup['latitude'].iloc[0], tmpdup['longitude'].iloc[0])[0] / 1000.
            dtime=(row['convtime']-tmpdup['convtime'].iloc[0]).total_seconds()
            dmag = row['mag'] - tmpdup['mag'].iloc[0]
            diffs = map('{:.2f}'.format, [dist, dtime, dmag])     
            #print(row)
            catstr= '{} {:.2f} {:.2f} {:.2f} {:.2f}'.format(row.iloc[0],
                                                               row.iloc[1],
                                                               row.iloc[2],
                                                               row.iloc[3],
                                                               row.iloc[4])
            dupline2 = catstr + ' ' + ' '.join(diffs) + '\n'
            duplines2.extend(dupline2)
    # generate dupline2
    catalog.set_index('time',inplace=True)
    catalog.to_csv(newcsv)
    with open('%s_duplicates.txt' % prefix, 'w') as dupfile:
        for dupline in duplines1:
            dupfile.write(dupline)
        for dupline in duplines2:
            dupfile.write(dupline)
    # delete duplicates1
    nodup=catalog[catalog['dup1']==0]
    tmpdup=catalog[catalog['dup1']>0].drop_duplicates(['dup1'])
    newcat=pd.concat([nodup,tmpdup]).sort_values('time')
    newcat.to_csv('remove_dup_{:.0f}_{:.0f}.csv'.format(timewindow,distwindow))
    # delete duplicates2
    nodup=catalog[catalog['dup2']==0]
    tmpdup=catalog[catalog['dup2']>0].drop_duplicates(['dup2'])
    newcat=pd.concat([nodup,tmpdup])
    newcat=pd.concat([nodup,tmpdup]).sort_values('time')
    newcat.to_csv('remove_dup_{}_{}.csv'.format(16,100))
    return thresh1dupes, thresh2dupes

@printstatus('Finding possible duplicates')
def list_duplicates(catalog, prefix, timewindow=2, distwindow=15,
                    magwindow=None, minmag=-5, locfilter=None):
    """Make a list of possible duplicate events."""
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T'))
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')
    catalog = catalog[catalog['mag'] >= minmag]
    if locfilter:
        catalog = catalog[catalog['place'].str.contains(locfilter, na=False)]
    cat = catalog[['time', 'convtime', 'id', 'latitude', 'longitude', 'depth',
                   'mag']].copy()
    cat.loc[:, 'time'] = [qcu.to_epoch(x) for x in cat['time']]

    duplines1 = [('Possible duplicates using %ss time threshold and %skm '
                  'distance threshold\n') % (timewindow, distwindow),
                 '***********************\n'
                 'date time id latitude longitude depth magnitude '
                 '(distance) (Δ time) (Δ magnitude)\n']
    duplines2 = [('\n\nPossible duplicates using 16s time threshold and 100km '
                  'distance threshold\n'),
                 '***********************\n'
                 'date time id latitude longitude depth magnitude '
                 '(distance) (Δ time) (Δ magnitude)\n']
    sep = '-----------------------\n'

    thresh1dupes, thresh2dupes = 0, 0
    for event in cat.itertuples():

        trimdf = cat[cat['convtime'].between(event.convtime, event.convtime
                 + pd.Timedelta(seconds=16), inclusive=False)]

        if len(trimdf) != 0:
            for tevent in trimdf.itertuples():
                dist = gps2dist_azimuth(event.latitude, event.longitude,
                            tevent.latitude, tevent.longitude)[0] / 1000.
                if dist < 100:
                    dtime = (event.convtime - tevent.convtime).total_seconds()
                    dmag = event.mag - tevent.mag
                    diffs = map('{:.2f}'.format, [dist, dtime, dmag])

                    dupline1 = ' '.join([str(x) for x in event[1:]]) + ' ' +\
                               ' '.join(diffs) + '\n'
                    dupline2 = ' '.join([str(x) for x in tevent[1:]]) + '\n'
                    duplines2.extend((sep, dupline1, dupline2))

                    thresh2dupes += 1

                    if (dist < distwindow) and (abs(dtime) < timewindow):
                        duplines1.extend((sep, dupline1, dupline2))
                        thresh1dupes += 1

            continue

    with open('%s_duplicates.txt' % prefix, 'w') as dupfile:
        for dupline in duplines1:
            dupfile.write(dupline)
        for dupline in duplines2:
            dupfile.write(dupline)

    return thresh1dupes, thresh2dupes


@printstatus('Mapping earthquake locations')
def map_detecs(catalog, prefix, minmag=-5, mindep=-50, title=''):
    """Make scatter plot of detections with magnitudes (if applicable)."""
    catalog = catalog[(catalog['mag'] >= minmag)
                      & (catalog['depth'] >= mindep)].copy()

    if len(catalog) == 0:
        print('\nCatalog contains no events deeper than %s.' % mindep)
        return

    # define map bounds
    lllat, lllon, urlat, urlon, _, _, _, clon = qcu.get_map_bounds(catalog)

    plt.figure(figsize=(12, 7))
    mplmap = plt.axes(projection=ccrs.PlateCarree(central_longitude=clon))
    #mplmap.set_extent([lllon, urlon, lllat, urlat], ccrs.PlateCarree())
    mplmap.coastlines('50m', facecolor='none')

    # if catalog has magnitude data
    if not catalog['mag'].isnull().all():
        bins = [0, 5, 6, 7, 8, 15]
        binnames = ['< 5', '5-6', '6-7', '7-8', r'$\geq$8']
        binsizes = [10, 25, 50, 100, 400]
        bincolors = ['g', 'b', 'y', 'r', 'r']
        binmarks = ['o', 'o', 'o', 'o', '*']
        catalog.loc[:, 'maggroup'] = pd.cut(catalog['mag'], bins,
                                            labels=binnames)

        for i, label in enumerate(binnames):
            mgmask = catalog['maggroup'] == label
            rcat = catalog[mgmask]
            lons, lats = list(rcat['longitude']), list(rcat['latitude'])
            if len(lons) > 0:
                mplmap.scatter(lons, lats, s=binsizes[i], marker=binmarks[i],
                               c=bincolors[i], label=binnames[i], alpha=0.8,
                               zorder=10, transform=ccrs.PlateCarree())

        plt.legend(loc='lower left', title='Magnitude')

    # if catalog does not have magnitude data
    else:
        lons, lats = list(catalog['longitude']), list(catalog['latitude'])
        mplmap.scatter(lons, lats, s=15, marker='x', c='r', zorder=10)

    shape_feature = ShapelyFeature(Reader(china_border_file).geometries(),
                                   ccrs.PlateCarree(), edgecolor='black',
                                   facecolor='none')
    shape_feature1 = ShapelyFeature(Reader(fault_file).geometries(),
                                   ccrs.PlateCarree(), edgecolor='red',
                                   facecolor='none')        
    mplmap.add_feature(shape_feature)
    mplmap.add_feature(shape_feature1)

    plt.title(title, fontsize=20)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    if mindep != -50:
        plt.savefig('%s_morethan%sdetecs.png' % (prefix, mindep), dpi=300)
    else:
        plt.savefig('%s_mapdetecs.png' % prefix, dpi=300)

    plt.close()


@printstatus('Mapping earthquake density')
def map_detec_nums(catalog, prefix, title='', numcolors=16, rmin=77, rmax=490,
                   minmag=-5, pltevents=True):
    """Map detections and a grid of detection density. rmax=510 is white,
    rmin=0 is black.
    """
    # generate bounds for map
    mask = catalog['mag'] >= minmag

    lllat, lllon, urlat, urlon, gridsize, hgridsize, _, clon = \
        qcu.get_map_bounds(catalog[mask])

    catalog = qcu.add_centers(catalog, gridsize)
    groupedlatlons, _, cmax = qcu.group_lat_lons(catalog, minmag=minmag)

    # print message if there are no detections with magnitudes above minmag
    if cmax == 0:
        print("No detections over magnitude %s" % minmag)

    # create color gradient from light red to dark red
    colors = qcu.range2rgb(rmin, rmax, numcolors)

    # put each center into its corresponding color group
    colorgroups = list(np.linspace(0, cmax, numcolors))
    groupedlatlons.loc[:, 'group'] = np.digitize(groupedlatlons['count'],
                                                 colorgroups)

    # create map
    plt.figure(figsize=(12, 7))
    mplmap = plt.axes(projection=ccrs.PlateCarree(central_longitude=clon))
#    mplmap.set_extent([lllon, urlon, lllat, urlat], ccrs.PlateCarree())
    mplmap.coastlines('50m')
    shape_feature = ShapelyFeature(Reader(china_border_file).geometries(),
                                   ccrs.PlateCarree(), edgecolor='black',
                                   facecolor='none')
   
    shape_feature1 = ShapelyFeature(Reader(fault_file).geometries(),
                                   ccrs.PlateCarree(), edgecolor='red',
                                   facecolor='none')    
    
    mplmap.add_feature(shape_feature)
    mplmap.add_feature(shape_feature1)
    
    plt.title(title, fontsize=20)
    plt.subplots_adjust(left=0.01, right=0.9, top=0.95, bottom=0.05)

    # create color map based on rmin and rmax
    cmap = LinearSegmentedColormap.from_list('CM', colors)._resample(numcolors)

    # make dummy plot for setting color bar
    colormesh = mplmap.pcolormesh(colors, colors, colors, cmap=cmap, alpha=1,
                                  vmin=0, vmax=cmax)

    # format color bar
    cbticks = [x for x in np.linspace(0, cmax, numcolors+1)]
    cbar = plt.colorbar(colormesh, ticks=cbticks)
    cbar.ax.set_yticklabels([('%.0f' % x) for x in cbticks])
    cbar.set_label('# of detections', rotation=270, labelpad=15)

    # plot rectangles with color corresponding to number of detections
    for center, _, cgroup in groupedlatlons.itertuples():
        minlat, maxlat = center[0]-hgridsize, center[0]+hgridsize
        minlon, maxlon = center[1]-hgridsize, center[1]+hgridsize
        glats = [minlat, maxlat, maxlat, minlat]
        glons = [minlon, minlon, maxlon, maxlon]

        color = colors[cgroup-1]

        qcu.draw_grid(glats, glons, color, alpha=0.8)

    # if provided, plot detection epicenters
    if pltevents and not catalog['mag'].isnull().all():
        magmask = catalog['mag'] >= minmag
        lons = list(catalog['longitude'][magmask])
        lats = list(catalog['latitude'][magmask])
        mplmap.scatter(lons, lats, c='k', s=7, marker='x', zorder=5)
    elif catalog['mag'].isnull().all():
        lons = list(catalog['longitude'])
        lats = list(catalog['latitude'])
        mplmap.scatter(lons, lats, c='k', s=7, marker='x', zorder=5)

    plt.savefig('%s_eqdensity.png' % prefix, dpi=300)
    plt.close()


@printstatus('Making histogram of given parameter')
def make_hist(catalog, param, binsize, prefix, title='', xlabel='',
              countlabel=False, maxval=None):
    """Plot histogram grouped by some parameter."""
    paramlist = catalog[pd.notnull(catalog[param])][param].tolist()
    minparam, maxparam = min(paramlist), max(paramlist)
    paramdown = qcu.round2bin(minparam, binsize, 'down')
    paramup = qcu.round2bin(maxparam, binsize, 'up')
    numbins = int((paramup-paramdown) / binsize)
    labelbuff = float(paramup-paramdown) / numbins * 0.5

    diffs = [abs(paramlist[i+1]-paramlist[i]) for i in range(len(paramlist))
             if i+1 < len(paramlist)]
    diffs = [round(x, 1) for x in diffs if x > 0]

    plt.figure(figsize=(10, 6))
    plt.title(title, fontsize=20)
    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel('Count', fontsize=14)

    if param == 'ms':
        parambins = np.linspace(paramdown, paramup, numbins+1)
        plt.xlim(paramdown, paramup)
    else:
        parambins = np.linspace(paramdown, paramup+binsize,
                                numbins+2) - binsize/2.
        plt.xlim(paramdown-binsize/2., paramup+binsize/2.)

    phist = plt.hist(paramlist, parambins, alpha=0.7, color='b', edgecolor='k')
    maxbarheight = max([phist[0][x] for x in range(numbins)] or [0])
    labely = maxbarheight / 50.

    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.11)

    if maxval:
        plt.xlim(xmax=maxval)
    plt.ylim(0, maxbarheight*1.1+0.1)

    # put count numbers above the bars if countlabel=True
    if countlabel:
        for i in range(numbins):
            plt.text(phist[1][i]+labelbuff, phist[0][i]+labely,
                     '%0.f' % phist[0][i], size=12, ha='center')

    if maxval:
        plt.savefig('%s_zoom%shistogram.png' % (prefix, param), dpi=300)
    else:
        plt.savefig('%s_%shistogram.png' % (prefix, param), dpi=300)

    plt.close()


@printstatus('Making histogram of given time duration')
def make_time_hist(catalog, timelength, prefix, title=''):
    """Make histogram either by hour of the day or by date."""
    timelist = catalog['time']

    plt.figure(figsize=(10, 6))
    plt.title(title, fontsize=20)
    plt.ylabel('Count', fontsize=14)

    if timelength == 'hour':
        lons = np.linspace(-180, 180, 25).tolist()
        hours = np.linspace(-12, 12, 25).tolist()

        tlonlist = catalog.loc[:, ['longitude', 'time']]
        tlonlist.loc[:, 'rLon'] = qcu.round2lon(tlonlist['longitude'])

        tlonlist.loc[:, 'hour'] = [int(x.split('T')[1].split(':')[0])
                                   for x in tlonlist['time']]
        tlonlist.loc[:, 'rhour'] = [x.hour + hours[lons.index(x.rLon)]
                                    for x in tlonlist.itertuples()]

        tlonlist.loc[:, 'rhour'] = [x+24 if x < 0 else x-24 if x > 23 else x
                                    for x in tlonlist['rhour']]

        hourlist = tlonlist.rhour.tolist()
        hourbins = np.linspace(-0.5, 23.5, 25)

        plt.hist(hourlist, hourbins, alpha=1, color='b', edgecolor='k')
        plt.xlabel('Hour of the Day', fontsize=14)
        plt.xlim(-0.5, 23.5)

    elif timelength == 'day':
        daylist = [x.split('T')[0] for x in timelist]
        daydf = pd.DataFrame({'date': daylist})
        daydf['date'] = daydf['date'].astype('datetime64[ns]')
        daydf = daydf.groupby([daydf['date'].dt.year.rename('year'),
                               daydf['date'].dt.month.rename('month'),
                               daydf['date'].dt.day.rename('day')]).count()

        eqdates = daydf.index.tolist()
        counts = daydf.date.tolist()

        eqdates = [date(x[0], x[1], x[2]) for x in eqdates]
        minday, maxday = min(eqdates), max(eqdates)

        plt.bar(eqdates, counts, alpha=1, color='b', width=1)
        plt.xlabel(u'日期', fontsize=14)
        plt.xlim(minday, maxday)

    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.11)

    plt.savefig('%s_%shistogram.png' % (prefix, timelength), dpi=300)
    plt.close()


@printstatus('Graphing mean time separation')
def graph_time_sep(catalog, prefix):
    """Make bar graph of mean time separation between events by date."""
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')
    catalog.loc[:, 'dt'] = catalog.convtime.diff().astype('timedelta64[ns]')
    catalog.loc[:, 'dtmin'] = catalog['dt'] / pd.Timedelta(minutes=1)

    mindate = catalog['convtime'].min()
    maxdate = catalog['convtime'].max()

    fig = plt.figure(figsize=(10, 6))
    axfull = fig.add_subplot(111)
    axfull.set_ylabel('Time separation (min)', fontsize=14, labelpad=20)
    axfull.spines['top'].set_color('none')
    axfull.spines['bottom'].set_color('none')
    axfull.spines['left'].set_color('none')
    axfull.spines['right'].set_color('none')
    axfull.tick_params(labelcolor='w', top='off', bottom='off',
                       left='off', right='off')

    if maxdate - mindate < pd.Timedelta(days=1460):
        # time separation between events
        fig.add_subplot(311)
        plt.plot(catalog['convtime'], catalog['dtmin'], alpha=1, color='b')
        plt.xlabel('Date')
        plt.title('Time separation between events')
        plt.xlim(mindate, maxdate)
        plt.ylim(0)

        # maximum monthly time separation
        fig.add_subplot(312)
        month_max = catalog.resample('1M', on='convtime').max()['dtmin']
        months = month_max.index.map(lambda x: x.strftime('%Y-%m')).tolist()
        months = [date(int(x[:4]), int(x[-2:]), 1) for x in months]
        plt.bar(months, month_max.tolist(), color='b', alpha=1, width=31,
                edgecolor='k')
        plt.xlabel('Month')
        plt.title('Maximum event separation by month')
        plt.xlim(mindate - pd.Timedelta(days=15),
                 maxdate - pd.Timedelta(days=16))

        # median monthly time separation
        fig.add_subplot(313)
        month_med = catalog.resample('1M', on='convtime').median()['dtmin']
        plt.bar(months, month_med.tolist(), color='b', alpha=1, width=31,
                edgecolor='k')
        plt.xlabel(u'月')
        plt.title('Median event separation by month')
        plt.tight_layout()
        plt.xlim(mindate - pd.Timedelta(days=15),
                 maxdate - pd.Timedelta(days=16))

    else:
        # time separation between events
        fig.add_subplot(311)
        plt.plot(catalog['convtime'], catalog['dtmin'], alpha=1, color='b')
        plt.xlabel('Date')
        plt.title('Time separation between events')
        plt.xlim(mindate, maxdate)
        plt.ylim(0)

        # maximum yearly time separation
        fig.add_subplot(312)
        year_max = catalog.resample('1Y', on='convtime').max()['dtmin']
        years = year_max.index.map(lambda x: x.strftime('%Y')).tolist()
        years = [date(int(x[:4]), 1, 1) for x in years]
        plt.bar(years, year_max.tolist(), color='b', alpha=1, width=365,
                edgecolor='k')
        plt.xlabel('Year')
        plt.title('Maximum event separation by year')
        plt.xlim(mindate - pd.Timedelta(days=183),
                 maxdate - pd.Timedelta(days=183))

        # median yearly time separation
        fig.add_subplot(313)
        year_med = catalog.resample('1Y', on='convtime').median()['dtmin']
        plt.bar(years, year_med.tolist(), color='b', alpha=1, width=365,
                edgecolor='k')
        plt.xlabel('Year')
        plt.title('Median event separation by year')
        plt.tight_layout()
        plt.xlim(mindate - pd.Timedelta(days=183),
                 maxdate - pd.Timedelta(days=183))

    plt.savefig('%s_timeseparation.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing median magnitude by time')
def med_mag(catalog, prefix):
    """Make a bar graph of median event magnitude by year."""
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')

    mindate = catalog['convtime'].min()
    maxdate = catalog['convtime'].max()

    if maxdate - mindate < pd.Timedelta(days=1460):
        month_max = catalog.resample('1M', on='convtime').max()['mag']
        months = month_max.index.map(lambda x: x.strftime('%Y-%m')).tolist()
        months = [date(int(x[:4]), int(x[-2:]), 1) for x in months]
        month_medmag = catalog.resample('1M', on='convtime').median()['mag']

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        ax.tick_params(bottom='off')
        plt.bar(months, month_medmag.tolist(), color='b', edgecolor='k',
                alpha=1, width=31)
        plt.xlabel('Month', fontsize=14)
        plt.ylabel('Magnitude', fontsize=14)
        plt.title('Monthly Median Magnitude', fontsize=20)
        plt.xlim(min(months) - pd.Timedelta(days=15),
                 max(months) + pd.Timedelta(days=15))
    else:
        year_max = catalog.resample('1Y', on='convtime').max()['mag']
        years = year_max.index.map(lambda x: x.strftime('%Y')).tolist()
        years = [date(int(x[:4]), 1, 1) for x in years]
        year_medmag = catalog.resample('1Y', on='convtime').median()['mag']

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        ax.tick_params(bottom='off')
        plt.bar(years, year_medmag.tolist(), color='b', edgecolor='k', alpha=1,
                width=365)
        plt.xlabel('Year', fontsize=14)
        plt.ylabel('Magnitude', fontsize=14)
        plt.title('Yearly Median Magnitude', fontsize=20)
        plt.xlim(min(years) - pd.Timedelta(days=183),
                 max(years) - pd.Timedelta(days=183))
    
    plt.savefig('%s_medianmag.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing magnitude completeness')
def cat_mag_comp(catalog, prefix, magbin=0.1):
    """Plot catalog magnitude completeness."""
    catalog = catalog[pd.notnull(catalog['mag'])]
    mags = np.array(catalog['mag'])
    mags = np.around(mags, 1)

    minmag, maxmag = min(min(mags), 0), max(mags)

    mag_centers = np.arange(minmag, maxmag + 2*magbin, magbin)
    cdf = np.zeros(len(mag_centers))

    for idx in range(len(cdf)):
        cdf[idx] = np.count_nonzero(
            ~np.isnan(mags[mags >= mag_centers[idx]-0.001]))

    mag_edges = np.arange(minmag - magbin/2., maxmag+magbin, magbin)
    g_r, _ = np.histogram(mags, mag_edges)
    idx = list(g_r).index(max(g_r))

    mc_est = mag_centers[idx]

    try:
        mc_est, bvalue, avalue, lval, mc_bins, std_dev = qcu.WW2000(mc_est,
            mags, magbin)
    except:
        mc_est = mc_est + 0.3
        mc_bins = np.arange(0, maxmag + magbin/2., magbin)
        bvalue = np.log10(np.exp(1))/(np.average(mags[mags >= mc_est])
                                      - (mc_est-magbin/2.))
        avalue = np.log10(len(mags[mags >= mc_est])) + bvalue*mc_est
        log_l = avalue-bvalue*mc_bins
        lval = 10.**log_l
        std_dev = bvalue/sqrt(len(mags[mags >= mc_est]))

    plt.figure(figsize=(8, 6))
    plt.scatter(mag_centers[:len(g_r)], g_r, edgecolor='r', marker='o',
                facecolor='none', label='Incremental')
    plt.scatter(mag_centers, cdf, c='k', marker='+', label='Cumulative')
    plt.axvline(mc_est, c='r', linestyle='--', label='Mc = %2.1f' % mc_est)
    plt.plot(mc_bins, lval, c='k', linestyle='--',
             label='B = %1.3f%s%1.3f' % (bvalue, u'\u00B1', std_dev))

    ax1 = plt.gca()
    ax1.set_yscale('log')
    max_count = np.amax(cdf) + 100000
    ax1.set_xlim([minmag, maxmag])
    ax1.set_ylim([1, max_count])
    plt.title('Frequency-Magnitude Distribution', fontsize=18)
    plt.xlabel('Magnitude', fontsize=14)
    plt.ylabel('Log10 Count', fontsize=14)
    plt.legend(numpoints=1)

    plt.savefig('%s_catmagcomp.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing magnitude versus time for each earthquake')
def graph_mag_time(catalog, prefix):
    """Plot magnitudes vs. origin time."""
    catalog = catalog[pd.notnull(catalog['mag'])]

    times = catalog['convtime'].copy()
    mags = catalog['mag'].copy()

    plt.figure(figsize=(10, 6))
    plt.xlabel('Date', fontsize=14)
    plt.ylabel('Magnitude', fontsize=14)
    plt.plot_date(times, mags, alpha=0.7, markersize=2, c='b')
    plt.xlim(min(times), max(times))
    plt.title('Magnitude vs. Time', fontsize=20)

    plt.savefig('%s_magvtime.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing event count by date and magnitude')
def graph_mag_count(catalog, prefix):
    """Graph event count grouped by magnitude and by date."""
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')

    mindate, maxdate = catalog['convtime'].min(), catalog['convtime'].max()
    bincond = maxdate - mindate < pd.Timedelta(days=1460)
    barwidth = 31 if bincond else 365
    timedelt = pd.Timedelta(days=barwidth/2.)

    minbin = qcu.round2bin(catalog['mag'].min()-0.1, 1, 'down')
    maxbin = qcu.round2bin(catalog['mag'].max()+0.1, 1, 'up')
    bins = np.arange(minbin, maxbin+0.1, 1)

    catalog.loc[:, 'magbin'] = pd.cut(catalog['mag'], bins=bins, right=True)
    maggroups = catalog['magbin'].sort_values().unique()

    fig, axlist = plt.subplots(len(maggroups), sharex=True)
    fig.set_size_inches(10, 14, forward=True)
    for i, mbin in enumerate(maggroups):

        trimcat = catalog[catalog['magbin'] == mbin]

        if len(trimcat) == 0:
            continue

        datelist = [x.split('T')[0] for x in trimcat['time']]
        datedf = pd.DataFrame({'date': datelist})
        datedf['date'] = datedf['date'].astype('datetime64[ns]')

        datedf = datedf.groupby([datedf['date'].dt.year.rename('year'),
                 datedf['date'].dt.month.rename('month')]).count() if bincond\
                 else datedf.groupby([datedf['date'].dt.year.rename('year')]).count()

        evdates = datedf.index.tolist()
        counts = datedf.date.tolist()

        evdates = [date(x[0], x[1], 1) if bincond else date(x, 1, 1)
                   for x in evdates]

        axlist[i].bar(evdates, counts, alpha=1, color='b', width=barwidth,
                      edgecolor='k')
        axlist[i].set_ylabel('%d-%d' % (bins[i], bins[i+1]), fontsize=10)
        plt.xlim(mindate - timedelt, maxdate - timedelt)
        plt.ylim(0, max(counts))
        axlist[i].get_yaxis().set_label_coords(-0.1, 0.5)

    plt.xlabel('Date', fontsize=14)
    for ax in axlist[:-1]:
        ax.xaxis.set_ticks_position('none')

    plt.savefig('%s_magtimecount.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing cumulative moment release')
def cumul_moment_release(catalog, prefix):
    """Graph cumulative moment release."""
    catalog = catalog[pd.notnull(catalog['mag'])]
#    catalog.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
#                                  for x in catalog['time'].tolist()]
#    catalog.loc[:, 'convtime'] = catalog['convtime'].astype('datetime64[ns]')
    times = catalog['convtime']

    minday, maxday = min(times), max(times)

    mag0 = 10.**((3/2.)*(catalog['mag']+10.7))
    mag0 = mag0 * 10.**(-7)
    cumulmag0 = np.cumsum(mag0)

    plt.figure(figsize=(10, 6))
    plt.plot(times, cumulmag0, 'k-')
    plt.xlabel('Date', fontsize=14)
    plt.ylabel(r'Cumulative Moment Release (N$\times$m)', fontsize=14)
    plt.xlim(minday, maxday)
    plt.ylim(0)
    plt.title('Cumulative Moment Release', fontsize=20)

    colors = ['r', 'm', 'c', 'y', 'g']
    largesteqs = catalog.sort_values('mag').tail(5)
    for i, eq in enumerate(largesteqs.itertuples()):
        plt.axvline(x=eq.time, color=colors[i], linestyle='--')

    plt.savefig('%s_cumulmomentrelease.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing cumulative event types')
def graph_event_types(catalog, prefix):
    """Graph number of cumulative events by type of event."""
    typedict = {}

    for evtype in catalog['type'].unique():
        typedict[evtype] = (catalog['type'] == evtype).cumsum()

    plt.figure(figsize=(12, 6))

    for evtype in typedict:
        plt.plot_date(catalog['convtime'], typedict[evtype], marker=None,
                      linestyle='-', label=evtype)

    plt.yscale('log')
    plt.legend()
    plt.xlim(min(catalog['convtime']), max(catalog['convtime']))

    plt.xlabel('Date', fontsize=14)
    plt.ylabel('Cumulative number of events', fontsize=14)
    plt.title('Cumulative Event Type', fontsize=20)

    plt.savefig('%s_cumuleventtypes.png' % prefix, dpi=300)
    plt.close()


@printstatus('Graphing possible number of duplicate events')
def cat_dup_search(catalog, prefix):
    """Graph possible number of duplicate events given various distances
    and time differences.
    """
    epochtimes = [qcu.to_epoch(row.time) for row in catalog.itertuples()]
    tdifsec = np.asarray(abs(np.diff(epochtimes)))

    lat1 = np.asarray(catalog.latitude[:-1])
    lon1 = np.asarray(catalog.longitude[:-1])
    lat2 = np.asarray(catalog.latitude[1:])
    lon2 = np.asarray(catalog.longitude[1:])
    ddelkm = [gps2dist_azimuth(lat1[i], lon1[i], lat2[i], lon2[i])[0] / 1000.
              for i in range(len(lat1))]

    diffdf = pd.DataFrame({'tdifsec': tdifsec, 'ddelkm': ddelkm})

    kmlimits = [1, 2, 4, 8, 16, 32, 64, 128, 256]
    tmax = 16
    dtime = 0.05
    timebins = np.arange(0, tmax+dtime/2, dtime)

    numevents = np.empty([len(kmlimits), len(timebins)-1])

    for jdx in range(len(kmlimits)):

        cat_subset = diffdf[diffdf.ddelkm <= kmlimits[jdx]]

        for idx in range(len(timebins)-1):

            numevents[jdx][idx] = cat_subset[cat_subset.tdifsec.between(
                timebins[idx], timebins[idx+1])].count()[0]

    totmatch = np.transpose(np.cumsum(np.transpose(numevents), axis=0))

    plt.figure(figsize=(10, 6))
    for idx in range(len(kmlimits)):
        times = timebins[1:]
        matches = totmatch[idx]
        lab = str(kmlimits[idx]) + ' km'
        plt.plot(times, matches, label=lab)

    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel('Possible duplicate events', fontsize=14)
    plt.xlim(0, tmax)
    plt.ylim(0, np.amax(totmatch)+0.5)
    plt.legend(loc=2, numpoints=1)
    plt.title(('Cumulative number of events within X seconds\n'
               'and Z km (Z specified in legend)'), fontsize=20)

    plt.savefig('%s_catdupsearch.png' % prefix, dpi=300)
    plt.close()

def create_figures_new(to_plot=None, db = None, catalog_file = None):
    #from shutil import copy2
    #prefix = '.'.join(args.specifyfile.split('.')[:-1])

    catalog = network
    startyear, endyear = map(int, [start_year, end_year])
    prefix = '%s%s-%s' % (catalog, start_year, end_year) if catalog else\
              'preferred%s-%s' % (start_year, end_year)
    #print(filename)
    #print(start_year, end_year)
    if db is None:          
        filename = pathname + "\\" + catalog_file
        datadf=qcu.get_local_data(filename, prefix, startyear=int(start_year),
                    endyear=int(end_year), minmag=min_mag, maxmag=max_mag)
    else:
        datadf=qcu.get_db_data(db, prefix, startyear=int(start_year),
                    endyear=int(end_year), minmag=min_mag, maxmag=max_mag)

    #print(datadf)
    if len(datadf) == 0:
        sys.stdout.write(('Catalog has no data available for that time period.'
                          ' Quitting...\n'))
        sys.exit()

    timewindow = time_window
    distwindow = dist_window
    
    if db is None: 
        datadf = datadf.sort_values(by='time').reset_index(drop=True)
        datadf.loc[:, 'ms'] = datadf['time'].str[-4:-1].astype('float')
        datadf.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
                                      for x in datadf['time'].tolist()]
        datadf.loc[:, 'convtime'] = datadf['convtime'].astype('datetime64[ns]')
    else:
        datadf = datadf.sort_values(by='time').reset_index(drop=True)
        datadf.loc[:, 'ms'] = datadf['time'].str[-7:-1].astype('float')/1000
        datadf.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
                                      for x in datadf['time'].tolist()]
        datadf.loc[:, 'convtime'] = datadf['convtime'].astype('datetime64[ns]')
    
    os.chdir(pathname)
    dup1, dup2 = my_list_duplicates(datadf, prefix, timewindow=timewindow,
                                 distwindow=distwindow)
    basic_cat_sum(datadf, prefix, dup1, dup2, timewindow, distwindow)
    largest_ten(datadf, prefix)

    # generate figures
    if (to_plot is None or to_plot['Detect Location']=='true'):
        map_detecs(datadf, prefix, title='Detection locations')
    if(to_plot is None or to_plot['Deeper than 50km']=='true'):
        map_detecs(datadf, prefix, mindep=50, title='Detections deeper than 50km')
    if(to_plot is None or to_plot['Detect Density']=='true'):
        map_detec_nums(datadf, prefix, title='Detection density')
    if(to_plot is None or to_plot['Mag Hist']=='true'):
        make_hist(datadf, 'mag', 0.1, prefix, xlabel='Magnitude',
                  title='Magnitude histogram')
    if(to_plot is None or to_plot['Depth Hist']=='true'):
        make_hist(datadf, 'depth', 1, prefix, xlabel='Depth (km)',
                  title='Depth histogram')
    if(to_plot is None or to_plot['Zoomed Depth Hist']=='true'):
        make_hist(datadf, 'depth', 0.5, prefix, maxval=20, xlabel='Depth (km)',
                  title='Zoomed depth histogram')
#    make_hist(datadf, 'ms', 20, prefix, xlabel='Milliseconds',
#              title='Histogram of milliseconds')
    if(to_plot is None or to_plot['Events per Hour']=='true'):
        make_time_hist(datadf, 'hour', prefix, title='Events per Hour of the Day')
    if(to_plot is None or to_plot['Events per Day']=='true'):
        make_time_hist(datadf, 'day', prefix, title='Events per Day')
    if(to_plot is None or to_plot['Time Seperation']=='true'):
        graph_time_sep(datadf, prefix)
    if(to_plot is None or to_plot['Mag Completeness']=='true'):
        cat_mag_comp(datadf, prefix)
    if(to_plot is None or to_plot['Mag Time']=='true'):
        graph_mag_time(datadf, prefix)
    if(to_plot is None or to_plot['Cumu Moment Release']=='true'):
        cumul_moment_release(datadf, prefix)
    if(to_plot is None or to_plot['Event Types']=='true'):
        graph_event_types(datadf, prefix)
    if(to_plot is None or to_plot['Num of Duplicates']=='true'):
        cat_dup_search(datadf, prefix)
    if(to_plot is None or to_plot['Med Mag']=='true'):
        med_mag(datadf, prefix)
    if(to_plot is None or to_plot['Mag Count']=='true'):
        graph_mag_count(datadf, prefix)

    return pathname,prefix


def create_figures(to_plot=None):
    """Generate and save all relevant figures and text files."""
    parser = argparse.ArgumentParser()

    parser.add_argument('catalog', nargs='?', type=str,
                        help='pick which catalog to download data from; to \
                        download data from all catalogs, use "preferred"')
    parser.add_argument('startyear', nargs='?', type=int,
                        help='pick starting year')
    parser.add_argument('endyear', nargs='?', type=int,
                        help='pick end year (to get a single year of data, \
                        enter same year as startyear)')

    parser.add_argument('-mr', '--magrange', type=float, nargs=2,
                        default=[-5, 12],
                        help='give the magnitude range for downloading data \
                        (default range is from -5 to 12)')
    parser.add_argument('-tw', '--timewindow', type=float, default=2,
                        help='change time window for finding duplicates \
                        (default is 2 seconds)')
    parser.add_argument('-dw', '--distwindow', type=float, default=15,
                        help='change distance window for finding duplicates \
                        (default is 15 kilometers)')
    parser.add_argument('-sf', '--specifyfile', type=str,
                        help='specify existing .csv file to use')
    parser.add_argument('-fd', '--forcedownload', action='store_true',
                        help='forces downloading of data even if .csv file \
                        exists')
    myargs=[]
    myargs.append(network)
    myargs.append(start_year)
    myargs.append(end_year)
    myargs.append('-sf {}/{}'.format(pathname,catalog_file))
    myargs.append('-tw {}'.format(time_window))
    myargs.append('-dw {}'.format(dist_window))
    myargs=' '.join(myargs)
    args = parser.parse_args(shlex.split(myargs))

    minmag, maxmag = args.magrange

    if args.specifyfile is None:
        #print('Nonew')

        if not args.catalog:
            sys.stdout.write('No catalog specified. Exiting...\n')
            sys.exit()
        elif not args.startyear:
            sys.stdout.write('No starting year specified. Exiting...\n')
            sys.exit()
        elif not args.endyear:
            sys.stdout.write('No ending year specified. Exiting...\n')
            sys.exit()

        catalog = args.catalog.lower()
        startyear, endyear = map(int, [args.startyear, args.endyear])
        download = args.forcedownload

        prefix = '%s%s-%s' % (catalog, startyear, endyear) if catalog else\
                  'preferred%s-%s' % (startyear, endyear)

        if download:
            try:
                os.makedirs(prefix)
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            datadf = qcu.get_data(catalog, prefix, startyear=startyear,
                endyear=endyear, minmag=minmag, maxmag=maxmag)
        else:
            # Python 2
            try:
                try:
                    datadf = pd.read_csv('%s/%s.csv' % (prefix, prefix))
                except IOError:
                    try:
                        os.makedirs(prefix)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                    datadf = qcu.get_data(catalog, prefix,
                        startyear=startyear, endyear=endyear, minmag=minmag,
                        maxmag=maxmag)
            # Python 3
            except:
                try:
                    datadf = pd.read_csv('%s/%s.csv' % (prefix, prefix))
                except FileNotFoundError:
                    try:
                        os.makedirs(prefix)
                    except OSError as exception:
                        if exception.errno != errno.EEXIST:
                            raise
                    datadf = qcu.get_data(catalog, prefix,
                        startyear=startyear, endyear=endyear, minmag=minmag,
                        maxmag=maxmag)

    else:
        from shutil import copy2
        #prefix = '.'.join(args.specifyfile.split('.')[:-1])
        #print(args.specifyfile)

        catalog = args.catalog.lower()
        startyear, endyear = map(int, [args.startyear, args.endyear])
        prefix = '%s%s-%s' % (catalog, startyear, endyear) if catalog else\
                  'preferred%s-%s' % (startyear, endyear)
#        try:
#            os.makedirs(prefix)
#        except OSError as exception:
#            if exception.errno != errno.EEXIST:
#                raise

        #datadf = pd.read_csv(args.specifyfile)
        print(args.specifyfile)
        datadf=qcu.get_local_data(args.specifyfile, prefix, startyear=startyear,
                        endyear=endyear, minmag=minmag, maxmag=maxmag)
        #copy2(args.specifyfile, prefix)

    if len(datadf) == 0:
        sys.stdout.write(('Catalog has no data available for that time period.'
                          ' Quitting...\n'))
        sys.exit()

    timewindow = args.timewindow
    distwindow = args.distwindow

    datadf = datadf.sort_values(by='time').reset_index(drop=True)
    datadf.loc[:, 'ms'] = datadf['time'].str[-4:-1].astype('float')
    datadf.loc[:, 'convtime'] = [' '.join(x.split('T')).split('.')[0]
                                  for x in datadf['time'].tolist()]
    datadf.loc[:, 'convtime'] = datadf['convtime'].astype('datetime64[ns]')
    print(pathname)
    os.chdir(pathname)
    dup1, dup2 = my_list_duplicates(datadf, prefix, timewindow=timewindow,
                                 distwindow=distwindow)
    basic_cat_sum(datadf, prefix, dup1, dup2, timewindow, distwindow)
    largest_ten(datadf, prefix)

    # generate figures
    if (to_plot is None or to_plot['Detect Location']=='true'):
        map_detecs(datadf, prefix, title='Detection locations')
    if(to_plot is None or to_plot['Deeper than 50km']=='true'):
        map_detecs(datadf, prefix, mindep=50, title='Detections deeper than 50km')
    if(to_plot is None or to_plot['Detect Density']=='true'):
        map_detec_nums(datadf, prefix, title='Detection density')
    if(to_plot is None or to_plot['Mag Hist']=='true'):
        make_hist(datadf, 'mag', 0.1, prefix, xlabel='Magnitude',
                  title='Magnitude histogram')
    if(to_plot is None or to_plot['Depth Hist']=='true'):
        make_hist(datadf, 'depth', 1, prefix, xlabel='Depth (km)',
                  title='Depth histogram')
    if(to_plot is None or to_plot['Zoomed Depth Hist']=='true'):
        make_hist(datadf, 'depth', 0.5, prefix, maxval=20, xlabel='Depth (km)',
                  title='Zoomed depth histogram')
#    make_hist(datadf, 'ms', 20, prefix, xlabel='Milliseconds',
#              title='Histogram of milliseconds')
    if(to_plot is None or to_plot['Events per Hour']=='true'):
        make_time_hist(datadf, 'hour', prefix, title='Events per Hour of the Day')
    if(to_plot is None or to_plot['Events per Day']=='true'):
        make_time_hist(datadf, 'day', prefix, title='Events per Day')
    if(to_plot is None or to_plot['Time Seperation']=='true'):
        graph_time_sep(datadf, prefix)
    if(to_plot is None or to_plot['Mag Completeness']=='true'):
        cat_mag_comp(datadf, prefix)
    if(to_plot is None or to_plot['Mag Time']=='true'):
        graph_mag_time(datadf, prefix)
    if(to_plot is None or to_plot['Cumu Moment Release']=='true'):
        cumul_moment_release(datadf, prefix)
    if(to_plot is None or to_plot['Event Types']=='true'):
        graph_event_types(datadf, prefix)
    if(to_plot is None or to_plot['Num of Duplicates']=='true'):
        cat_dup_search(datadf, prefix)
    if(to_plot is None or to_plot['Med Mag']=='true'):
        med_mag(datadf, prefix)
    if(to_plot is None or to_plot['Mag Count']=='true'):
        graph_mag_count(datadf, prefix)

    return pathname,prefix


def generate_html(pathname,prefix,to_show=None):
    """Generate an HTML file containing all of the generated images and text
    files."""
    catalog = prefix[:-9].upper()
    startyear = prefix[-9:-5]
    endyear = prefix[-4:]

    os.chdir(pathname)
    with open('{0}_catalogsummary.txt'.format(prefix)) as sumfile:
        catsum = '\t\t' + '\t\t'.join(sumfile.readlines())
    with open('{0}_largestten.txt'.format(prefix)) as tenfile:
        largest = '\t\t' + '\t\t'.join(tenfile.readlines())
    with open('{0}_duplicates.txt'.format(prefix)) as dupfile:
        duplist = '\t\t' + '\t\t'.join(dupfile.readlines())

    #toc = ('# Report for {1} catalog from {2} to {3}\n'
    toc = ('## 地震目录 {1} 质量报告 {2} — {3}\n'       
           '### 内容\n').format(prefix, catalog, startyear, endyear)
    contents={}
    contents['mulujichuxinxi']=('- [目录基本信息](#catsum)\n'
           '    - [最大事件](#bigevs)\n')
    contents['dizhenhuodong']=('- [地震活动图](#seismap)\n'
           '- [地震密度图](#densmap)\n')
    contents['zhenyuanshendufenbu']='- [震源深部分布图](#dephist)\n'
    contents['dizhenpincitongji']=('- [事件频率](#evfreq)\n'
           '- [每小时的事件频率](#hrevfreq)\n'
           '- [事件之间的事件间隔](#tmsep)\n')
    contents['zhenjitongji']=('- [震级分布图](#magdist)\n'
           '    - [全部震级](#allmag)\n'
           '    - [震级直方图](#allmaghist)\n'
           '    - [分级直方图](#magevcount)\n'
           '    - [中等震级](#medmag)\n'
           '    - [目录完整性](#magcomp)\n')
    contents['leijidizhenjushifang']='- [累积地震矩释放](#cumulrel)\n'
    contents['dizhenleixingtongji']='- [事件类型-频率](#evtypes)\n'
    contents['shaixuanmulu']=('- [重复事件扫描](#dupgraph)\n'
           '- [可能的重复事件](#duplist)\n---\n')
    for key in contents:
        if(to_show is None or to_show[key]=="true"):
           toc=toc+contents[key]
    mdstring=''
    to_show['mulujichuxinxi']="true"
    to_show['dizhenhuodong']="true"
    to_show['zhenyuanshendufenbu']="true"
    to_show['zhenjitongji']="true"
    to_show['dizhenpincitongji']="true"
    to_show['chongfushijian']="true"
    to_show['leijidizhenjushifang']="true"
    to_show['shaixuanmulu']="true"
    to_show['dizhenleixingtongji']="true"
    md_contents={}
    md_contents['mulujichuxinxi'] = ('### Basic Catalog Summary <a name="catsum"></a>\n---\n'
                '{4}\n'
                '#### Largest Events <a name="bigevs"></a>\n---\n'
                '{5}\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['dizhenhuodong'] = ('### Seismicity Map <a name="seismap"></a>\n---\n'
                '<img width="50%" src="{0}_mapdetecs.png">\n'
                '### Seismicity Density Plot <a name="densmap"></a>\n---\n'
                '<img width="50%" src="{0}_eqdensity.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['zhenyuanshendufenbu'] = ('### Depth Distribution <a name="dephist"></a>\n---\n'
                '<img width="50%" src="{0}_depthhistogram.png">\n'
                '<img width="50%" src="{0}_zoomdepthhistogram.png">\n'
                '<img width="50%" src="{0}_morethan50detecs.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['dizhenpincitongji'] = ('### Event Frequency <a name="evfreq"></a>\n---\n'
                '<img width="50%" src="{0}_dayhistogram.png">\n'
                '### Hourly Event Frequency <a name="hrevfreq"></a>\n---\n'
                '<img width="50%" src="{0}_hourhistogram.png">\n'
                '### Inter-Event Temporal Spacing <a name="tmsep"></a>\n---\n'
                '<img width="50%" src="{0}_timeseparation.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['zhenjitongji'] = ('## Magnitude Distribution <a name="magdist"></a>\n'
                '### All Magnitudes <a name="allmag"></a>\n---\n'
                '<img width="50%" src="{0}_magvtime.png">\n'
                '### All Magnitudes Histogram <a name="allmaghist"></a>\n---\n'
                '<img width="50%" src="{0}_maghistogram.png">\n'
                '### Magnitude & Event Count <a name="magevcount"></a>\n---\n'
                '<img width="50%" src="{0}_magtimecount.png">\n'
                '### Median Magnitudes <a name="medmag"></a>\n---\n'
                '<img width="50%" src="{0}_medianmag.png">\n'
                '### Overall Completeness <a name="magcomp"></a>\n---\n'
                '<img width="50%" src="{0}_catmagcomp.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['leijidizhenjushifang'] = ('### Cumulative Moment Release <a name="cumulrel"></a>\n---\n'
                '<img width="50%" src="{0}_cumulmomentrelease.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['dizhenleixingtongji'] = ('### Event Type Frequency <a name="evtypes"></a>\n---\n'
                '<img width="50%" src="{0}_cumuleventtypes.png">\n').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    md_contents['chongfushijian'] = ('### Finding Duplicate Events <a name="dupgraph"></a>\n---\n'
                '<img width="50%" src="{0}_catdupsearch.png">\n'
                '### Possible Duplicate Events <a name="duplist"></a>\n---\n'
                '{6}').format(prefix, catalog, startyear, endyear, catsum, largest,
                        duplist)
    for key in md_contents:
        if(to_show is None or to_show[key]=="true"):
           mdstring=mdstring+md_contents[key]

    html = markdown.markdown(toc + mdstring)

    with open('{0}_report.html'.format(prefix), 'w') as htmlfile:
        htmlfile.write(html)

def qcinit():
    to_plot={}
    to_plot['Detect Location']=to_show['mulujichuxinxi']
    to_plot['Deeper than 50km']=to_show['zhenyuanshendufenbu']
    to_plot['Detect Density']=to_show['dizhenpincitongji']
    to_plot['Mag Hist']=to_show['zhenjitongji']
    to_plot['Depth Hist']=to_show['zhenyuanshendufenbu']
    to_plot['Zoomed Depth Hist']=to_show['zhenyuanshendufenbu']
    to_plot['Events per Hour']=to_show['dizhenpincitongji']
    to_plot['Events per Day']=to_show['dizhenpincitongji']
    to_plot['Time Seperation']=to_show['dizhenpincitongji']
    to_plot['Mag Completeness']=to_show['zhenjitongji']
    to_plot['Cumu Moment Release']=to_show['leijidizhenjushifang']
    to_plot['Num of Duplicates']=to_show['chongfushijian']
    to_plot['Med Mag']=to_show['zhenjitongji']
    to_plot['Mag Count']=to_show['zhenjitongji']
    to_plot['Mag Time']=to_show['zhenjitongji']
    to_plot['Event Types']=to_show['dizhenleixingtongji']  
    
    return to_plot

if __name__ == '__main__':
## picture generation will be guided selected output ##
    to_plot={}
    to_plot['Detect Location']=to_show['mulujichuxinxi']
    to_plot['Deeper than 50km']=to_show['zhenyuanshendufenbu']
    #to_plot['Detect Density']=False
    to_plot['Detect Density']=to_show['dizhenpincitongji']
    to_plot['Mag Hist']=to_show['zhenjitongji']
    to_plot['Depth Hist']=to_show['zhenyuanshendufenbu']
    to_plot['Zoomed Depth Hist']=to_show['zhenyuanshendufenbu']
    to_plot['Events per Hour']=to_show['dizhenpincitongji']
    to_plot['Events per Day']=to_show['dizhenpincitongji']
    to_plot['Time Seperation']=to_show['dizhenpincitongji']
    to_plot['Mag Completeness']=to_show['zhenjitongji']
    to_plot['Cumu Moment Release']=to_show['leijidizhenjushifang']
    to_plot['Num of Duplicates']=to_show['chongfushijian']
    to_plot['Med Mag']=to_show['zhenjitongji']
    to_plot['Mag Count']=to_show['zhenjitongji']
    to_plot['Mag Time']=to_show['zhenjitongji']
    to_plot['Event Types']=to_show['dizhenleixingtongji']
##
    try:
        pathname,prefix = create_figures(to_plot) #画图
        generate_html(pathname,prefix,to_show) #生成html
    except (KeyboardInterrupt, SystemError):
        sys.stdout.write('\nProgram canceled. Exiting...\n')
        sys.exit()
