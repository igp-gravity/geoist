#!/usr/bin/env python
import os
import time
import errno
import argparse
from datetime import datetime
from math import radians, cos, sqrt

import numpy as np
import pandas as pd

###############################################################################
###############################################################################
###############################################################################

def eq_dist(lat1, lon1, lat2, lon2):
    """Calculate equirectangular distance between two points."""
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])
    xval = (lon2 - lon2) * cos(0.5 * (lat2+lat1))
    yval = lat2 - lat1
    eqd = 6371 * sqrt(xval*xval + yval*yval)

    return eqd


def to_epoch(ogtime):
    """Convert from catalog time format to Unix/epoch time."""
    fstr = '%Y-%m-%dT%H:%M:%S.%fZ'
    epoch = datetime(1970, 1, 1)
    epochtime = (datetime.strptime(ogtime, fstr) - epoch).total_seconds()

    return epochtime


def detect_duplicates(catfile, otwindow=16, distwindow=100, diffnet=True):
    """Return possible duplicate events within a catalog."""
    dupefilename = '%s_duplicates.csv' % '.'.join(catfile.split('.')[:-1])
    cat = pd.read_csv(catfile).sort_values(by='time')
    cat.loc[:, 'time'] = [to_epoch(time) for time in cat['time'].tolist()]

    splitlist = []

    for i in range(len(cat)):
        e2ix = cat[(abs(cat['time'] - cat.ix[i]['time']) <= otwindow)
                    & (cat.index > i)].index.values

        if len(e2ix) != 0:
            C = np.array([eq_dist(cat.ix[i]['latitude'],
                                  cat.ix[i]['longitude'],
                                  cat.ix[x]['latitude'],
                                  cat.ix[x]['longitude'])
                          for x in e2ix])

            inds = np.where(C < distwindow)[0]

            if len(inds) != 0:
                for ind in inds:
                    e1id = cat.ix[i]['id']
                    e1net = cat.ix[i]['net']
                    e2id = cat.ix[e2ix[ind]]['id']
                    e2net = cat.ix[e2ix[ind]]['net']

                    if diffnet:
                        if (e1id != e2id) and (e1net != e2net):
                            splitlist.append([str(e1id), str(e2id)])
                    else:
                        if (e1id != e2id):
                            splitlist.append([str(e1id), str(e2id)])

    outdf = pd.DataFrame(splitlist, columns=['id1', 'id2'])
    outdf.to_csv(dupefilename, index=False)


###############################################################################
###############################################################################
###############################################################################


def main():
    """Main function. Defines command line arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument('catalogfile', nargs='?', type=str, default='',
                        help='pick which catalog .csv file to find \
                        duplicates in')

    parser.add_argument('-tw', '--timewindow', nargs='?', default=16, type=int,
                        help='set time window for finding duplicate events in \
                        seconds; default is 16')
    parser.add_argument('-dw', '--distancewindow', nargs='?', default=100,
                        type=int, help='set distance window for finding \
                        duplicate events in kilometers; default is 100')
    parser.add_argument('-df', '--differentnetwork', action='store_true',
                        help='make it so that splits are only found between \
                        different networks')

    args = parser.parse_args()

    catfile = args.catalogfile
    otwindow, distwindow = args.timewindow, args.distancewindow

    detect_duplicates(catfile, otwindow=otwindow, distwindow=distwindow,
                      diffnet=args.differentnetwork)

if __name__ == "__main__":

    main()
