# -*- coding: utf-8 -*-
"""
 Name        : igrf.py
 Created on  : 2018/09/11 17:00
 Author      : Steve Chen <chenshi@cea-igp.ac.cn>
 Affiliation : Institute of Geophysics, CEA.
 Version     : 0.1.0
 Copyright   : Copyright (C) 2018-2020 GEOIST Development Team. All Rights Reserved.
 License     : Distributed under the MIT License. See LICENSE.txt for more info.
 Github      : https://igp-gravity.github.io/
 Description : Application for ***.
"""


#This code is based on the Geomag software, with information given below. It was
#translated into Python from the Geomag code.
"""
| This program, originally written in FORTRAN, was developed using subroutines
| written by   : A. Zunde
|                USGS, MS 964, Box 25046 Federal Center, Denver, Co.  80225
|                and
|                S.R.C. Malin & D.R. Barraclough
|                Institute of Geological Sciences, United Kingdom.

| Translated
| into C by    : Craig H. Shaffer
|                29 July, 1988

| Rewritten by : David Owens
|                For Susan McLean

| Maintained by: Stefan Maus
| Contact      : stefan.maus@noaa.gov
|                National Geophysical Data Center
|                World Data Center-A for Solid Earth Geophysics
|                NOAA, E/GC1, 325 Broadway,
|                Boulder, CO  80303
"""
#stdlib imports
from math import sin
from math import cos
from math import sqrt
from math import atan2
import os
#third party imports
import numpy as np
# local imports t

class IGRF():
    """
    IGRF field calculation

    This produces two datasets. The first is an IGRF dataset for the area of
    interest, defined by some input magnetic dataset. The second is the IGRF
    corrected form of that input magnetic dataset.

    To do this, the input dataset must be reprojected from its local projection
    to degrees, where the IGRF correction will take place. This is done within
    this class.

    Attributes
    ----------
    parent : parent
        reference to the parent routine
    indata : dictionary
        dictionary of input datasets
    outdata : dictionary
        dictionary of output datasets

    Parameters
    ----------
    altmin : Double
        Minimum height of selected model.
    altmax : Double array
        array of MAXMOD Maximum height of model.
    maxalt : Double
        Maximum height of selected model.
    d : float
        Declination of the field from the geographic north (deg).
    sdate : float
        start date inputted
    ddot : float
        annual rate of change of decl. (arc-min/yr)
    alt : float
        altitude above WGS84 Ellipsoid
    epoch : list
        list of MAXMOD epoch of model.
    latitude : float
        Latitude.
    longitude : float
        Longitude.
    gh : numpy array
        Schmidt quasi-normal internal spherical harmonic coeff.
        Schmidt quasi-normal internal spherical harmonic coeff.
        Coefficients of resulting model.
        Coefficients of rate of change model.
    i : float
        Inclination (deg).
    idot : float
        Rate of change of i (arc-min/yr).
    igdgc : int
        Flag for geodetic or geocentric coordinate choice.
    irec_pos : int array
        array of MAXMOD Record counter for header
    fileline : int
        Current line in file (for errors)
    max1 : list, int
        array of MAXMOD Main field coefficient.
    max2 : list, int
        array of MAXMOD Secular variation coefficient.
    max3 : list, int
        array of MAXMOD Acceleration coefficient.
    minyr : float
        Min year of all models
    maxyr : float
        Max year of all models
    yrmax : list, float
        array of MAXMOD  Max year of model.
    yrmin : list, float
        array of MAXMOD  Min year of model.
    """
    def __init__(self, parent=None):
        self.indata = {}
        self.outdata = {}
        self.reportback = {}
        MAXDEG = 13
        MAXCOEFF = (MAXDEG*(MAXDEG+2)+1)

        self.gh = np.zeros([4, MAXCOEFF])
        self.d = 0
        self.f = 0
        self.h = 0
        self.i = 0
        self.dtemp = 0
        self.ftemp = 0
        self.htemp = 0
        self.itemp = 0
        self.x = 0
        self.y = 0
        self.z = 0
        self.xtemp = 0
        self.ytemp = 0
        self.ztemp = 0
        self.xdot = 0
        self.ydot = 0
        self.zdot = 0
        self.fdot = 0
        self.hdot = 0
        self.idot = 0
        self.ddot = 0
        self.dsb_alt = 0
        self.date = ''

    def pnt(self, latitude, longitude, alt, Model='IGRF12'):
        """
        Settings Dialog. This is the main entrypoint into this routine. It also
        contains the main IGRF code.
        """
# Variable declaration
# Control variables
        if Model == 'IGRF12':
          mdf = open(os.path.abspath(__file__).rpartition('\\')[0]+'\\IGRF12.cof')
        else:
          mdf = open(os.path.abspath(__file__).rpartition('\\')[0]+'\\IGRF11.cof')
                     
        #mdf = open(__file__.rpartition('\\')[0]+'\\IGRF12.cof')
        modbuff = mdf.readlines()
        fileline = -1                            # First line will be 1
        model = []
        epoch = []
        max1 = []
        max2 = []
        max3 = []
        yrmin = []
        yrmax = []
        altmin = []
        altmax = []
        irec_pos = []
# First model will be 0
        for i in modbuff:
            fileline += 1  # On new line
            if i[:3] == '   ':
                i2 = i.split()
                model.append(i2[0])
                epoch.append(float(i2[1]))
                max1.append(int(i2[2]))
                max2.append(int(i2[3]))
                max3.append(int(i2[4]))
                yrmin.append(float(i2[5]))
                yrmax.append(float(i2[6]))
                altmin.append(float(i2[7]))
                altmax.append(float(i2[8]))
                irec_pos.append(fileline)
                

        # Pick model
        sdate = 2018
        yrmax = np.array(yrmax)
        modelI = sum(yrmax < sdate)
        igdgc = 1
        #print(yrmax, sdate,modelI)
        #print(max2)
        if max2[modelI] == 0:
            self.getshc(modbuff, 1, irec_pos[modelI], max1[modelI], 0)
            self.getshc(modbuff, 1, irec_pos[modelI+1], max1[modelI+1], 1)
            nmax = self.interpsh(sdate, yrmin[modelI], max1[modelI],
                                 yrmin[modelI+1], max1[modelI+1], 2)
            nmax = self.interpsh(sdate+1, yrmin[modelI], max1[modelI],
                                 yrmin[modelI+1], max1[modelI+1], 3)
        else:
            self.getshc(modbuff, 1, irec_pos[modelI], max1[modelI], 0)
            self.getshc(modbuff, 0, irec_pos[modelI], max2[modelI], 1)
            nmax = self.extrapsh(sdate, epoch[modelI], max1[modelI],
                                 max2[modelI], 2)
            nmax = self.extrapsh(sdate+1, epoch[modelI], max1[modelI],
                                 max2[modelI], 3)
#
#
## Do the first calculations
            self.shval3(igdgc, latitude, longitude, alt, nmax, 3)
            self.dihf(3)

        self.reportback['Int']=str(self.f)
        self.reportback['Inc']=str(np.rad2deg(self.i))
        self.reportback['Dec']=str(np.rad2deg(self.d))
        #print(self.reportback)
        return [self.f, np.rad2deg(self.i), np.rad2deg(self.d)]

    def getshc(self, file, iflag, strec, nmax_of_gh, gh):
        """
        Reads spherical harmonic coefficients from the specified model into an
        array.

        | FORTRAN: Bill Flanagan, NOAA CORPS, DESDIS, NGDC, 325 Broadway,
        | Boulder CO.  80301
        | C: C. H. Shaffer, Lockheed Missiles and Space Company, Sunnyvale CA

        Parameters
        ----------
        file : file
            reference to a file object
        iflag :
            Flag for SV equal to 1 or not equal to 1 for designated read
            statements
        strec : int
            Starting record number to read from model
        nmax_of_gh : int
            Maximum degree and order of model

        Returns
        -------
        gh : list
            Schmidt quasi-normal internal spherical harmonic coefficients
        """
        ii = -1
        cnt = 0

        for nn in range(1, nmax_of_gh+1):
            for _ in range(nn+1):
                cnt += 1
                tmp = file[strec+cnt]
                tmp = tmp.split()
                m = int(tmp[1])

                if iflag == 1:
                    g = float(tmp[2])
                    hh = float(tmp[3])
                else:
                    g = float(tmp[4])
                    hh = float(tmp[5])

                ii = ii + 1
                self.gh[gh][ii] = g

                if m != 0:
                    ii = ii + 1
                    self.gh[gh][ii] = hh

        return

    def extrapsh(self, date, dte1, nmax1, nmax2, gh):
        """
        Extrapolates linearly a spherical harmonic model with a
        rate-of-change model.

        | FORTRAN : A. Zunde, USGS, MS 964, box 25046 Federal Center, Denver,
        | CO. 80225
        | C : C. H. Shaffer, Lockheed Missiles and Space Company, Sunnyvale CA

        Parameters
        ----------
        date : float
            date of resulting model (in decimal year)
        dte1 : float
            date of base model
        nmax1 : int
            maximum degree and order of base model
        gh  : numpy array
            Schmidt quasi-normal internal spherical harmonic coefficients of
            base model and rate-of-change model
        nmax2 : int
            maximum degree and order of rate-of-change model

        Returns
        -------
        gh : numpy array
            Schmidt quasi-normal internal spherical harmonic coefficients
        nmax : int
            maximum degree and order of resulting model
        """
        factor = date - dte1
        if nmax1 == nmax2:
            k = nmax1 * (nmax1 + 2)
            nmax = nmax1
        else:
            if nmax1 > nmax2:
                k = nmax2 * (nmax2 + 2)
                l = nmax1 * (nmax1 + 2)
                for ii in range(k, l):
                    self.gh[gh][ii] = self.gh[0][ii]

                nmax = nmax1
            else:
                k = nmax1 * (nmax1 + 2)
                l = nmax2 * (nmax2 + 2)
                for ii in range(k, l):
                    self.gh[gh][ii] = factor * self.gh[1][ii]

                nmax = nmax2

        for ii in range(k):
            self.gh[gh][ii] = self.gh[0][ii] + factor * self.gh[1][ii]

        return nmax

    def interpsh(self, date, dte1, nmax1, dte2, nmax2, gh):
        """
        Interpolates linearly, in time, between two spherical harmonic
        models.

        | FORTRAN : A. Zunde, USGS, MS 964, box 25046 Federal Center, Denver,
        | CO. 80225
        | C : C. H. Shaffer, Lockheed Missiles and Space Company, Sunnyvale CA

        Parameters
        ----------
        date : float
            date of resulting model (in decimal year)
        dte1 : float
            date of earlier model
        nmax1 : int
            maximum degree and order of earlier model
        gh : numpy array
            Schmidt quasi-normal internal spherical harmonic coefficients of
            earlier model and internal model
        dte2 : float
            date of later model
        nmax2 : int
            maximum degree and order of later model

        Returns
        -------
        gh : numpy array
            coefficients of resulting model
        nmax : int
            maximum degree and order of resulting model
        """
        factor = (date - dte1) / (dte2 - dte1)
        if nmax1 == nmax2:
            k = nmax1 * (nmax1 + 2)
            nmax = nmax1
        else:
            if nmax1 > nmax2:
                k = nmax2 * (nmax2 + 2)
                l = nmax1 * (nmax1 + 2)
                for ii in range(k, l):
                    self.gh[gh][ii] = self.gh[0][ii] + factor*(-self.gh[0][ii])
                nmax = nmax1
            else:
                k = nmax1 * (nmax1 + 2)
                l = nmax2 * (nmax2 + 2)
                for ii in range(k, l):
                    self.gh[gh][ii] = factor * self.gh[1][ii]

                nmax = nmax2

        for ii in range(k):
            self.gh[gh][ii] = self.gh[0][ii] + factor*(self.gh[1][ii] -
                                                       self.gh[0][ii])

        return nmax

    def shval3(self, igdgc, flat, flon, elev, nmax, gh):
        """
        Calculates field components from spherical harmonic (sh) models.

        Based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,
        report no. 71/1, institute of geological sciences, U.K.

        | FORTRAN : Norman W. Peddie, USGS, MS 964, box 25046 Federal Center,
        | Denver, CO. 80225
        | C : C. H. Shaffer, Lockheed Missiles and Space Company, Sunnyvale CA

        Parameters
        ----------
        igdgc : int
            indicates coordinate system used set equal to 1 if geodetic, 2 if
            geocentric
        latitude : float
            north latitude, in degrees
        longitude : float
            east longitude, in degrees
        elev : float
            WGS84 altitude above ellipsoid (igdgc=1), or radial distance from
            earth's center (igdgc=2)
        a2,b2 : float
            squares of semi-major and semi-minor axes of the reference spheroid
            used for transforming between geodetic and geocentric coordinates
            or components
        nmax : int
            maximum degree and order of coefficients

        Returns
        -------
        x : float
            northward component
        y : float
            eastward component
        z : float
            vertically-downward component
        """

        sl = np.zeros(14)
        cl = np.zeros(14)
        p = np.zeros(119)
        q = np.zeros(119)
        earths_radius = 6371.2
        dtr = np.pi/180.0
        a2 = 40680631.59            # WGS84
        b2 = 40408299.98            # WGS84
        r = elev
        argument = flat * dtr
        slat = sin(argument)
        if (90.0 - flat) < 0.001:
            aa = 89.999            # 300 ft. from North pole
        elif (90.0 + flat) < 0.001:
            aa = -89.999        # 300 ft. from South pole
        else:
            aa = flat

        argument = aa * dtr
        clat = cos(argument)
        argument = flon * dtr
        sl[1] = sin(argument)
        cl[1] = cos(argument)

        if gh == 3:
            self.x = 0
            self.y = 0
            self.z = 0
        elif gh == 4:
            self.xtemp = 0
            self.ytemp = 0
            self.ztemp = 0

        sd = 0.0
        cd = 1.0
        l = 0
        n = 0
        m = 1
        npq = int((nmax * (nmax + 3)) / 2)
        if igdgc == 1:
            aa = a2 * clat * clat
            bb = b2 * slat * slat
            cc = aa + bb
            argument = cc
            dd = sqrt(argument)
            argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc
            r = sqrt(argument)
            cd = (elev + dd) / r
            sd = (a2 - b2) / dd * slat * clat / r
            aa = slat
            slat = slat * cd - clat * sd
            clat = clat * cd + aa * sd

        ratio = earths_radius / r
        argument = 3.0
        aa = sqrt(argument)
        p[1] = 2.0 * slat
        p[2] = 2.0 * clat
        p[3] = 4.5 * slat * slat - 1.5
        p[4] = 3.0 * aa * clat * slat
        q[1] = -clat
        q[2] = slat
        q[3] = -3.0 * clat * slat
        q[4] = aa * (slat * slat - clat * clat)
        for k in range(1, npq+1):
            if n < m:
                m = 0
                n = n + 1
                argument = ratio
                power = n + 2
                rr = pow(argument, power)
                fn = n

            fm = m
            if k >= 5:
                if m == n:
                    argument = (1.0 - 0.5/fm)
                    aa = sqrt(argument)
                    j = k - n - 1
                    p[k] = (1.0 + 1.0/fm) * aa * clat * p[j]
                    q[k] = aa * (clat * q[j] + slat/fm * p[j])
                    sl[m] = sl[m-1] * cl[1] + cl[m-1] * sl[1]
                    cl[m] = cl[m-1] * cl[1] - sl[m-1] * sl[1]
                else:
                    argument = fn*fn - fm*fm
                    aa = sqrt(argument)
                    argument = ((fn - 1.0)*(fn-1.0)) - (fm * fm)
                    bb = sqrt(argument)/aa
                    cc = (2.0 * fn - 1.0)/aa
                    ii = k - n
                    j = k - 2 * n + 1
                    p[k] = (fn + 1.0) * (cc * slat/fn * p[ii]-bb/(fn-1.0)*p[j])
                    q[k] = cc * (slat * q[ii] - clat/fn * p[ii]) - bb * q[j]

            if gh == 3:
                aa = rr * self.gh[2][l]
            elif gh == 4:
                aa = rr * self.gh[3][l]

            if m == 0:
                if gh == 3:
                    self.x = self.x + aa * q[k]
                    self.z = self.z - aa * p[k]
                elif gh == 4:
                    self.xtemp = self.xtemp + aa * q[k]
                    self.ztemp = self.ztemp - aa * p[k]
                else:
                    print("\nError in subroutine shval3")

                l = l + 1
            else:
                if gh == 3:
                    bb = rr * self.gh[2][l+1]
                    cc = aa * cl[m] + bb * sl[m]
                    self.x = self.x + cc * q[k]
                    self.z = self.z - cc * p[k]
                    if clat > 0:
                        self.y = (self.y + (aa*sl[m] - bb*cl[m])*fm *
                                  p[k]/((fn + 1.0)*clat))
                    else:
                        self.y = self.y + (aa*sl[m] - bb*cl[m])*q[k]*slat
                    l = l + 2
                elif gh == 4:
                    bb = rr * self.gh[3][l+1]
                    cc = aa * cl[m] + bb * sl[m]
                    self.xtemp = self.xtemp + cc * q[k]
                    self.ztemp = self.ztemp - cc * p[k]
                    if clat > 0:
                        self.ytemp = (self.ytemp + (aa*sl[m] - bb*cl[m])*fm *
                                      p[k]/((fn + 1.0) * clat))
                    else:
                        self.ytemp = (self.ytemp + (aa*sl[m] - bb*cl[m]) *
                                      q[k]*slat)
                    l = l + 2

            m = m + 1

        if gh == 3:
            aa = self.x
            self.x = self.x * cd + self.z * sd
            self.z = self.z * cd - aa * sd
        elif gh == 4:
            aa = self.xtemp
            self.xtemp = self.xtemp * cd + self.ztemp * sd
            self.ztemp = self.ztemp * cd - aa * sd

    def dihf(self, gh):
        """
        Computes the geomagnetic d, i, h, and f from x, y, and z.

        | FORTRAN : A. Zunde, USGS, MS 964, box 25046 Federal Center, Denver,
        | CO. 80225
        | C : C. H. Shaffer, Lockheed Missiles and Space Company, Sunnyvale CA

        Parameters
        ----------
        x : float
            northward component
        y : float
            eastward component
        z : float
            vertically-downward component

        Returns
        -------
        d : float
            declination
        i : float
            inclination
        h : float
            horizontal intensity
        f : float
            total intensity
        """
        sn = 0.0001

        if gh == 3:
            x = self.x
            y = self.y
            z = self.z
            h = self.h
            f = self.f
            i = self.i
            d = self.d
        else:
            x = self.xtemp
            y = self.ytemp
            z = self.ztemp
            h = self.htemp
            f = self.ftemp
            i = self.itemp
            d = self.dtemp

        for _ in range(2):
            h2 = x*x + y*y
            argument = h2
            h = sqrt(argument)       # calculate horizontal intensity
            argument = h2 + z*z
            f = sqrt(argument)      # calculate total intensity
            if f < sn:
                d = np.nan        # If d and i cannot be determined,
                i = np.nan        # set equal to NaN
            else:
                argument = z
                argument2 = h
                i = atan2(argument, argument2)
                if h < sn:
                    d = np.nan
                else:
                    hpx = h + x
                    if hpx < sn:
                        d = np.pi
                    else:
                        argument = y
                        argument2 = hpx
                        d = 2.0 * atan2(argument, argument2)

        if gh == 3:
            self.h = h
            self.f = f
            self.i = i
            self.d = d
        else:
            self.htemp = h
            self.ftemp = f
            self.itemp = i
            self.dtemp = d

if __name__ == "__main__":
 m1=IGRF()
 m1.pnt(45,105,100)