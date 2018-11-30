#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Module for Specific Catalogue Parsers
"""

import math as mt
import .Catalogue as Cat
import .AsciiTools as AT
import .CatUtils as CU

#-----------------------------------------------------------------------------------------

class Database(Cat.Database):

  def __init__(self, Name=[], Info=[]):

    super(Database, self).__init__(Name=Name, Info=Info)

  #---------------------------------------------------------------------------------------

  def ImportIsf(self, file_name):

    def _GetHeader(f):
      Header = f.readline().strip()
      IscType = f.readline().strip()
      return IscType

    def _SplitEvents(f):
      Event = ''
      for Line in f:
        Event += Line
      Event = Event.split('Event ')
      return Event

    def _SplitBlock(Event):
      Id = []
      LocBlock = []
      MagBlock = []
      if Event:
        Rows = Event.split('\n')
        Id = Rows[0].split()[0]
        for R in Rows[1:]:
          if  '(' not in R:
            if (len(R) == 136) and ('Date' not in R):
              LocBlock.append(_GetLocation(R))
            if (len(R) == 38) and ('Magnitude' not in R):
              MagBlock.append(_GetMagnitude(R))
          if '#PRIME' in R:
            LocBlock[-1]['Prime'] = True
      return Id, LocBlock, MagBlock

    def _GetLocation(IsfLoc):
      L = {}
      if IsfLoc:
        L['Year'] = IsfLoc[0:4].strip(' ')
        L['Month'] = IsfLoc[5:7].strip(' ')
        L['Day'] = IsfLoc[8:10].strip(' ')
        L['Hour'] = IsfLoc[10:13].strip(' ')
        L['Minute'] = IsfLoc[14:16].strip(' ')
        L['Second'] = IsfLoc[17:22].strip(' ')
        L['Latitude'] = IsfLoc[36:44].strip(' ')
        L['Longitude'] = IsfLoc[45:54].strip(' ')
        L['Depth'] = IsfLoc[71:76].strip(' ')
        L['SecError'] = IsfLoc[24:29].strip(' ')
        L['DepError'] = IsfLoc[78:82].strip(' ')
        L['LocCode'] = IsfLoc[118:127].strip(' ')
      return L

    def _GetMagnitude(IsfMag):
      M = {}
      if IsfMag:
        M['MagType'] = IsfMag[0:5].strip(' ')
        M['MagSize'] = IsfMag[6:10].strip(' ')
        M['MagError'] = IsfMag[11:14].strip(' ')
        M['MagCode'] = IsfMag[20:29].strip(' ')
      return M

    # Open ISF file
    with open(file_name, 'r') as f:

      _GetHeader(f)
      Event = _SplitEvents(f)

      for E in Event:
        Id, LocBlock, MagBlock = _SplitBlock(E)
        if Id:
          self.AddEvent(Id, LocBlock, MagBlock)

      f.close()
      return

    # Warn user if model file does not exist
    print 'File not found.'

  #---------------------------------------------------------------------------------------

  def ImportNdk(self, file_name):

    def _GetLocation(HypoStr):
      L = {}
      if HypoStr:
        L['Year'] = HypoStr[5:9].strip(' ')
        L['Month'] = HypoStr[10:12].strip(' ')
        L['Day'] = HypoStr[13:15].strip(' ')
        L['Hour'] = HypoStr[16:18].strip(' ')
        L['Minute'] = HypoStr[19:21].strip(' ')
        L['Second'] = HypoStr[22:26].strip(' ')
        L['Latitude'] = HypoStr[27:33].strip(' ')
        L['Longitude'] = HypoStr[34:41].strip(' ')
        L['Depth'] = HypoStr[42:47].strip(' ')
        L['LocCode'] = HypoStr[0:4].strip(' ')
      return L

    def _GetId(CmtStr):
      I = []
      if CmtStr:
        I = CmtStr[0:16].strip(' ')
      return I

    def _GetExponent(CmtStr):
      E = []
      if CmtStr:
        E = float(CmtStr[0:2])
      return E

    def _GetMagnitude(CmtStr, E):
      M = {}
      if CmtStr:
        # Hanks & Kanamori equation
        Moment = float(CmtStr[49:56])*(10.**E)
        Mw = ((2./3.)*mt.log10(Moment))-10.7

        M['MagSize'] = float(format(Mw,'.2f'))
        M['MagType'] = 'MW'
        M['MagError'] = 0.
        M['MagCode'] = 'GCMT-NDK'
      return M

    def _GetMomentTensor(CmtStr):
      T = {}
      if CmtStr:
        T['Mrr'] = float(CmtStr[3:9].strip(' '))
        T['Mtt'] = float(CmtStr[16:22].strip(' '))
        T['Mpp'] = float(CmtStr[29:35].strip(' '))
        T['Mrt'] = float(CmtStr[42:48].strip(' '))
        T['Mrp'] = float(CmtStr[55:61].strip(' '))
        T['Mtp'] = float(CmtStr[68:74].strip(' '))
      return T

    # Open NDK file
    with open(file_name, 'r') as f:

      self.Tensor = []

      while True:

        HypoStr = f.readline().strip('\n')
        Cmt1Str = f.readline().strip('\n')
        Cmt2Str = f.readline().strip('\n')
        Cmt3Str = f.readline().strip('\n')
        Cmt4Str = f.readline().strip('\n')

        L = _GetLocation(HypoStr)
        I = _GetId(Cmt1Str)
        E = _GetExponent(Cmt3Str)
        T = _GetMomentTensor(Cmt3Str)
        M = _GetMagnitude(Cmt4Str, E)

        if I:
          self.AddEvent(I, L, M)
          self.Tensor.append(T)
        else:
          break

      f.close()
      return

    # Warn user if model file does not exist
    print 'File not found.'

  #---------------------------------------------------------------------------------------

  def ImportIgn(self, file_name):

    def _GetId(EventStr):
      I = []
      if EventStr:
        I = EventStr[0:12].strip(' ')
      return I

    def _GetLocation(EventStr):
      L = {}
      if EventStr:
        L['Year'] = EventStr[23:27].strip(' ')
        L['Month'] = EventStr[20:22].strip(' ')
        L['Day'] = EventStr[17:19].strip(' ')
        L['Hour'] = EventStr[34:36].strip(' ')
        L['Minute'] = EventStr[37:39].strip(' ')
        L['Second'] = EventStr[40:42].strip(' ')
        L['Latitude'] = EventStr[42:57].strip(' ')
        L['Longitude'] = EventStr[57:72].strip(' ')
        L['Depth'] = EventStr[72:87].strip(' ')
        L['LocCode'] = 'IGN'
      return L

    def _GetMagnitude(EventStr):
      M = {}
      if EventStr:
        M['MagSize'] = EventStr[102:117].strip(' ')
        M['MagError'] = 0.
        M['MagCode'] = 'IGN'
        if EventStr[131] == '1':
          M['MagType'] = 'MDs'
        if EventStr[131] == '2':
          M['MagType'] = 'MbLg'
        if EventStr[131] == '3':
          M['MagType'] = 'mb'
        if EventStr[131] == '4':
          M['MagType'] = 'mbLg'
        if EventStr[131] == '5':
          M['MagType'] = 'Mw'
        if EventStr[131] == '6':
          M['MagType'] = 'MLv'
        if EventStr[131] == '7':
          M['MagType'] = 'mb'
        if EventStr[131] == '8':
          M['MagType'] = 'mB'
        if EventStr[131] == '9':
          M['MagType'] = 'Mwp'
        if EventStr[131] == '10':
          M['MagType'] = 'MwmB'
        if EventStr[131] == '11':
          M['MagType'] = 'MwMwp'
      return M

    def _GetLog(EventStr):
      O = []
      if EventStr:
        O = 'INFO({0});'.format(EventStr[135:-1])
      return O

    # Open IGN file
    with open(file_name, 'r') as f:

      # Read header
      Header = f.readline().strip('\n')

      while True:

        EventStr = f.readline().strip('\n')

        I = _GetId(EventStr)
        L = _GetLocation(EventStr)
        M = _GetMagnitude(EventStr)
        O = _GetLog(EventStr)

        if I:
          self.AddEvent(I, L, M, O)
        else:
          break

      f.close()
      return

    # Warn user if model file does not exist
    print 'File not found.'

  #---------------------------------------------------------------------------------------

  def ImportIscGemExt(self, file_name):
  
    # Defining the headers/estruture
    Header = ['Id','LocCode','MagCode','','Year','Month','Day','Hour','Minute','Second',
              'SecError','Longitude','Latitude','','','','Depth','DepError','MagSize',
              'MagError']

    def _GetMagCode(MagCode):
        try:
           MC = MagCode.split('|')[1]
        except:
           MC = MagCode
        return MC
    
    # Open/importing ISC_GEM_EXT csv file
    tab = AT.AsciiTable()
    tab.Import(file_name, header=Header,
                         delimiter=',',
                         skipline=1,
                         comment='#',
                         dtype='s')

    if tab.data:
      for I,D in enumerate(tab.data):
            if 'Id' in D.keys():
               I = D['Id']
            else:
               I += 1            
            # Changing the Magcode [Identifier]
            D['MagCode'] = _GetMagCode(D['MagCode'])
            L = CU.LocationInit()
            M = CU.MagnitudeInit()
      
            for K in tab.header:
                if K in L:
                    L[K] = D[K]
                if K in M:
                    M[K] = D[K]
            if 'Log' in D.keys():
                O = D['Log']
            else:
                O = ''
            self.AddEvent(I, L, M, O)
