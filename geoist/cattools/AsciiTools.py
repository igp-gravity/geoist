#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
A simple tool to manipulate data from/to ascii files.
"""

import copy as cp
import numpy as np
import fnmatch as fnm

#-----------------------------------------------------------------------------------------

class AsciiTable():

  def __init__(self, header=[]):

    if header:
      self.header = header
    else:
      self.header = []

    self.data = []

  #---------------------------------------------------------------------------------------

  def AddElement(self, data=[]):
    """
    Add an element (with header's format) to the data structure.
    Element can be empty or filled with data.
    """

    newitem = {}

    for i, key in enumerate(self.header):

      if not data:
        newitem[key] = []
      else:
        newitem[key] = data[i]

    self.data.append(newitem)

  #---------------------------------------------------------------------------------------

  def AddKey(self, key, data=[], index=-1):
    """
    Add header key at given position (default is last element).
    Data structure can optionally be inflated.
    """

    # Default index simply appends
    if index == -1:
      index = len(self.header)

    self.header.insert(index, key)

    # Loop over data
    for i, item in enumerate(self.data):

      # Check value types
      if not data or type(data) != list:
        element = data
      else:
        element = data[i]

      # Add element at corresponding key
      self.data[i][key] = element

  #---------------------------------------------------------------------------------------

  def RemoveKey(self, key):
    """
    Remove a given key from header and data structure.
    """

    # Remove from header
    i = self.header.index(key)
    self.header.pop(i)

    # Remove from data
    for i, item in enumerate(self.data):
      self.data[i].pop(key)

  #---------------------------------------------------------------------------------------

  def RenameKey(self, old_key, new_key):
    """
    Rename a given key in the header and the data structure.
    """

    # Rename header's key
    i = self.header.index(old_key)
    self.header[i] = new_key

    # Rename key in data structure
    for i, item in enumerate(self.data):
      self.data[i][new_key] = self.data[i].pop(old_key)

  #---------------------------------------------------------------------------------------

  def Replace(self, key, old_value, new_value):
    """
    Replace occurences of a key give value.
    If old_value is '*' it replaces all values.
    """

    # Loop over data
    for i, item in enumerate(self.data):

      # Replace all keys
      if old_value == '*':
        self.data[i][key] = new_value
      # Replace matching values only
      else:
        if self.data[i][key] == old_value:
          self.data[i][key] = new_value

  #---------------------------------------------------------------------------------------

  def Size(self):
    """
    Method to return size of the data matrix.
    """

    enum = len(self.data)
    hnum = len(self.header)

    return [enum, hnum]

  #---------------------------------------------------------------------------------------

  def Import(self, ascii_file,
                   header=[],
                   dtype='float',
                   delimiter=',',
                   skipline=0,
                   comment='#',
                   empty=[]):
    """
    Method to import data from ascii file (tabular)
    """

    self.header = []
    self.data = []

    # Open input ascii file
    with open(ascii_file, 'r') as f:

      # Ignore initial line(s) if necessary
      for i in range(0, skipline):
        f.readline()

      # Import header (skip comments)
      if not header:
        while 1:
          line = f.readline()
          if line[0] != comment: break
        header = line.strip().split(delimiter)

      # Removing empty fields from header
      for h in header:
        if h != '':
          self.header.append(h)

      # Loop over lines
      for line in f:

         # Skip comments, if any
        if line[0] != comment:
          value = line.strip().split(delimiter)

          # Loop over data values
          data = []
          for i, h in enumerate(header):

            # Skip empty header fields
            if h != '':

              # Data type(s) switch
              if type(dtype) == list:
                dtp = dtype[i]
              else:
                dtp = dtype

              # Check for empty elements
              if not value[i]:
                value[i] = empty

              data.append(_CastValue(value[i],dtp))

          self.AddElement(data)

      f.close()
      return

    # Warn user if model file does not exist
    print('File not found.')
  #---------------------------------------------------------------------------------------
  def ExportEQT(self, ascii_file,
                   write_header='no',
                   delimiter=' '):
    """
    Method to export data object into an ascii file.
    """

    try:
      with open(ascii_file, 'w') as f:

        # Write header
        if write_header == 'yes':
          header = delimiter.join(self.header)
          f.write(header + '\n')

        # Write data (loop over rows)
        for i, item in enumerate(self.data):
          data = [_CastValue(item[j],'s') for j in self.header]
          data = delimiter.join(data)

          if i < (self.Size()[0]-1):
            f.write(data + '\n')
          else:
            f.write(data)

        f.close()

    except:
      # Warn user if model file does not exist
      print('File not found.')

  #---------------------------------------------------------------------------------------


  #---------------------------------------------------------------------------------------

  def Export(self, ascii_file,
                   write_header='yes',
                   delimiter=','):
    """
    Method to export data object into an ascii file.
    """

    try:
      with open(ascii_file, 'w') as f:

        # Write header
        if write_header == 'yes':
          header = delimiter.join(self.header)
          f.write(header + '\n')

        # Write data (loop over rows)
        for i, item in enumerate(self.data):
          data = [_CastValue(item[j],'s') for j in self.header]
          data = delimiter.join(data)

          if i < (self.Size()[0]-1):
            f.write(data + '\n')
          else:
            f.write(data)

        f.close()

    except:
      # Warn user if model file does not exist
      print('File not found.')

  #---------------------------------------------------------------------------------------

  def Append(self, new_table):
    """
    Method to merge two data structures consecutively.
    (Header structure must be identical).
    """

    if self.header == new_table.header:
      for i in range(0,new_table.Size()[0]):
        self.data.append(new_table.data[i])

    else:
      print('Error: headers do not match...')

  #---------------------------------------------------------------------------------------

  def Extract(self, key, dtype='float'):
    """
    Method to extract data values by key.
    Data type can be specified.
    """

    values = []

    for item in self.data:
      value = _CastValue(item[key], dtype)
      values.append(value)

    return values

  #---------------------------------------------------------------------------------------

  def Filter(self, key, filter_key):
    """
    Method to filter the data table by key value.
    Value can be a string to amtch (* and ? allowed)
    or a numerical range (as a list of floats).
    In output it is returned a new table.
    """

    NewTab = AsciiTable(self.header)

    # String matching
    if type(filter_key) is str:

      for item in self.data:
        if fnm.fnmatch(item[key],filter_key):
           NewTab.data.append(item)

    # Filter by value
    if type(filter_key) is list:

      for item in self.data:

        if not _isNaN(item[key]):
          ik = _CastValue(item[key])

          if ik >= filter_key[0] and ik <= filter_key[1]:
            NewTab.data.append(item)

    return NewTab

#-----------------------------------------------------------------------------------------

def _CastValue(value, dtype='float'):
  """
  Private method to recast variables.
  """

  # Check for empty fields or nans
  if not _isEmpty(value) and not _isNaN(value):

    # Data casting
    if dtype in ['Int','int','I','i']:
      value = int(value)

    if dtype in ['Float','float','F','f']:
      value = float(value)

    if dtype in ['String','string','S','s']:
      value = str(value)

  else:

    # Using default strings
    if dtype in ['String','string','S','s']:

      if _isEmpty(value):
        value = ''

      if _isNaN(value):
        value = 'nan'

  return value

#-----------------------------------------------------------------------------------------

def _isNaN(number):
  """
  Simple private method to check if a number is NaN.
  It returns a boolean evaluation.
  """
  return number != number

def _isEmpty(number):
  """
  Simple private method to check if a variable (list) is empty.
  It returns a boolean evaluation.
  """
  return (number == [] or number == '')