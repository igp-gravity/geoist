"""
EXAMPLE 3 - ISF CATALOGUE EXPLORATION
"""
from os.path import dirname
from geoist.cattools import Parsers as Par

#-----------------------------------------------------------------------------------------
# Isf parsing
pathname = dirname(__file__)
Db = Par.Database('ISC-HOM')
Db.ImportIsf(pathname+'/data/isc-rev-africa.isf')

#-----------------------------------------------------------------------------------------
# Basic Exploration and filtering

LocCode = ['ISC','NEIC','CSEM']
MagCode = ['GCMT','HRVD','NEIC','ISC','IDC']
MagType = ['MW','MS','ML']

Key, Freq = Db.KeyStat('LocCode')
Key, Freq = Db.KeyStat('MagCode')
Key, Freq = Db.KeyStat('MagType')

Db.Filter('LocCode',LocCode)
Db.Filter('MagCode',MagCode)
Db.Filter('MagType',MagType)

Db.Filter('LocCode',LocCode,Best=1)
Db.Filter('MagCode',MagCode,Best=1)
Db.Filter('MagType',MagType,Best=1)

#-----------------------------------------------------------------------------------------
# Save as CSV
Db.Export(pathname+'/data/isc-rev-africa-select.csv')

# Save as binary
Db.Dump(pathname+'/data/isc-rev-africa-select.bin')
