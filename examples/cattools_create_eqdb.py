"""
MANUAL CREATION OF AN EARTHQUAKE DATABASE
"""

from geoist.cattools import Catalogue as Cat
Db = Cat.Database('First Test','NOTE: Just a test')

L = [{'Year': 1960},
     {'Year': 1961, 'Month': 12, 'Latitude': 10., 'Longitude': 20.},
     {'Year': 1962, 'Month': 12, 'Day': 3, 'Hour': 5, 'Minute': 20, 'Second': 10}]

M = [{'MagCode': 'AAA', 'MagSize':5, 'MagError': 0.1, 'MagType':'Mw'},
     {'MagCode': 'XXX', 'MagSize':7, 'MagError': 0.2, 'MagType':'ML'}]

# Creating an new empty catalogue item
Db.AddEvent('E001')

# Creating a new item with just Location information
Db.AddEvent('E002', Location=L)

# Creating a new item with Location and Magnitude information
Db.AddEvent('E003', Location=L, Magnitude=M)

# Adding new information to an existing item
Db.AddEvent('E001', L, [], Append=True)
Db.AddEvent('E002', Magnitude=M, Append=True)

# Remove an existing item (by ID)
Db.DelEvent('E003')

# Visualize item information
Db.Print('E002')
