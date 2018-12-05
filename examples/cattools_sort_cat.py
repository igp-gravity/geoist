"""
EXAMPLE 8 - Sort catalogue
"""
from os.path import dirname
from geoist.cattools import Catalogue as Cat
from geoist.cattools import Exploration as Exp

#-----------------------------------------------------------------------------------------
# Import catalogue
pathname = dirname(__file__)
Db1 = Cat.Database()
Db1.Load(pathname+'/data/isc-rev-africa-select.bin')

#-----------------------------------------------------------------------------------------
# Randomly shuffling catalogue (for testing)

import numpy as np
Ind = np.random.randint(0,Db1.Size(),Db1.Size())

Db2 = Cat.Database('Unsorted')
for I in Ind:
  Db2.Events.append(Db1.Events[I])

#-----------------------------------------------------------------------------------------
# Sorting again

Db2.Sort()
print('Cat Sort!')
