"""
EXAMPLE 5 - CATALOGUE MERGING
"""
from os.path import dirname
from geoist.cattools import Catalogue as Cat
from geoist.cattools import Selection as Sel
from geoist.cattools import MagRules as MR

#-----------------------------------------------------------------------------------------
# Import Catalogues
pathname = dirname(__file__)
Db1 = Cat.Database()
Db1.Load(pathname+'/data/isc-rev-africa-select.bin')

Db2 = Cat.Database()
Db2.Load(pathname+'/data/isc-gem-v3.bin')

#-----------------------------------------------------------------------------------------
# Duplicate findings

# Between Catalogues
Db3, Log = Sel.MergeDuplicate(Db1,Db2,Twin=60.,Swin=50.,Log=1, Owrite=False)

# Within a catalogue
Log = Sel.MergeDuplicate(Db1,Twin=60.,Swin=50.,Log=1)

#-----------------------------------------------------------------------------------------
# Magnitude conversion

# Apply to all agency
Sel.MagConvert(Db1,'*',['Ms','MS'],'Mw',MR.Ms_Mw_Lin_DiGiacomo2015)

# Apply to single agency
Sel.MagConvert(Db1,'ISC','Ms','Mw',MR.Ms_Mw_Exp_DiGiacomo2015)
