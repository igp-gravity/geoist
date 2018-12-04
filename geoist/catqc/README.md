Catalog Statistics and Comparison
=================================

Generate plots and statistics for single catalogs and for catalog-catalog comparisons. `QCreport.py` generates a catalog summary and figures concerning various data within the catalog. It is run by typing `QCreport.py <catalog> <startyear> <endyear>`; e.g. `QCreport.py us 2010 2012`. The resulting figures and data can be found in the generated folder with the format `<catalog><startyear>-<endyear>`; e.g. `us2010-2012`.

`QCmulti.py` generates figures comparing various data with two different catalogs. It is run by typing `QCmulti.py <catalog1> <catalog2> <startyear> <endyear>`; e.g. `QCmulti.py nc us 2010 2012`. The resulting figures and data can be found in the generated folder with the format `<catalog1>-<catalog1>_<startyear>-<endyear>`; e.g. `nc-us_2010-2012`.

Tested and working in Python 2.7 and 3.6.

Required Python packages
------------------------
1. [NumPy](http://www.numpy.org)
2. [Matplotlib](https://matplotlib.org)
3. [Pandas](http://pandas.pydata.org)
4. [Cartopy](http://scitools.org.uk/cartopy)
5. [ObsPy](https://www.obspy.org)
6. [Markdown](https://python-markdown.github.io)
