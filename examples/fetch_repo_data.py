# coding: utf-8
# download data from Github repo 

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from geoist.others import fetch_data

data = fetch_data.fetch_gra_data()  #dataframe

lon = data["Longitude"].values
lat = data["Latitude"].values
bathymetry = data["Elevation"].values

# Make a plot of the decimated data using Cartopy
plt.figure()
plt.axis('scaled')
ax = plt.axes(projection=ccrs.Mercator())
ax.set_title("Block Median Bathymetry")
# Plot the bathymetry as colored circles.
plt.scatter(lon, lat, c=bathymetry, s=5, transform=ccrs.PlateCarree())
plt.colorbar().set_label("meters")
# Use a utility function to setup the tick labels and land feature
fetch_data.setup_cn_bathymetry_map(ax)
plt.tight_layout()
plt.show()
