#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

#-----------------------------------------------------------------------------------------

class GeoMap:
  '''
  INFO:

  Map boundary edges order:
  [LeftLowerLon,LeftLowerLat,UpperRightLon,UpperRightLat]

  Background type:
    'none'
    'etopo'
    'esri' --> background source

  Background sources available for 'esri':
    ESRI_Imagery_World_2D (MapServer)
    ESRI_StreetMap_World_2D (MapServer)
    I3_Imagery_Prime_World (GlobeServer)
    NASA_CloudCover_World (GlobeServer)
    NatGeo_World_Map (MapServer)
    NGS_Topo_US_2D (MapServer)
    Ocean_Basemap (MapServer)
    USA_Topo_Maps (MapServer)
    World_Imagery (MapServer)
    World_Physical_Map (MapServer)
    World_Shaded_Relief (MapServer)
    World_Street_Map (MapServer)
    World_Terrain_Base (MapServer)
    World_Topo_Map (MapServer)
  '''

  #---------------------------------------------------------------------------------------

  def __init__(self, Cfg=[]):

    # Defaults (example setting)
    if not Cfg:
      self._cfg = {'Bounds': [7., 36., 19., 48.],
                   'FigSize': [6., 6.],
                   'Background': ['esri','World_Terrain_Base',1500],
                   'Grid': [5., 5.]}
    else:
      self._cfg = Cfg

    self._zo = 1

  #---------------------------------------------------------------------------------------

  def BasePlot(self):

    plt.figure(figsize = (self._cfg['FigSize'][0],
                          self._cfg['FigSize'][1]))

    # Basemap
    self._map = Basemap(self._cfg['Bounds'][0],
                        self._cfg['Bounds'][1],
                        self._cfg['Bounds'][2],
                        self._cfg['Bounds'][3],
                        resolution = 'l',
                        projection = 'tmerc',
                        epsg = 3857)

    # Background land
    if self._cfg['Background'][0] == 'color':
      self._map.drawlsmask(land_color = self._cfg['Background'][1],
                           ocean_color = self._cfg['Background'][2],
                           grid = 1.25,
                           lakes = True)

    if self._cfg['Background'][0] == 'etopo':
      self._map.etopo(zorder = self._zo)

    if self._cfg['Background'][0] == 'esri':
      self._map.arcgisimage(service = self._cfg['Background'][1],
                            xpixels = self._cfg['Background'][2],
                            dpi = 300,
                            zorder = self._zo)

    if self._cfg['Background'][0] == 'relief':
      self._map.shadedrelief()

  #---------------------------------------------------------------------------------------

  def DrawGrid(self):

    # Parallels and meridians
    parallels = np.arange(-90, 90, self._cfg['Grid'][1])
    meridians = np.arange(0, 360., self._cfg['Grid'][0])

    self._zo += 1
    self._map.drawparallels(parallels, labels = [1,0,0,0],
                            fontsize = 14, weight = 'normal',
                            linewidth = 0.5,
                            zorder = self._zo)
    self._zo += 1
    self._map.drawmeridians(meridians, labels = [0,0,0,1],
                            fontsize = 14, weight = 'normal',
                            linewidth = 0.5,
                            zorder = self._zo)

  #---------------------------------------------------------------------------------------

  def DrawBounds(self):

    # Boundaries and lines
    self._zo += 1
    self._map.drawcoastlines(linewidth = 0.8,
                             zorder = self._zo)
    self._zo += 1
    self._map.drawstates(linewidth = 0.8,
                         zorder = self._zo)
    self._zo += 1
    self._map.drawcountries(linewidth = 0.8,
                            zorder = self._zo)
    self._zo += 1
    self._map.drawrivers(linewidth = 0.1,
                         color = 'b',
                         zorder = self._zo)
    """
    self._zo += 1
    self._map.drawmapboundary(linewidth = 2,
                              color = 'k',
                              zorder = self._zo)
    """

  #---------------------------------------------------------------------------------------

  def Title(self, string, Set=['bold','k',18]):

    plt.title(string, weight = Set[0],
                      color = Set[1],
                      fontsize = Set[2])

  #---------------------------------------------------------------------------------------

  def PointPlot(self, Lon, Lat, Label=[], Set=['o','y',5,1]):

    x, y = self._map(Lon, Lat)

    self._zo += 1
    self._map.plot(x, y, Set[0],
                         color = Set[1],
                         markersize = Set[2],
                         markeredgewidth = Set[3],
                         label = Label,
                         zorder = self._zo)

  #---------------------------------------------------------------------------------------

  def LabelPlot(self, Lon, Lat, Label, Set=['normal','k',14]):

    x, y = self._map(Lon, Lat)

    # If only one label provided, convert to list
    if isinstance(Label, str):
      x = [x]
      y = [y]
      Label = [Label]

    self._zo += 1
    for i, string in enumerate(Label):
      plt.text(x[i], y[i], string, weight = Set[0],
                                   color = Set[1],
                                   fontsize = Set[2],
                                   zorder = self._zo)

  #---------------------------------------------------------------------------------------

  def AreaPlot(self, Lon, Lat, Set=['y',1,'k',1]):

    x, y = self._map(Lon, Lat)

    if Set[0]:
      self._zo += 1
      plt.fill(x, y, color = Set[0],
                     alpha = Set[1],
                     zorder = self._zo)
    if Set[2]:
      self._zo += 1
      plt.plot(x, y, Set[2],
                     linewidth = Set[3],
                     zorder = self._zo)

  #---------------------------------------------------------------------------------------

  def MeshPlot(self, Lon, Lat, Elev, Cmap=[], Clim=[], Mesh=True):

    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator

    x, y = self._map(Lon, Lat)
    z = Elev

    if not Cmap:
      Cmap = cm.jet
    # cmap.set_under('w', alpha=0.)

    if not Clim:
      Clim = [z.min(), z.max()]

    levels = MaxNLocator(nbins=16).tick_values(Clim[0], Clim[1])
    norm = BoundaryNorm(levels, ncolors = Cmap.N, clip=True)

    if not Mesh:
      self._zo += 1
      h = plt.scatter(x, y, c = z,
                            s = 20,
                            marker = 's',
                            cmap = Cmap,
                            vmin = Clim[0],
                            vmax = Clim[1],
                            lw = 0,
                            alpha = 1.,
                            zorder = self._zo)

    else:
      self._zo += 1
      z = z[:-1, :-1]
      h = plt.pcolormesh(x, y, z, cmap = Cmap,
                                  norm = norm,
                                  vmin = Clim[0],
                                  vmax = Clim[1],
                                  lw = 0,
                                  alpha = 1.,
                                  zorder = self._zo)

    clb = plt.gcf().colorbar(h, orientation = 'vertical')
    clb.outline.set_linewidth(1)
    clb.ax.tick_params(labelsize=14)
    clb.set_label('Spectral Acceleration ($g$)', size=12)

  #---------------------------------------------------------------------------------------

  def ShapeFile(self, ShpFile, Name, Color='k'):

    # NOTE: this doesn't always work with polygons,
    # better to use the function in crd_tool

    self._zo += 1
    self._map.readshapefile(ShpFile,
                            Name,
                            linewidth = 1.5,
                            drawbounds = True,
                            color = Color,
                            zorder = self._zo)

  #---------------------------------------------------------------------------------------

  def Legend(self, Location=[]):

    self._zo += 1
    if Location:
      l = plt.legend(loc = Location, numpoints = 1)
    else:
      # Default outside
      l = plt.legend(bbox_to_anchor = (1.05, 1), 
                     loc = 2, borderaxespad = 0.,
                     numpoints = 1)
    l.set_zorder(self._zo)

  #---------------------------------------------------------------------------------------

  def Show(self):
 
    plt.show(block = False)

  #---------------------------------------------------------------------------------------

  def Close(self):

    plt.close('all')

  #---------------------------------------------------------------------------------------

  def SaveFig(self, OutFile, Dpi=150):

    plt.savefig(OutFile, bbox_inches = 'tight', dpi = Dpi)
