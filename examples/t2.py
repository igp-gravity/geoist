# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:03:33 2018

@author: chens
"""
import struct
import numpy as np
from osgeo import gdal
from matplotlib import pyplot as plt

class GrdData(object):
    """
    Grid Data Object
    Attributes
    ----------
    data : numpy masked array
        array to contain raster data
    tlx : float
        Top Left X coordinate of raster grid
    tly : float
        Top Left Y coordinate of raster grid
    xdim : float
        x-dimension of grid cell
    ydim : float
        y-dimension of grid cell
    nrofbands : int
        number of raster bands
    dataid : str
        band name or id
    rows : int
        number of rows for each raster grid/band
    cols : int
        number of columns for each raster grid/band
    nullvalue : float
        grid null or nodata value
    norm : dictionary
        normalized data
    gtr : tuple
        projection information
    wkt : str
        projection information
    units : str
        description of units to be used with color bars
    """
    def __init__(self):
        self.data = np.ma.array([])
        self.tlx = 0.0  # Top Left X coordinate
        self.tly = 0.0  # Top Left Y coordinate
        self.xdim = 1.0
        self.ydim = 1.0
        self.nrofbands = 1
        self.dataid = ''
        self.rows = -1
        self.cols = -1
        self.nullvalue = 1e+20
        self.norm = {}
        self.gtr = (0.0, 1.0, 0.0, 0.0, -1.0)
        self.wkt = ''
        self.units = ''
                 
def export_surfer(file_out, grd):
    """
    Export a surfer binary grid

    Parameters
    ----------
    data : grid Data
        dataset to export
    """
    fno = open(file_out, 'wb')
    xmin = grd.tlx
    xmax = grd.tlx + grd.cols*grd.xdim
    ymin = grd.tly - grd.rows*grd.ydim
    ymax = grd.tly
    bintmp = struct.pack('cccchhdddddd', b'D', b'S', b'B', b'B',
                         grd.cols, grd.rows,
                         xmin, xmax,
                         ymin, ymax,
                         np.min(grd.data),
                         np.max(grd.data))
    fno.write(bintmp)

    ntmp = 1.701410009187828e+38
    tmp = grd.data.astype('f')
    tmp = tmp.filled(ntmp)
    tmp = tmp[::-1]
    fno.write(tmp.tostring())

    fno.close()


def export_ascii(file_out, grd):
    """
    Export Ascii file

    Parameters
    ----------
    data : grid Data
        dataset to export
    """
    fno = open(file_out, 'w')

    xmin = grd.tlx
    ymin = grd.tly - grd.rows*grd.ydim

    fno.write("ncols \t\t\t" + str(grd.cols))
    fno.write("\nnrows \t\t\t" + str(grd.rows))
    fno.write("\nxllcorner \t\t\t" + str(xmin))
    fno.write("\nyllcorner \t\t\t" + str(ymin))
    fno.write("\ncellsize \t\t\t" + str(grd.xdim))
    fno.write("\nnodata_value \t\t" + str(grd.nullvalue))

    tmp = grd.data.filled(grd.nullvalue)

    for j in range(grd.rows):
        fno.write("\n")
        for i in range(grd.cols):
            fno.write(str(tmp[j, i]) + " ")

    fno.close()

if __name__ == "__main__":
    
 dat = [GrdData()]  #dat是一个列表对象存放Grid数据
 i = 0              #本例只处理一个grd数据，所以i=0即可

 #d1=gdal.Open("D:\demo\Demogrid.grd", gdal.GA_ReadOnly)
 d1=gdal.Open("D:\demo\GridFile.grd", gdal.GA_ReadOnly)
 rtmp=d1.GetRasterBand(1)
 r11=rtmp.ReadAsArray() 
 lonsdim = ((r11.max()-r11.min())/(r11.shape[1]-1))/2
 dat[i].data = r11    #把读入的数据存到GrdData的实例中
 nval = -9999.0       #已知的无数据特征值
 dat[i].data = np.ma.masked_equal(dat[i].data, nval)
 if dat[i].data.mask.size == 1:
   dat[i].data.mask = (np.ma.make_mask_none(dat[i].data.shape) +
                      dat[i].data.mask)

 dat[i].nrofbands = d1.RasterCount
 dat[i].tlx = 0
 dat[i].tly = 0
 dat[i].dataid = 'only for test'
 dat[i].nullvalue = nval
 dat[i].rows = r11.shape[0]
 dat[i].cols = r11.shape[1]
 dat[i].xdim = 1
 dat[i].ydim = 1 


 #r11[r11<-9000]=np.nan #如果有nodata的情况下，用这个函数
 r00=dat[0].data
 dat[0].data=r00*r00 #数据处理逻辑部分，r12为输出节点
 #保存本地文件
 export_ascii("D:\demo\GridFile1.grd", dat[0])  #保存后应该给下一个grd节点使用
 #export_surfer("D:\demo\GridFile2.grd", dat[0]) #二进制保存surfer格式grd文件
 
###画图输出用  
 plt.imshow(r00)     #显示绘图结果
 plt.imsave('D:\demo\GridFile.png',r00)  #保存绘图结果到本地，对接到运行模式的“显示运行窗口”模式下
