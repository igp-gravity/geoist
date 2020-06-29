# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 18:25:58 2019

@author: chens
"""

from tkinter import *
import os
#from simpledialog import simpledialog
from geoist.vis import gimodule

gimodule.maxwidth  = 140
# Since the interface now has to columns of buttons this must be wider

## Constants
programname = "MagTools APIs - Geomagnetic reference field models"
version = "0.1"

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

        mainLabel = gimodule.mainline(frame, "定向钻探(NWD)参考地磁场模型计算接口程序", version)
        
        gimodule.seperator_line(frame,10, "获得地磁场信息")
        covfit = gimodule.LauncherButton(frame,11,18,0,"最近地磁台",
        lambda:[covfit.run("covfit.py")], "API: http://0.0.0.0/magv1/nearestMagSta?")
        empcov = gimodule.LauncherButton(frame,11,18,1,"指定范围地磁台", lambda:
        [empcov.run("empcov.py")],	"API: http://0.0.0.0/magv1/selMagSta?")
        geocol	= gimodule.LauncherButton(frame,12,18,0,"地磁分量转换", lambda:
        [geocol.run("geocol.py")], "API: xyz2hdi/hdi2xyz")
        tc  = gimodule.LauncherButton(frame,12,18,1,"模型解算最优时变", lambda:
        [tc.run("tc.py")], "API: magts")

        gimodule.seperator_line(frame,40, "主磁场+岩石圈磁场：EMM2015地磁模型接口")    
        geogrid = gimodule.LauncherButton(frame,41,18,0,"模型解算单点",
        lambda: [geogrid.run("geogrid.py")], "API: emmpnt ")
        geoip = gimodule.LauncherButton(frame,41,18,1,"模型解算网格", lambda:
        [geoip.run("geoip.py")], "API: emmgrd")
        geoegm = gimodule.LauncherButton(frame,42,18,0,"模型解算时间序列",
        lambda:[geoegm.run("geoegm.py")], "API: emmts")
        stokes = gimodule.LauncherButton(frame,42,18,1,"模型解算多点", lambda:
        [stokes.run("stokes.py")], "API: emmpnts")

        gimodule.seperator_line(frame,70, "主磁场1：IGRF12地磁模型接口")    
        geogrid = gimodule.LauncherButton(frame,71,18,0,"模型解算单点",
        lambda: [geogrid.run("geogrid.py")], "API: igrfpnt ")
        geoip = gimodule.LauncherButton(frame,71,18,1,"模型解算网格", lambda:
        [geoip.run("geoip.py")], "API: igrfgrd")
        geoegm = gimodule.LauncherButton(frame,72,18,0,"模型解算时间序列",
        lambda:[geoegm.run("geoegm.py")], "API: igrfts")
        stokes = gimodule.LauncherButton(frame,72,18,1,"模型解算多点", lambda:
        [stokes.run("stokes.py")], "API: igrfpnts")

        gimodule.seperator_line(frame,100, "主磁场2：WMM2015地磁模型接口")    
        geogrid = gimodule.LauncherButton(frame,101,18,0,"模型解算单点",
        lambda: [geogrid.run("geogrid.py")], "API: wmmpnt ")
        geoip = gimodule.LauncherButton(frame,101,18,1,"模型解算网格", lambda:
        [geoip.run("geoip.py")], "API: wmmgrd")
        geoegm = gimodule.LauncherButton(frame,102,18,0,"模型解算时间序列",
        lambda:[geoegm.run("geoegm.py")], "API: wmmts")
        stokes = gimodule.LauncherButton(frame,102,18,1,"模型解算多点", lambda:
        [stokes.run("stokes.py")], "API: wmmpnts")

        gimodule.seperator_line(frame,130, "电离层磁场：DIFI-4地磁模型接口")    
        geogrid = gimodule.LauncherButton(frame,131,18,0,"模型解算单点",
        lambda: [geogrid.run("geogrid.py")], "API: difipnt ")
        geoip = gimodule.LauncherButton(frame,131,18,1,"模型解算网格", lambda:
        [geoip.run("geoip.py")], "API: difigrd")
        geoegm = gimodule.LauncherButton(frame,132,18,0,"模型解算时间序列",
        lambda:[geoegm.run("geoegm.py")], "API: dififts")
        stokes = gimodule.LauncherButton(frame,132,18,1,"模型解算多点", lambda:
        [stokes.run("stokes.py")], "API: difipnts")                        
            
        gimodule.seperator_line(frame,gimodule.maxrow-2)
        button = Button(frame, text="退出", width=8, command=frame.quit)
        button.grid(row=gimodule.maxrow, column=0, sticky=W)
######################################################
## Initiate the program and start program loop
######################################################

root = Tk()
app = App(root)
root.title(programname)
root.mainloop()