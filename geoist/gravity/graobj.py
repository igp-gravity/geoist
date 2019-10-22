# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 11:15:41 2019 at ISM, Tokyo

Author:  Shi CHEN   @ IGP-CEA
         Bei ZHANGE @ IGP-CEA

Co-Author: Jiancang ZHUANG @ ISM-Tokyo

####################################################
            MIT license
        Copyright @ pyBACGS 2019 
         All right reserved 
####################################################

Module which contains the main classes of the pyBACGS
            
CLASS list:
    Meter 
    Agmeter
    Rgmeter
    Station
    AGstation
    Network
    Chanlist
    Loop
    Survey
    Campaign
"""

import os, sys, glob, json, time
import numpy as np
import scipy.linalg as slg
from scipy.sparse import csr_matrix
from collections import OrderedDict
import geoist.pfm.tide as tide
from datetime import date, datetime, timedelta
import h5py
from functools import wraps
import geoist.gravity.adjmethods as adj

###############################################################################

def timeit(prefix):		# the params can be showed.
    def decorator(func):	
        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time.time()
            func(*args, **kwargs)
            end = time.time()
            #print(start, end)
            print('%s: duration time: %fs' % (prefix, end - start))
        return wrapper
    return decorator

# Denfined by pyBACGS 
class Meter(object): 
    """Gravimeter
    
    Properties:
     mtype : LCR-G/LCR-D/Burris/CG4/CG5/CG6, 
     msn : G147 
     muser : IGP-CEA
     msf :
     status : good/stable/unstable
     sftable : pd.DataFrame
                    
    functions:
      - read_table
      - export_table
      - save_meter
    export:
        m1 = Meter()
        s = json.dumps(m1.__dict__)
    """    
    
    #_scalefactor = 1.0
    count = 0
    
    def __init__(self, gtype, gsn):
        """
        """
        self._mtype = gtype
        self._msn = gsn
        self._scalefactor = 1.0
        self.table_data = []
        Meter.count += 1
        
    @property
    def mtype(self):
        return self._mtype
    
    @property
    def msn(self):
        return self._msn   

    @property
    def muser(self):
        return self._muser       
    
    @property
    def msf(self):
        return self._scalefactor
    
    @msf.setter
    def msf(self, value):
        if not isinstance(value, float):
            raise ValueError('scalefactor must be an float!')
        if value < 0.1 or value > 1.9:
            raise ValueError('scalefactor must between 0.1 ~ 1.9!')
        self._scalefactor = value
 
    @property
    def status(self):
        return self._mstatus   

    @property
    def sftable(self):
        return self._sftable
           
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        """
        print('%s type gravimeter, the ID is %s'%(self.mtype,self.msn))       
        return "====Meter Print by BACGS===="
    
    @staticmethod
    def showcount():
        print(Meter.count)

    def read_table(self, filename):
        """read the scalefactor table if LCR gravimeter
        """
        self.table_data = []
        try:
            f = open(filename, mode = 'r')
            if len(self.msn) > 0:
                print('The table of {} meter will load'.format(self.msn))
            else:
                print('Please set the SN of meter!')
                return
            i = 0
            iflag = False
            for line in f:
                i += 1
                line = line.strip()  # Clean line
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == '#'): continue
                if (line[0:5] == '99999') or (line[0:5] == '66666'): break
                if (line[0:5] == '88888'): iflag = False
                if iflag:
                    vals=line.split()
                    self.table_data.append(vals)
                if (line.lower() == self.msn.lower()):
                    iflag = True
            f.close
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('check raw data file')
        except IndexError:
            print('check raw data file: possibly last line?')
        #if self.msn
        print('The {} values has been read.'.format(len(self.table_data)))
        

    def export_table(self, filename):
        pass
    
    def save_meter(self, filename):
        try:
            f = open(filename, mode = 'w')
            f.write(json.dumps(self.__dict__))
            f.close
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('check raw data file')
        except IndexError:
            print('check raw data file: possibly last line?')  

class Agmeter(Meter):
    """
     Absolute gravimeter information
    """
    pass

class Rgmeter(Meter):
    """
     Relative gravimeter information
    """    
    pass
        

class Station(object):
    """
     Specific Station for repeated gravity observation
    """    
    def __init__(self, name, sid, stype, lon, lat, elev):
        """
        """
        self._stype = stype
        self._sid = sid
        self._sname = name
        self._slon = lon
        self._slat = lat
        self._elev = elev
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        """
        print('%s type station, the name is %s, the ID is %s'%(self._stype,self._sname,self._sid))       
        return "====Station Print by BACGS===="    
    @property
    def slon(self):
        return self._slon  
    @slon.setter
    def slon(self, value):
        if not isinstance(value, float):
            raise ValueError('ref gravity value must be an float!')
        if value < 70.0 or value > 135.0:
            raise ValueError('Longtitude must between E70 ~ 135!')           
        self._slon = value  
    @property
    def slat(self):
        return self._slat  
    @slat.setter
    def slat(self, value):
        if not isinstance(value, float):
            raise ValueError('ref gravity value must be an float!')
        if value < 10.0 or value > 60.0:
            raise ValueError('Longtitude must between N10 ~ 60!')  
        self._slat = value  
        
class AGstation(Station):
    """
     derived class of Station 
     Adding Properties:
       ref_gra
       ref_h_east_gd
       ref_north_gd
       ref_v_gd    
    """

    @property
    def ref_gra(self):
        return self._ref_gra  
    @ref_gra.setter
    def ref_gra(self, value):
        if not isinstance(value, float):
            raise ValueError('ref gravity value must be an float!')
        self._ref_gra = value
        
    @property
    def ref_gra_err(self):
        return self._ref_gra_err  
    @ref_gra_err.setter
    def ref_gra_err(self, value):
        if not isinstance(value, float):
            raise ValueError('ref gravity value must be an float!') 
        self._ref_gra_err = value     
        
    @property
    def ref_h_east_gd(self):
        return self._ref_h_east_gd 
    @ref_h_east_gd.setter
    def ref_h_east_gd(self, value):
        if not isinstance(value, float):
            raise ValueError('ref H East gravity gradient must be an float!')
        self._ref_h_east_gd = value 

        
    @property
    def ref_north_gd(self):
        return self._ref_north_gd  
    @ref_north_gd.setter
    def ref_north_gd(self, value):
        if not isinstance(value, float):
            raise ValueError('ref H North gravity gradient must be an float!')
        self._ref_north_gd = value

        
    @property
    def ref_v_gd(self):
        return self._ref_v_gd    
    @ref_v_gd.setter
    def ref_v_gd(self, value):
        if not isinstance(value, float):
            raise ValueError('ref V gravity gradient value must be an float!')
        self._ref_v_gd = value 

        
class Network(object):
    """Gravity Network Configure
    
    Properties:
                    
    functions:
    - 
                
    """     

    __type = {1: 'region', 2:'nationwide', 3: 'temporary',4: 'other'}
    
    def __init__(self, name, ntype):
        """
        """
        self.ntype=ntype
        self._name=name
        self.sta_list = []
        self.sta_table = []
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('%s type network %s, which include %d stations'%(self._ntype,self._name,len(self.sta_list)))       
        return "====Network Print by BACGS===="
    @property
    def ntype(self):
        return self._ntype 
    
    @ntype.setter
    def ntype(self, value):
        if isinstance(value, float):
            raise ValueError('Network type cannot be an float!')
        if isinstance(value, int):
            if value < 1 or value > 4:
                raise ValueError('Network type should be 1-4!')
            self._ntype = self.__type[value]
        elif isinstance(value, str):
            list(self.__type.values()).index(value.lower())
            self._ntype = value.lower() 
        else:
            raise ValueError('Network type wrong!!!')
        
    @property
    def name(self):
        return self._name 
        
    def add_station(self,sta):
        if not isinstance(sta, Station):
           raise ValueError('Input data type must be a Station!') 
        self.sta_list.append(sta)
    
    def read_pnts(self, filename):
        """
        read a generic dataset file and multi channels
        """                 
        try:
            fh = open(filename, 'r')  
            #print("number of lines: %d"%len([1 for line in open(filename, 'r')]))
            i = 0 
            for line in fh:    
                i += 1
                line = line.strip()  # Clean line
                # Skip blank and comment lines
                if (not line) or (line[0] == '/') or (line[0] == '#'): continue
                if (line[0:5] == '99999') or (line[0:5] == '66666'): break
                vals=line.split()
                #print(vals)
                # fill object properties:
                self.sta_table.append(vals)
                self.add_station(Station(vals[1],vals[0],1,float(vals[3]),
                                         float(vals[2]),float(vals[4])))                                                                                    
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('Reading at line %d : check raw data file'%(i))
        except IndexError:
            print('Reading at line %d : check raw data file: possibly last line?'%(i))       

class Chanlist(object):   
    # tuple
#    metername=[]
#    netname=[] 
#    obstime=[]
#    stationid=[]  
#    obsval=[]
#    airpress=[]
#    temp=[]
    def __init__(self):
        """
        """
        self.metername=[]
        self.netname=[]
        self.obstime=[]
        self.stationid=[]
        self.obsval=[]
        self.airpress=[]
        self.temp=[]
        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('length of channel %d'%(len(self.metername)))       
        return "====Channel print by BACGS===="    
    
class Loop(Chanlist):
    """
     Loop mean the whole close route. 
     In general, one day job can be seen a Loop.

    """
    meter_list = []    
    #tide_grav = []
    #airpress_grav = []
    station_dic = OrderedDict() #Dictionary that remembers insertion order
    
    def __init__(self, meter, starttime):
        """
        """
        super(Loop,self).__init__() # call properties from the baseclass
        self.meter = meter
        self.starttime =starttime
        self.tide_grav = []
        self.airpress_grav = []
        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('length of channel %d'%(len(self.metername)))  
        return "====Loop print by BACGS====" 
    
    def add_recording(self, vals):
         """
           add muli data channels
         """
         self.metername.append(vals[0])
         self.netname.append(vals[1]) 
         self.obstime.append(int(vals[2])) 
         self.stationid.append(vals[3])   
         self.obsval.append(float(vals[4])) 
         self.airpress.append(float(vals[5])) 
         self.temp.append(float(vals[6]))          
         self.tide_grav.append(0.0)
         self.airpress_grav.append(0.0)

class Survey(Chanlist):
    """
    Survey reading dataset from recording files.
    The instance of Survey class must be setted before send to Campaign class.
    % 55555 一个闭合结束, one loop
    % 99999 整个测量结束, one survey
    % 44444 66666 一天内中断
    """
    #loop_list = [] 
    #survey_table = []

    def __init__(self, name, time_tag):
        """
        """
        super(Survey,self).__init__() # call properties from the baseclass
        self._name = name
        self._time_tag = time_tag
        self.meter_list = []
        self.loop_list = [] 
        self.survey_table = []
        self.meter_sf_index = 0
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('Number of gravimeter used in Survey: %d'%(len(self.meter_list)))       
        print('Number of loops in Survey: %d'%(len(self.loop_list)))  
        print('Number of recordings can be used: %d'%(len(self.survey_table)))
        return "====Survey print by BACGS====" 
    
    def add_meter(self, meter):
        if not isinstance(meter, Meter):
            raise ValueError('Input is not a instance of Meter class!')        
        self.meter_list.append(meter)
    @property
    def name(self):
        return self._name 
    @property
    def time_tag(self):
        return self._time_tag    
    
    @property
    def net(self):
        return self._net 
    
    @net.setter
    def net(self, value):
        if not isinstance(value, Network):
            raise ValueError('Input is not a instance of Network class!')
        self._net = value
    
    @staticmethod
    def tran_datetime(obstime): #eg. 1304070850
        if not isinstance(obstime, int):
            raise ValueError('Datetime must be a integer type with 10 bit!')
    
        timeval=float(obstime);
        iy=int(timeval/1E8)
        im=int(timeval/1E6)-iy*100
        idd=int(timeval/1E4)-iy*10000-im*100
        ih=int(timeval/100)-iy*1000000-im*10000-idd*100
        imin=timeval-iy*100000000-im*1000000-idd*10000-ih*100
        hour_val=ih+imin/60
        
        if iy > 50: 
            iy=1900+iy  
        else:
            iy=2000+iy
            
        day_fraction = ih / 24.0 + imin / 1440.0 + 0.0 / 86400.
        days_val = date(iy, im, idd).toordinal() #+ day_fraction
        #print(days_val)
        time_list = [iy, im, idd, ih, imin, days_val, hour_val]
        
        return time_list
    
    def _get_pnt_loc(self, sid):
        """
        Get location(lon, lat) from the Network instance by the ID of station
        """
        lon = 0.0
        lat = 0.0
        for snet in self.net.sta_list:
            if (snet._sid == sid):
                lon = snet.slon
                lat = snet.slat
                elev = snet._elev
                break
        return lon, lat, elev

    def _get_meter_sf(self, mname):
        """
        Get scalefactor of meter 
        """
        for sname in self.meter_list:
            if (sname.msn == mname):
                sfval = sname._scalefactor
                break
        return sfval
    
    def corr_aux_effect(self, alpha = 1.16, beta = -0.3):
        """Get gravity by table of LCG type meter
        Get earthtide and atomsphere effect 
        """
        j = 0
        for slp in self.loop_list:
            i = 0
            gra_list = []
            for sname in self.meter_list:
                if (sname.msn == slp.metername[0]):
                    if (len(sname.table_data) > 1):
                        table_list = sname.table_data
                        gra_list = self.read2mgal(slp.obsval, table_list)
                        break
             #Just for LCG meter
            for stime in slp.obstime:
                tt = self.tran_datetime(int(stime))
                #gdate0 = datetime(int(tt[0]), int(tt[1]), int(tt[2])) #20191013
                slocid = slp.stationid[i]
                slon, slat, selev = self._get_pnt_loc(slocid)
                gdate = datetime(int(tt[0]), int(tt[1]), int(tt[2]), 
                                 int(tt[3]), int(tt[4]), 00)
                gdate = gdate - timedelta(hours=8)
                g1 = tide.TideModel()
                gm, gs, g =g1.solve_longman(slat,slon,selev,gdate)
                slp.tide_grav[i] = -g #*alpha -1.0 ****
                slp.airpress_grav[i] = slp.airpress[i]*beta*0.001   
                sfval = self._get_meter_sf(slp.metername[i])
                if (len(gra_list) == 0):
                    slist = [slp.metername[i], slocid, tt[5], tt[6], slp.obsval[i], sfval,
                         slp.tide_grav[i], slp.airpress_grav[i], i+1]
                else:
                    slist = [slp.metername[i], slocid, tt[5], tt[6], gra_list[i], sfval,
                         slp.tide_grav[i], slp.airpress_grav[i], i+1]
                self.survey_table.append(slist)
                i += 1  
            j += 1
            
        return

    @staticmethod
    def read2mgal(read_list, table_list):
        """Transform the reading to gravity value using table_data of meter.      
        read2mgal is staticmethod
        Method : x1val*(table(j,2)+(dtd(i)-table(j,1))*table(j,3))
        Returns
        -------
        gra_list: the gravity values

        """
        gra_list = []
        if len(table_list) > 0:
            for reading in read_list:
                i = 0
                for table in table_list[:-1]:
                    i += 1
                    tt = [float(f) for f in table]
                    tt1 = [float(f) for f in table_list[i]]
                    if (reading >= tt[0]) and (reading < tt1[0]):
                        gval = tt[1] + (reading - tt[0])*tt[2]
                        gra_list.append(gval)
        return gra_list

    def read_survey_file(self, filename):
        if not isinstance(self.net, Network):
            raise ValueError('Property net must be set!')  
        if len(self.meter_list) < 1:
            raise ValueError('Property meter must more than 1 set!')              
        
        try:
            fh = open(filename, 'r')  
            #print("number of lines: %d"%len([1 for line in open(filename, 'r')]))
            i = 0 
            flag = 0
            for line in fh:    
                i += 1
                line = line.strip()  # Clean line
                # Skip blank and comment lines
                if (not line) or (line[0].lower() == 'a') or (line[0] == '#'): 
                   i -= 1
                   continue   

                vals=line.split()
                 
                if (i == 1): l1tmp = Loop(vals[0],vals[2])
                
                #print(l1tmp) 
                if (flag == 1): 
                    l1tmp = Loop(vals[0],vals[2])
                    flag = 0
                    
                if (line[0:5] == '55555') or (line[0:5] == '33333') or \
                (line[0:5] == '44444') or (line[0:5] == '66666'):
                    self.loop_list.append(l1tmp) #loop
                    flag = 1
                elif (line[0:5] == '99999'): 
                    self.loop_list.append(l1tmp)
                    break #loop
                else:    
                    l1tmp.add_recording(vals)

        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('Reading at line %d : check raw data file'%(i))
        except IndexError:
            print('Reading at line %d : check raw data file: possibly last line?'%(i))       
           

class Campaign(object):
    """
    经典单期平差，贝叶斯平差， 动态平差，混合平差，格值估计，
    type = 1.单期; 2.动态; 3.混合;
    eg: gravwork = Campaign('IGP201702',type = 1)        #1.初始化平差对象
        
        m1 = Meter('G147')                      #2.初始化重力仪
        m1.read_table('table_name')
        m1.set_scale = 1.000002
        
        n1 = Network('Huabei')                  #3.读取重力测网信息（点坐标，类型等）
        n1.read_netinfo('datsfilename')
        
        s1 = Survey('survey_name1','time_tag1') #4.初始化一次测量
        s1.set_net(n1)
        s1.add_meters(meter_List)
        s1.read_raw_data('datsfilename')
        s1.corr_earthtide(alpha = 1.16)
        s1.corr_atmosphere(beta = -0.3)
        print(s1)                               #输出测量基本信息
        
        s2 = Survey('survey_name2','time_tag2')
        ...
        
        stagrav1 = Station('BJT01')             #台站校正后的日均值数据-混合平差
        stagrav1.set_meter = mxx
        stagrav1.set_tide_free_data = dataframe1 
        
        ag = AGstation()
        
        gravwork.add_ag_sta(ag)                  #添加绝对点信息
        gravwork.add_surveys(survey_list)        #添加测量到平差任务
        gravwork.adj_method =  1                 #平差方法选择
        gravwork.pre_adj                         #准备平差矩阵（optional）
        
        gravwork.run_adj                         #5.运行平差
        
        gravwork.export_station()                #输出点值结果
        gravwork.export_drift()                  #输出漂移结果
        gravwork.export_error()
        
    """
        
    survey_dic = {}
    _adj_method = 1
    __adj_method = {1: 'cls', 2:'bay', 3: 'bay1',4: 'bay2'}
    #survey_list = []
    #agstation_list = []
    _mat_list = []
    _obs_list = []
    def __init__(self, name, camp_type):
        """
        """
        self._name = name
        self._type = camp_type
        self.survey_list = []
        self.agstation_list = [] 
        self._adj_method = 1
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('Number of Survey used in Campaign: %d'%(len(self.survey_list)))       
        return "====Campaign print by BACGS====" 
    
    @property
    def adj_method(self):
        return self._adj_method 
    
    @adj_method.setter
    def adj_method(self, value):
        if isinstance(value, float):
            raise ValueError('Adjustment method cannot be an float!')
        if isinstance(value, int):
            if value < 1 or value > 4:
                raise ValueError('Adjustment method should be 1-4!')
            self._adj_method = self.__adj_method[value]
        elif isinstance(value, str):
            list(self.__adj_method.values()).index(value.lower())
            self._adj_method = value.lower() 
        else:
            raise ValueError('Adjustment method wrong!!!')
    @property
    def mat_list(self):
        return self._mat_list 
    @property
    def obs_list(self):
        return self._obs_list 
            
    def add_ag_sta(self, ag):
        """add the absolute gravity station to the Campaign class
        """
        if not isinstance(ag, AGstation):
            raise ValueError('Input is not a instance of AGstation class!')        
        self.agstation_list.append(ag)

    def add_surveys(self, survey):
        """
        """
        if not isinstance(survey, Survey):
            raise ValueError('Input is not a instance of Survey class!')        
        self.survey_list.append(survey)
        
    def pre_adj(self):
        """
         checking survey information
         generated matrix 
        """
        try:
            #sur_tab = self.survey_list[0].survey_table #only one survey can be used
            self._mat_list, self._obs_list = self.survey_matrix(0)
            #print(_mat_list)
        except:
            return False
            print('Please check the survey data file')            
            
        return True
    
    @staticmethod
    def tran_obs_array(rec_vals, pnt_id, kt):
        """
        Matrix data type used Numpy ndarray. eg: np.ndarray([3,2])

        Matlab old version code likes:  
            [Ut,Ur,dcnos,dcnoe,hoursdd1,ett1,prr1,DYn]=
            tra2array(point_no21,point_no_ord,day_num1,data_num1,
            g1_val,et1_val,hour1_val,press_val1,kt2);
            
        (pnt_ids, pnt_id0, day_num, data_num, g1_val, et_val, 
                       hour_val, press_val, kt):
        """
        
        mlen = len(pnt_id) 
        nlen = len(rec_vals)
        pnt_no = get_2d_list_slice(rec_vals, 0, nlen,1,2)
        idx_num = get_2d_list_slice(rec_vals, 0, nlen,8,9)
        
        g_val = []
        et_val = []
        hour_val = []
        press_val = []
        day_num = []
        loop_num = []

        for reci in rec_vals:
            hour_val.append(reci[3])
            g_val.append(reci[4]*reci[5])
            et_val.append(reci[6])
            press_val.append(reci[7])
            day_num.append(reci[2])
            loop_num.append(reci[8])
        iprev = 1

        data_num = []
        for ii in idx_num[1:]:
            itmp = ii[0] 
            if (itmp > iprev): 
                iprev = itmp
            else:
                data_num.append(iprev)
                iprev = 0
        data_num.append(ii[0])
        #print(data_num) 
        Um3 = np.zeros([int(kt[0]),mlen])   
        Um2 = np.zeros([int(kt[0]),2])
        DYn = np.zeros(int(kt[0]))
        hoursdd1 = np.zeros([int(kt[0]),2])
        dcnos = []
        dcnoe = []
        k = 0
        #print(kt)
        
        for i in range(int(kt[1])):
            jtval=data_num[i]
            jt1=sum(data_num[0:i+1])
            #print('dd',i, jtval, jt1)
            for j in range(data_num[i]-1):
                jt2=jt1-jtval+j
                jt0=jt1-jtval+1            
                if (loop_num[jt2] < jtval):
                    #print(jt2, jt0)
                    DYn[k] = g_val[jt2+1] - g_val[jt2]
                    #print( g_val[jt2+1] , g_val[jt2])
                    ett1 = et_val[jt2+1] - et_val[jt2]
                    prr1 = press_val[jt2+1] - press_val[jt2]
                    hoursdd1[k,0] = hour_val[jt2] - hour_val[jt0]
                    hoursdd1[k,1] = hour_val[jt2+1] - hour_val[jt2]
                    ktemp = pnt_no[jt2]
                    jkk = pnt_id.index(ktemp[0])
                    #jkk = np.where(pnt_id == ktemp) #jkk=find(pnt_id == ktemp)
                    Um3[k,jkk]=Um3[k,jkk]-1
                    ktemp=pnt_no[jt2+1]
                    jkk = pnt_id.index(ktemp[0])
                    #print(jkk, ktemp[0]) 
                    #jkk = np.where(pnt_id == ktemp) #jkk=find(pnt_id==ktemp)        
                    Um3[k,jkk] = Um3[k,jkk]+1 
                    Um2[k,0] = ett1
                    Um2[k,1] = prr1
                    dcnos.append(pnt_no[jt2])
                    dcnoe.append(pnt_no[jt2+1])
                    k=k+1    
       # print(DYn) #Um2/Um3/hrs/DYn/data_num/day_num
        #print('dddddd' ,k)
        obs_mat = []
        obs_mat.append(Um2)
        obs_mat.append(Um3.astype(int))
        obs_mat.append(hoursdd1)
        obs_mat.append(DYn)
        obs_mat.append(data_num)
        obs_mat.append(loop_num)
        obs_mat.append(day_num)

        return obs_mat, dcnos, dcnoe
    
        
    def survey_matrix(self, index = 0):
        """
        Obs matrix likes as follows:
                   |Ud1  0      Ur1| |Xd1|  |Db1|  
              A X= |0    Ud2    Ur2| |Xd2|= |Db2|  
                   |Us1   0     0  | |XX|   |0| 
                   |0    Us2    0  |        |0|         
        """
        #obs_mat  = np.ndarray([2, 3])
        #time_mat = np.ndarray([2, 3])
        # Get all meter obs station list
        #get_2d_list_slice()
        #self.survey_list[0].survey_table
        #pnt_ids = ##
        obs_len = len(self.survey_list[index].survey_table)
        #kt1 =  #segment number
        #kt2 =  #days
        num_meter = len(self.survey_list[index].meter_list)
        #num_loop = len(self.survey_list[index].loop_list)
        loop_meter = []
        loop_time = []
        loop_obsmeter = []
        for ll in self.survey_list[index].loop_list:
            #print([int(r1/10000.) for r1 in ll.obstime])
            loop_time.append([int(r1/10000.) for r1 in ll.obstime])
            loop_obsmeter.append(ll.metername)
            loop_meter.append(ll.meter)

        loop_time1 = list(flatten(loop_time))
        loop_obsmeter1 = list(flatten(loop_obsmeter))
        #print('bbbbb',len(loop_time1))
        #print(len(loop_obsmeter1),len(loop_time1))
        #print(loop_obsmeter1)
        gravlen = np.zeros([num_meter,3])
        for ii in range(num_meter):
          loop_time2 = []
          mname = self.survey_list[index].meter_list[ii].msn
          meterx = get_2d_list_slice(self.survey_list[index].survey_table, 0,obs_len,0,1)
          #daysx = get_2d_list_slice(self.survey_list[index].survey_table, 0,obs_len,8,9)
          kt1 = meterx.count([self.survey_list[index].meter_list[ii].msn])
          #print(loop_meter)
          #print(loop_time)
          kt2 = loop_meter.count(self.survey_list[index].meter_list[ii].msn)
          gravlen[ii,0] = int(kt1 - kt2)
          for ll, tt in zip(loop_obsmeter1, loop_time1):
              if (ll.lower() == mname.lower()):
                  loop_time2.append(tt)

          kt3 = len(set(loop_time2)) #20191014
          gravlen[ii,1] = int(kt2)
          gravlen[ii,2] = int(kt3)
         #print(kt1,kt2,kt3)
        #print(gravlen)
        self._gravlen = gravlen.astype(int) #save to object
        ls_all = get_2d_list_slice(self.survey_list[index].survey_table, 0,obs_len,1,2)
        pnt_ids = []
        for ll in ls_all: pnt_ids.append(ll[0])
        pnt_id0 = list(set(pnt_ids)) #去重复点号
        pnt_id0.sort()
        #print(pnt_id0)
        print("the %d stations will be computed."%len(pnt_id0))
        
        m = len(self.agstation_list)
        n1 = len(pnt_id0) 
        n2 = sum(self._gravlen[:,2]) #20191014
        uag = np.zeros([m, n1 + n2], dtype= int) 
        dag = np.zeros(m, dtype= float)
        wag = np.zeros(m, dtype= float)
        #print(m, n1, n2)
        k_ind = 0
        #print('99999')
        for ii in range(m):
            k_ind = -1
            #print(self.agstation_list[ii]._sid)
            #print(pnt_id0)
            #print(k_ind)
            try:
                k_ind = pnt_id0.index(self.agstation_list[ii]._sid)
            except ValueError:
                print('the AG station {} information cannot be find in Survey'.format(self.agstation_list[ii]._sid))
            #print(self.agstation_list[ii]._sid)
            dag[ii] = self.agstation_list[ii]._ref_gra
            wag[ii] = self.agstation_list[ii]._ref_gra_err**2 #sigma

            #if (k_ind < 0):
            #    print('Error, the AG station information cannot be find in Survey')
            uag[ii,n2+k_ind]=1
        #print(uag)   
         #to-do list
        dcs = []
        dce = []

        for ii in range(num_meter):   
            
            #mat_list = self.tran_obs_array() 
            meter_name = self.survey_list[index].meter_list[ii].msn
            print(meter_name)
            #print('ddd',self.adj_method)
            if (self.adj_method.lower() == 'bay1'):
                meter_sf =self.survey_list[index].meter_list[ii].msf
                print(meter_sf)
            recs = [] #recording data by gravimeter
            for jj in self.survey_list[index].survey_table:
                if (jj[0] == meter_name):
                    recs.append(jj)
            #print(recs[0])
            #print(recs[2])
            obs_mat, dcnos, dcnoe = self.tran_obs_array(recs, pnt_id0, gravlen[ii,:]) #step 1   
            #obs_mat : Um2/Um3/hrs/DYn/data_num/loop_num/day_num
            maxhrs = 16.0 #the max working duration per day.
            isgood, hourmax = self.check_hs(obs_mat[2], maxhrs) 
            
            if not isgood:
                print('the max hour is %f, which is too long in one day!'%hourmax)
            else:
                print('the max work hour in the survey is %f '%hourmax)
            ud = self.gen_drift(obs_mat[2], obs_mat[4], obs_mat[6])  #obs_mat[5]  #step 2
            #print(int(gravlen[ii,1]))
            sm = self.smooth_matrix(int(gravlen[ii,2]), 2) #order equals to two. #step 3
            #print(sm)
            print('segment g(mGal) max = %f/ min =%f'%(max(obs_mat[3]),min(obs_mat[3])))
            print('segment time(hr) max = %f/ min =%f'%(max(obs_mat[2][:,1]),min(obs_mat[2][:,1])))
            um2 = obs_mat[0]
            um3 = obs_mat[1]
            if (self.adj_method.lower() == 'bay1'):
                dyn = obs_mat[3]/meter_sf
                dyn2 = - um2[:,0] - um2[:,1]
            else:
                dyn = obs_mat[3] - um2[:,0] - um2[:,1]
                dyn2 = obs_mat[3] - um2[:,0] - um2[:,1]
            #print(ii)
            if (ii == 0):
                udd = ud
                uss = sm
                dbb = dyn
                urr = um3
                dbb2 = dyn2
            else:
                udd = slg.block_diag(udd,ud)               
                uss = slg.block_diag(uss,sm)                
                dbb = np.hstack([dbb, dyn])
                urr = np.vstack([urr, um3])
                dbb2 = np.hstack([dbb2, dyn2])
                
            dcs.append(dcnos)
            dce.append(dcnoe)

        mat_list =[]
        mat_list.append(udd)
        mat_list.append(uss)
        mat_list.append(dbb)
        mat_list.append(urr)
        mat_list.append(uag)
        mat_list.append(dag)
        mat_list.append(wag)
        mat_list.append(dbb2) #only be used when adj_method = bay1
        #print(uag.shape)
        dno_list =[]
        dno_list.append(pnt_id0) #slon, slat, selev = self._get_pnt_loc(slocid)
        dno_list.append(dcs)
        dno_list.append(dce)
        
        return mat_list, dno_list 
    
    @staticmethod
    def smooth_matrix(length, order):
        """Compute discrete derivative operators.
        Computes the discrete approximation L to the derivative operator
        of order d on a regular grid with n points, i.e. L is (n-d)-by-n.
           L is stored as a matrix.    
           eg. L_mat = Campaign.smooth_matrix(6,2)
        Modified from Matlab code by Per Christian Hansen, IMM, 02/05/98.
        
        """
        if not isinstance(order, int):
            raise ValueError('Order must be an integer!')
        if order < 0:
            raise ValueError('Order must be nonnegative!')
        # Zero'th derivative.
        if (order == 0):
            smooth_mat = np.eye(length)
            return smooth_mat 
        # Compute smooth_mat.         
        c = np.concatenate((np.array([-1, 1]),np.zeros(order-1))).astype(int)
        nd = length - order
        for i in range(1, order):
            c = np.append([0],c[0 : order]) -  np.append(c[0 : order], [0])

        smooth_mat = np.zeros([nd, length])
        
        for i in range(order + 1):
            tmp_mat = np.zeros([nd, length])
            for j in range(nd):
                tmp_mat[j,j+i] = c[i]
                
            smooth_mat = smooth_mat + tmp_mat
            
        return smooth_mat.astype(int)
    
    @staticmethod
    def gen_drift(hrs, data_num, day_num):
        """Generate drift matrix from the observed time recording
           eg.  d_mat = Campaign.gen_drift(*,*,*,*)
        
        """
        kt0 = sum(data_num)-len(data_num) #有效测段数
        kt2 = len(set(day_num)) #len(data_num) #测量天数 //测量记录数
        uarray = np.zeros([kt0, kt2])
        #print(data_num)
        #print(day_num)
        #print(len(data_num))
        #print(kt0, kt2)
        k = 0
        ii = 0
        for i in range(len(data_num)): #kt2) : #total days
            for j in range(data_num[i]-1): #data(i)-1 total segment in each loop
                uarray[k, ii] = hrs[k,1]
                if (j == (data_num[i]-2)):
                    kk = sum(data_num[:i+1]) -1
                    #print(kk,ii,j)
                    if (i<len(data_num)-1):
                        if (day_num[kk] != day_num[kk+1]):
                            #print(kk, day_num[kk+2])
                            ii +=1
                #print(kk, data_num[i], ii)
                k = k+1
        #print(uarray.shape)
        return uarray

    @staticmethod
    def check_hs(hrs, Dlength):
        
        is_good = False
        hr_max = max(hrs[:,0]+hrs[:,1])
        if (hr_max <= Dlength):
            is_good = True 
        else:
            print(max(hrs[:,0], min(hrs[:,0])))  
        return is_good, hr_max
    
    @timeit('runadj')
    def run_adj(self, filename = '', method = 1, maxiter = 1000):
        #import geoist.gravity.adjmethods as adj
        if (len(self.mat_list)<1):
            raise ValueError('Please run pre_adj() to generate Matrix!')

        if self.adj_method == 'cls':
            # initial adjustment class
            print(self.adj_method)

            xinit = 0.01
            xopt = adj.Clsadj.goadj(self.mat_list, self._gravlen, xinit) 
            print('The optimization has finished. AIC value is = %f'%xopt.fun)
            for ii in xopt.x :
                print(np.sqrt(np.exp(ii))*1000, 'uGal')
            xx, err, res = adj.Clsadj.result(xopt.x,self.mat_list, self._gravlen)
            dlen = len(self._gravlen)
            self.survey_dic['weight_SD_uGal'] = (np.sqrt(np.exp(xopt.x))*1000).tolist()
            self.survey_dic['drift_uGal_hr'] = np.squeeze(np.array(xx[0:dlen]*1000)).tolist()
            self.survey_dic['drift_err_uGal_hr'] = (err[0:dlen]*1000).tolist()
            self.survey_dic['staid'] = self.obs_list[0]
            self.survey_dic['gvalue_mGal'] = np.squeeze(np.array(xx[dlen:])).tolist()
            self.survey_dic['gerror_mGal'] = err[dlen:].tolist()
            self.survey_dic['obsres_mGal'] = np.squeeze(np.array(res)).tolist()
            
        elif self.adj_method == 'bay':
            print(self.adj_method)
            #print(self._gravlen)

            xinit = 0.01
            dinit = 1.0 #initial drift SD value
            xopt = adj.Bayadj.goadj(self.mat_list, self._gravlen, xinit, dinit, method, maxiter)

            print('The optimization has finished. ABIC value is = %f'%xopt.fun)
            #print(xopt.x)
            for ii in xopt.x :
                print(np.sqrt(np.exp(ii))*1000, 'uGal')
            xx, err, res = adj.Bayadj.result(xopt.x,self.mat_list, self._gravlen)
            dlen = len(self._gravlen)
            dobs = len(self.obs_list[0])
            self.survey_dic['weight_SD_uGal'] = (np.sqrt(np.exp(xopt.x))*1000).tolist()
            self.survey_dic['drift_uGal_hr'] = np.squeeze(np.array(xx[dobs:]*1000)).tolist()
            self.survey_dic['drift_err_uGal_hr'] = (err[dobs:]*1000).tolist()
            self.survey_dic['staid'] = self.obs_list[0]
            self.survey_dic['gvalue_mGal'] = np.squeeze(np.array(xx[:dobs])).tolist()
            self.survey_dic['gerror_mGal'] = err[:dobs].tolist()
            self.survey_dic['obsres_mGal'] = np.squeeze(np.array(res)).tolist()

        elif self.adj_method == 'bay1':
            print(self.adj_method)
            xinit = 0.01
            dinit = 1.0  #initial drift SD value
            sfinit = 1.0
            kstart = self.survey_list[0].meter_sf_index #不用标定格值的仪器数
            xopt = adj.Bayadj1.goadj1(self.mat_list, self._gravlen,
                                     xinit, dinit, sfinit, kstart)
            print('The optimization has finished. ABIC value is = %f'%xopt.fun)
            for ii in xopt.x :
                print(np.sqrt(np.exp(ii)))
            xx, err, res = adj.Bayadj1.result1(xopt.x,self.mat_list, self._gravlen)
            dlen = len(self._gravlen)
            dobs = len(self.obs_list[0])
            self.survey_dic['weight_SD_uGal'] = (np.sqrt(np.exp(xopt.x))).tolist()
            self.survey_dic['drift_uGal_hr'] = np.squeeze(np.array(xx[0:dlen]*1000)).tolist()
            self.survey_dic['drift_err_uGal_hr'] = (err[0:dlen]*1000).tolist()
            self.survey_dic['staid'] = self.obs_list[0]
            self.survey_dic['gvalue_mGal'] = np.squeeze(np.array(xx[:dobs])).tolist()
            self.survey_dic['gerror_mGal'] = err[:dobs].tolist()
            self.survey_dic['obsres_mGal'] = np.squeeze(np.array(res)).tolist()

        elif self.adj_method == 'bay2':
            print(self.adj_method) #TO-DO list
        else:
            print(self.adj_method) #TO-DO list

        if (len(filename)>1):
            try:
                f = open(filename, mode = 'w')
                f.write(json.dumps(self.survey_dic))
                f.close
            except IOError:
                print('No file : %s' %(filename))
            except ValueError:
                print('check raw data file')
            except IndexError:
                print('check raw data file: possibly last line?')
        else:
            print('please check survey_dic results')
            #return self.survey_dic # need remove @

    def save_mat_hdf5(self, filename):
        """save matrix to HDF5 file with Binary
           TO-DO by Bei
        """
        h = h5py.File(filename, 'w')
        h.create_dataset('udd', data = self.mat_list[0])
        h.create_dataset('uss', data = self.mat_list[1])
        h.create_dataset('ubb', data = self.mat_list[2])
        h.create_dataset('urr', data = self.mat_list[3])
        h.create_dataset('uag', data = self.mat_list[4])
        h.create_dataset('dag', data = self.mat_list[5])
        h.create_dataset('wag', data = self.mat_list[6])
        h.create_dataset('glen', data = self._gravlen)
        h.close()
    def save_mat_npz(self, filename):
        """save matrix to NPY/NPZ format file before transfrom to sparse with CSR format
           TO-DO by Bei
        """
        pass

    def save_sparse_mat_hdf5(self, filename):
        """save matrix to HDF5 file before transfrom to sparse with CSR format
           TO-DO by Bei
        """
        #b = csr_matrix(arr)
        #csr_matrix((data,indices,indptr),shape=(4,5)).toarray()
        pass
    def save_mat_json(self, filename): 
        """Save matrix to JSON file with Ascii
                mat_list:udd/uss/dbb/urr/
                obs_list:dcs/dce
        """
        save_dict = {}
        save_dict['udd'] = self.mat_list[0].tolist()
        save_dict['uss'] = self.mat_list[1].tolist()
        save_dict['ubb'] = self.mat_list[2].tolist()
        save_dict['urr'] = self.mat_list[3].tolist()
        save_dict['uag'] = self.mat_list[4].tolist()
        save_dict['dag'] = self.mat_list[5].tolist()
        save_dict['wag'] = self.mat_list[6].tolist()
        save_dict['sid'] = self.obs_list[0]
        save_dict['dcs'] = self.obs_list[1]
        save_dict['dce'] = self.obs_list[2]
        save_dict['glen'] = self._gravlen
        try:
            f = open(filename, mode = 'w')
            f.write(json.dumps(save_dict))
            f.close
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('check raw data file')
        except IndexError:
            print('check raw data file: possibly last line?')   
            
    def export_camp_json(self, filename):
        """Export instance of Campaign class to JSON file
        """        
        try:
            f = open(filename, mode = 'w')
            f.write(json.dumps(self, default=lambda o: getattr(o, '__dict__', str(o))))
            f.close
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('check raw data file')
        except IndexError:
            print('check raw data file: possibly last line?')          


def get_2d_list_slice(matrix, start_row, end_row, start_col, end_col):
    return [row[start_col:end_col] for row in matrix[start_row:end_row]]

def flatten(items): 
    """Yield items from any nested iterable; see REF."""
    from collections.abc import Iterable
    for x in items: 
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)): 
            yield from flatten(x) 
        else: 
            yield x

if __name__ == '__main__':
    
    m1 = Meter('LCR','G147')
    #m1.read_table('./data/table1.dat')
    m2 = Meter('LCR','G570')
    #m2.read_table('./data/table1.dat')
    n1 = Network('NorthChina',1)
    n1.read_pnts('./data/hball8.txt')
    print(n1)
    s1 = Survey('HBtest', '200901')
    s1.add_meter(m1)
    s1.add_meter(m2)
    s1.net = n1
    #s1.meter_sf_index = 1 # if select bay1, please check this parameter
    s1.read_survey_file('./data/simd1nf.G147') #sd1wn12.G147
    s1.read_survey_file('./data/simd2dn.G570')
    s1.corr_aux_effect()
    print(s1)         
    ag = AGstation('地球所','11000220','A', 116.31, 39.946, 51.9)
    ag.ref_gra = 0.0
    ag.ref_gra_err = 2.1E-3 
    print(ag)
    gravwork = Campaign('IGP201702', 1)
    gravwork.add_ag_sta(ag)                  #添加绝对点信息
    gravwork.add_surveys(s1)        #添加测量到平差任务
    print(gravwork)
    gravwork.adj_method = 1 #1:cls ; 2:Baj; 3:Baj1
    if gravwork.pre_adj():
        print(len(gravwork.mat_list[0]))
        gravwork.run_adj('./data/grav_baj2.txt', 3, 1000) #1:simplex 2:BFGS
    #aa = json.load(open('./data/grav_baj.txt'))
    #gravwork.export_camp_json('./data/gravwork.json') #保存对象示例到磁盘
    #gravwork.save_mat_json('./data/gravmat.json') #保存平差矩阵到磁盘
    #gravwork.save_mat_hdf5('./data/gravmat1.h5')    