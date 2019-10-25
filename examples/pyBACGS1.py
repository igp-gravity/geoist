# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 20:06:55 2019

@author: chens
"""
import json
import numpy as np
#local packages
import geoist.gravity.graobj as gg
import geoist.gravity.adjmethods as adj

m1 = gg.Meter('CG-5','C098')
m1.msf = 1.00009
m2 = gg.Meter('CG-5','C099')
m2.msf = 1.000637
m3 = gg.Meter('CG-5','C097')
m3.msf = 1.000163
#m3.read_table('./data/table.dat')
m4 = gg.Meter('CG-5','C0981')
m4.msf = 1.00009
n1 = gg.Network('NorthChina',1)
n1.read_pnts('./data/SDQSX8.DZJ')
print(n1)
s1 = gg.Survey('DQSW+SX', '201508')
s1.add_meter(m1)
s1.add_meter(m2)
s1.add_meter(m3)
s1.add_meter(m4)
s1.net = n1
s1.read_survey_file('./data/QSCW201508p.098')
s1.read_survey_file('./data/QSCW201508p.099')
s1.read_survey_file('./data/SXCW971508p.ori')
s1.read_survey_file('./data/SXCW981508p.ori', 'C0981')


s1.corr_aux_effect()
print(s1)
#查找一个测量工程中某个测点号对应的坐标
slon, slat, selev = s1._get_pnt_loc('11003902')
#查找一个测量工程中某台重力仪的格值
sf_val = s1._get_meter_sf('C097')

ag1 = gg.AGstation('白山洞绝对','11014121','A', 116.169, 40.018, 212.5)
ag1.ref_gra = 1110.54453
ag1.ref_gra_err = 5.0E-3 
#ag2 = gg.AGstation('张家口绝对','13072402','A', 114.902, 40.830, 859.  )
#ag2.ref_gra = 961.0019
#ag2.ref_gra_err = 5.0E-3 
#ag3 = gg.AGstation('代县地震台','35017301','A', 112.933, 39.050, 820.  )
#ag3.ref_gra = 757.52477
#ag3.ref_gra_err = 5.0E-3 
#ag4 = gg.AGstation('阳原地震台','10114700','A', 114.1493, 40.1256, 943.70 )
#ag4.ref_gra = 868.0904
#ag4.ref_gra_err = 5.0E-3 
#print(ag1,ag2,ag3,ag4)
gravwork = gg.Campaign('IGP201604', 1)
gravwork.add_ag_sta(ag1)                  #添加绝对点信息 可以添加多次
#gravwork.add_ag_sta(ag2)                  #添加绝对点信息 可以添加多次
#gravwork.add_ag_sta(ag3)                  #添加绝对点信息 可以添加多次
#gravwork.add_ag_sta(ag4)                  #添加绝对点信息 可以添加多次
gravwork.add_surveys(s1)        #添加测量到平差任务
print(gravwork)
#开始平差pre_adj是完成从观测文件重，生成平差矩阵的
gravwork.adj_method = 2 #1:cls ; 2:Baj; 3:Baj1
#写法1：平差结果导出到txt，json格式
#if gravwork.pre_adj():
#    #print(len(gravwork.mat_list[0]))
#    gravwork.run_adj('./data/grav_baj.txt')
#
#grav_dict = json.load(open('./data/grav_baj.txt')) #加载结果

#写法2：run_adj后，调用survey_dic
gravwork.pre_adj()
gravwork.run_adj('./data/grav_cls.txt',3) #也可以这样运行，没有文件名时直接返回结果
#grav_dict1 = gravwork.survey_dic

#写法3：如果不用run_adj函数可以这样写，直接用adj里面的平差算法
#xinit = 0.01
#xopt = adj.Clsadj.goadj(gravwork.mat_list, gravwork._gravlen, xinit)
#print('The optimization has finished. AIC value is = %f'%xopt.fun)
#for ii in xopt.x :
#    print(np.sqrt(ii)*1000, 'uGal')
#xx, err, res = adj.Clsadj.result(xopt.x,gravwork.mat_list, gravwork._gravlen)
