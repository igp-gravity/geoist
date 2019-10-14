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

m1 = gg.Meter('CG-5','C097')
m2 = gg.Meter('CG-5','C098')
m3 = gg.Meter('LCR','G570')
n1 = gg.Network('NorthChina',1)
n1.read_pnts('./data/QSCW912018.DZJ')
print(n1)
s1 = gg.Survey('HBDQSCW', '201604')
s1.add_meter(m1)
s1.add_meter(m2)
s1.add_meter(m3)
s1.net = n1
s1.read_survey_file('./data/QSCW201603.097')
s1.read_survey_file('./data/QSCW201603.098')
s1.read_survey_file('./data/QSCW201604.570')
s1.corr_aux_effect()
print(s1)
#查找一个测量工程中某个测点号对应的坐标
slon, slat, selev = s1._get_pnt_loc('11003902')
#查找一个测量工程中某台重力仪的格值
sf_val = s1._get_meter_sf('C097')

ag = gg.AGstation('香山外','11003902','A', 116.150, 40.025, 151.9   )
ag.ref_gra = 0.0
ag.ref_gra_err = 2.1E-3 
print(ag)
gravwork = gg.Campaign('IGP201604', 1)
gravwork.add_ag_sta(ag)                  #添加绝对点信息 可以添加多次
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
gravwork.run_adj() #也可以这样运行，没有文件名时直接返回结果
grav_dict1 = gravwork.survey_dic

#写法3：如果不用run_adj函数可以这样写，直接用adj里面的平差算法
#xinit = 0.01
#xopt = adj.Clsadj.goadj(gravwork.mat_list, gravwork._gravlen, xinit)
#print('The optimization has finished. AIC value is = %f'%xopt.fun)
#for ii in xopt.x :
#    print(np.sqrt(ii)*1000, 'uGal')
#xx, err, res = adj.Clsadj.result(xopt.x,gravwork.mat_list, gravwork._gravlen)
