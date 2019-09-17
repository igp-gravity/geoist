# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 22:25:52 2019


@Author:  Shi CHEN   @ IGP-CEA
         Bei ZHANGE @ IGP-CEA

Co-Author: Jiancang ZHUANG @ ISM-Tokyo

####################################################
            MIT license
        Copyright @ pyBACGS 2019 
         All right reserved 
####################################################

Module which contains the main classes of the pyBACGS
            
CLASS list:
    Genmodel 

"""

class Genmodel(object):
    """simulated Test Model Dataset
    
    Properties:
                    
    functions:
    - extract_subset
                
    """       
    method = []
    campaign = [] 
    def __init__(self, fname, gscale):
        """
        """
        self.mtype=fname
        self.mscale=gscale
        
    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('%s仪器的格值为%d'%(self.mtype,self.mscale))       
        return "Print by BACGS"


def gen_obs_noise():
    pass
