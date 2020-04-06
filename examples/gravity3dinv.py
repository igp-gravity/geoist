# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 21:02:14 2020

@author: chens
"""
 
from geoist.pfm.inv3d import ModelBuilder, RegularMesh, GravField
from geoist.pfm.inv3d import Dimension2, Dimension, Volume, HexSource

mb = ModelBuilder()    
t1mesh = RegularMesh('grid1')
mb.set_field(GravField('test1.grav','mGal', t1mesh))

dims = Dimension(2,4,8)
vols = Volume(1,2,3,4,5,6)
mb.set_source(HexSource('hex1.src','km'),dims, vols)    
    
mb.set_problem('prism', '3D')

m1 = mb.builder('model_test1')
print(m1)
weights = {'dx':23.54,'dy':860.,'dz':60.28,'depth':1.,'obs':1.69,'refer':0.67,'bound':1.0}
m1.set_weight(weights)
m1.source.gen_single_body(1.0)
m1.set_depth_constraint()
m1.run()     # generate model 

m1.gen_kernel()
m1.forward()

m1.model.min_density = 0.
m1.model.max_density = 1.
orig_obs = m1.model.obs_data.copy()
m1.model.walsh_transform()
m1.model.abic_val
list(m1.model._weights.keys())
m1.model.do_linear_solve()
m1.model.calc_abic()

#export model information to JSON file
m1.save_to_json('./data/m1json.txt')

mb2 = ModelBuilder()
mb2.set_field(GravField('test2.grav','mGal', RegularMesh('grid2')))
smooth_components = ['dx','dy','dxy','dt']
weights = {'obs':1,'dx':1,'dy':1,'dxy':1,'dt':1}

dims = Dimension2(20,20)
vols = Volume(-50000, 50000, -50000, 50000, 1000, 11000)
mb2.set_source(HexSource('hex2.src','km'),dims, vols)    
    
mb2.set_problem('Tesseroid', '4D')
m2 = mb2.builder('model_test2')    
print(m2)
m2.set_weight(weights)
m2.set_smooth(smooth_components)
m2.run()
m2.save_to_json('./data/m2json.txt')
#m2 = mb.builder_ts()
#m1.field.load_from_json()