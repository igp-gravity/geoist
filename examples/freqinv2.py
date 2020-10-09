# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 18:29:09 2020
D:\MyResearch\abic_inversion\freqinv

Gravity inversion in frequency domain
@author: chens
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# local imports
from geoist.inversion import freqinv
from geoist.pfm import prism
from geoist.vis import giplt
#from pyevtk.hl import gridToVTK
import datetime
import pathlib

if __name__ == '__main__':
    # define model volume and grid shape

    work_dir = pathlib.Path('/home/zhangb/work/people/zhangbei/freq')
    nzyx = [16,16,16] # 源网格剖分，分别是z,y,x方向网格数
    source_volume = [-3000, 3000, -3000, 3000, 100, 2000] # 源的分布范围，分别是 xmin,xmax,ymin,ymax,zmin,zmax
    obs_area = [-3000,3000,-3000,3000] # 观测场的范围，分别是xmin,xmax,ymin,ymax
    nobsyx = [64,64] # 观测场的格点数，分别是y,x方向的格点个数

    model_density = np.zeros(tuple(nzyx)) # 正演用的密度模型
    model_density[10:12,5:6,7:8] = 2000

    refer_densities = [] # 反演用的参考密度模型(reference)列表，列表中每一个成员是一个参考密度模型。
    weights = {'refers':[1]} # 权重，这里只有reference的权重，没有光滑之类的。reference 权重的个数需要和reference的个数相等。
    refer_density = np.zeros(tuple(nzyx)) # 设置一个全0的reference,相当于最小模型约束
    refer_density = refer_density.ravel()
    refer_densities.append(refer_density) # 将全0的reference加入到reference列表里。

    small_model = freqinv.FreqInvModel(nzyx=nzyx,
                                       source_volume=source_volume,
                                       nobsyx=nobsyx,
                                       obs_area=obs_area,
                                       model_density=model_density,
                                       weights=weights,
                                       refer_densities=refer_densities)
    small_model.gen_mesh(height=1.)  # 生成网格
    small_model.gen_kernel() # 生成核矩阵
    small_model.forward(update=True) # 正演
                                     # 如果不给密度源，则使用small_model.model_density作为密度源.
                                     # 如果update=True,更新small_model.freq和small_model.obs_field为正演结果。
                                     # 如果update=False,则把计算结果返回，small_model.freq 和small_model.obs_field 不发生变化。
    true_freq = small_model.freq     # 备份真实的频率域重力
    small_model.obs_freq = small_model.freq # small_model.obs_freq,观测到的频率域重力。
    # add noise
    real_noise = 0.1*np.max(np.abs(small_model.obs_freq.real))
    imag_noise = 0.1*np.max(np.abs(small_model.obs_freq.imag))
    small_model.obs_freq += real_noise + imag_noise*1j
    print(np.max(small_model._model_density))
    print(len(small_model._model_density))
    print(np.conj(small_model.FX).shape,small_model.FY.shape)
    st = datetime.datetime.now()
    small_model.do_linear_solve_quiet() # 反演,结果保存在small_model.solution里。
    ed = datetime.datetime.now()
    print("inversion use time: ",ed-st)

# check solution
    freq,recover = small_model.forward(small_model.solution.real) # 使用反演的结果进行正演。
    print('field recover misfit: ',np.linalg.norm(recover - small_model.obs_field)/np.linalg.norm(small_model.obs_field)) #反演结果正演重力场和真实重力场相对误差。
    print('freq recover misfit: ',np.linalg.norm(freq - true_freq)/np.linalg.norm(freq)) # 反演结果正演频率域重力场和真实值相对误差。
    print('freq(abs) recover misfit: ',np.linalg.norm(np.abs(freq)
                                                      - np.abs(true_freq))/np.linalg.norm(np.abs(freq))) # 反演结果正演得到的频率域重力的模和真实值相对误差。
    angle_diff = (np.angle(freq) - np.angle(true_freq)) % (np.pi*2)
    print('freq(angle) recover misfit: ',np.linalg.norm(angle_diff)/np.linalg.norm(np.angle(freq))) # 频率域幅角和真实值的相对误差(现在的算法有问题）。
    rhs = small_model.C*small_model.W*freq.T.ravel()
    rhs = small_model.dW.T*rhs
    rhs = freqinv.kron_matvec([np.conj(small_model.FY.T),
                               np.conj(small_model.FX.T)],
                              rhs.reshape(small_model.nz,-1)).ravel()
    print('rhs recover misfit: ',np.linalg.norm(rhs - small_model.rhs)/np.linalg.norm(small_model.rhs)) # 方程右端项的检查。

# save density model to vtk, comment out this section if gridToVTK not installed.

#    nz,ny,nx = nzyx
#    arr = small_model.solution.real.copy()
#    xs = np.zeros((nz+1,ny+1,nx+1)) + np.array(small_model.mesh.get_xs()).reshape(1,1,-1)
#    ys = np.zeros((nz+1,ny+1,nx+1)) + np.array(small_model.mesh.get_ys()).reshape(1,-1,1)
#    zs = np.zeros((nz+1,ny+1,nx+1)) + np.array(small_model.mesh.get_zs()).reshape(-1,1,1)
#    gridToVTK(str(work_dir/"results"),xs,ys,zs,cellData={"recover":arr.reshape(nz,ny,nx),
#                                             "orig":model_density})
#

    field0 = prism.gz(small_model.xp, small_model.yp, small_model.zp, small_model.mesh)
    shape = (nobsyx[0],nobsyx[1])

# plot field
    fig = plt.figure(figsize=(16, 8))
    # original field
    axes = fig.subplots(2,2)
#    plt.axis('scaled')
    ca = axes[0][0]
    plt.sca(ca)
    levels = giplt.contourf(small_model.yp * 0.001, small_model.xp * 0.001, field0, shape, 15)
    cb = plt.colorbar()
    giplt.contour(small_model.yp * 0.001, small_model.xp * 0.001, field0, shape,
                levels, clabel=False, linewidth=0.1)
    ca.set_title('gravity anomlay')

    # field from frequency
    ca = axes[0][1]
    plt.sca(ca)
#    levels = giplt.contourf(small_model.yp * 0.001, small_model.xp * 0.001,
#    small_model.obs_field, shape, 15)
    levels = giplt.contourf(small_model.yp * 0.001, small_model.xp * 0.001, small_model.obs_field, shape, 15)
    cb = plt.colorbar()
    giplt.contour(small_model.yp * 0.001, small_model.xp * 0.001, small_model.obs_field, shape,
                levels, clabel=False, linewidth=0.1)
    ca.set_title('gz by freq')

    # field recovered from solution
    ca = axes[1][0]
    plt.sca(ca)
#    levels = giplt.contourf(small_model.yp * 0.001, small_model.xp * 0.001, recover, shape, 15)
    levels = giplt.contourf(small_model.yp * 0.001, small_model.xp * 0.001, recover, shape, 15)
    cb = plt.colorbar()
    giplt.contour(small_model.yp * 0.001, small_model.xp * 0.001, recover, shape,
                levels, clabel=False, linewidth=0.1)
    ca.set_title('recovered')
    ca = axes[1][1]
    ca.axis('off')
    plt.savefig(work_dir/pathlib.Path('field.jpg'))

# plot frequency
    fig = plt.figure(figsize=(16, 16))
    axes = fig.subplots(2,2)

    # spectrum observation
    freq = true_freq
    np.savetxt(work_dir/'true_freq.txt',freq)
    # mag
    ca = axes[0][0]
    ca.imshow(np.abs(np.fft.fftshift(freq)))
    ca.set_title('observation freq mag')
    # phase
    ca = axes[0][1]
    ca.imshow(np.angle(freq))
    ca.set_title('observation freq phase')
    print(freq.shape)

    # spectrum calculated
    freq,recover = small_model.forward(small_model.solution.real)
    np.savetxt(work_dir/'recover_freq.txt',freq)
    print(freq.shape)

    # mag
    ca = axes[1][0]
    ca.imshow(np.abs(np.fft.fftshift(freq)))
    ca.set_title('calculated freq mag')
    # phase
    ca = axes[1][1]
    ca.imshow(np.angle(freq))
    ca.set_title('calculated freq phase')
    plt.savefig(work_dir/pathlib.Path('frequency.jpg'))
# save data
#    df1 = pd.DataFrame(columns=['x','y','z','g_orig','g_recover'])
#    df1['x'] = small_model.yp
#    df1['y'] = small_model.xp
#    df1['z'] = small_model.zp
#    df1['g_orig'] = field0
#    df1['g_recover'] = small_model.obs_field
#    df1.to_csv('freq.txt')
#
