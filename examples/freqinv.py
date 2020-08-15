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
from matplotlib import image 
# local imports
from geoist import gridder
from geoist.inversion import geometry
from geoist.pfm import prism,pftrans
from geoist.vis import giplt
from geoist.pfm import giutils
from geoist.pfm.giconstants import G, SI2EOTVOS, CM, T2NT, SI2MGAL

## 1. generate a model
def gen_grav(shape = (100, 100)):

    model = [geometry.Prism(-1000, 1000, -1000, 1000, 100, 2000, {'density': 1000}),
            geometry.Prism(-4000, -3000, -4000, -3000, 100, 2000, {'density': -900}),
            geometry.Prism(2000, 4000, 3000, 4000, 100, 2000, {'density': 1300})]
    xp, yp, zp = gridder.regular((-5000, 5000, -5000, 5000), shape, z=0)
    field0 = prism.gz(xp, yp, zp, model) 
    #prism.potential(xp, yp, zp, model)
    #field1 = giutils.contaminate(field0, 0.05, percent = True)
    fieldfreq = gzfreq(xp, yp, zp, shape, model)
    # fields = [prism.gx(xp, yp, zp, model),
    #           prism.gy(xp, yp, zp, model),
    #           prism.gz(xp, yp, zp, model),
    #           prism.gxx(xp, yp, zp, model),
    #           prism.gxy(xp, yp, zp, model),
    #           prism.gxz(xp, yp, zp, model),
    #           prism.gyy(xp, yp, zp, model),
    #           prism.gyz(xp, yp, zp, model),
    #           prism.gzz(xp, yp, zp, model)]
    
    return field0, fieldfreq, xp, yp, zp

## 2. look the spetrum with AMP and phase
def gzfreq(xp, yp, zp, shape, prisms, dens=None):
    """
    Calculates the :math:`g_z` gravity acceleration component in frequency domain.
    Returns:

    * res : array
        The field calculated on xp, yp, zp

    """
    if xp.shape != yp.shape or xp.shape != zp.shape:
        raise ValueError("Input arrays xp, yp, and zp must have same length!")
    size = len(xp)
    res = np.zeros(size, dtype=np.float)
    for prism in prisms:
        if prism is None or ('density' not in prism.props and dens is None):
            continue
        if dens is None:
            density = prism.props['density']
        else:
            density = dens
        x1, x2 = prism.x1, prism.x2
        y1, y2 = prism.y1, prism.y2
        z1, z2 = prism.z1, prism.z2
        r1 = _gzfreq(xp, yp, zp, x1, x2, y1, y2, z1, z2, shape)
        res = res + r1 *2*np.pi*G *density* SI2MGAL
    #res *= 
    return res

def _gzfreq(x, y, data, x1, x2, y1, y2, z1, z2, shape):
    
    x0 = (x1+x2)/2.0
    y0 = (y1+y2)/2.0
    a = np.abs(x1-x2)/2.0
    b = np.abs(y1-y2)/2.0
    nx, ny = shape
    dx = (x.max() - x.min())/(nx - 1)
    dy = (y.max() - y.min())/(ny - 1)
    print(x0,y0,a,b,z1,z2)
    
    nx, ny = shape
    # Pad the array with the edge values to avoid instability
    padded, padx, pady = _pad_data(data, shape)
    kx, ky = _fftfreqs(x, y, shape, padded.shape)
    kz = np.sqrt(kx**2 + ky**2)
    kxm = np.ma.array(kx, mask= kx==0)
    kym = np.ma.array(ky, mask= ky==0)
    kzm = np.ma.array(kz, mask= kz==0)    
    # print(kx.shape,ky.shape,kz.shape)
    # plt.imshow(kx)
    # plt.show()
    # plt.imshow(ky)
    # plt.show()
    complex1 = 0+1j
    ker1 = (np.exp(-kz*z2)-np.exp(-kz*z1))/kzm
    ker2 = np.sin(ky*b)/kym
    ker3 = np.sin(kx*a)/kxm
    keruv = -4*ker1*ker2*ker3
    
    keruv[kxm.mask] = -4*a*np.sin(ky[kxm.mask]*b)*ker1[kxm.mask]/ky[kxm.mask]
    keruv[kym.mask] = -4*b*np.sin(kx[kym.mask]*a)*ker1[kym.mask]/kx[kym.mask]
    keruv[kzm.mask] = 4*a*b*(z2-z1)
    nxe, nye = padded.shape
    
    M_left=(nxe-nx)/2+1
    M_right=M_left+nx-1
    N_down=(nye-ny)/2+1
    N_up=N_down+ny-1    
    
    XXmin=x.min()-dx*(M_left-1)
    XXmax=x.max()+dx*(nxe-M_right)
    YYmin=y.min()-dy*(N_down-1)
    YYmax=y.max()+dy*(nye-N_up)    
    
    keruv = keruv*np.exp(-ky*y0*complex1)*np.exp(-kx*x0*complex1) #*np.exp(kz*data[0])
    #keruv = kz * keruv 
    keruv = keruv*np.exp(complex1*((x.max()+x.min())*kx/2+(y.max()+y.min())*ky/2))*np.exp(complex1*((XXmin-XXmax)*kx/2+(YYmin-YYmax)*ky/2))/dx/dy

    kshif = np.fft.fftshift(keruv)
    plt.imshow(np.log(np.abs(kshif)))
    plt.show()
    plt.imshow(np.angle(kshif))
    plt.show()
    
    res = np.real(np.fft.ifft2(keruv))
    print(res.max(),res.min())
    res0 = res[padx: padx + nx, pady: pady + ny].ravel()
    return res0
 

def _pad_data(data, shape):
    n = _nextpow2(np.max(shape))
    nx, ny = shape
    padx = (n - nx)//2
    pady = (n - ny)//2
    padded = np.pad(data.reshape(shape), ((padx, padx), (pady, pady)),
                       mode='edge')
    return padded, padx, pady

def _nextpow2(i):
    buf = np.ceil(np.log(i)/np.log(2))
    return int(2**buf)

def _fftfreqs(x, y, shape, padshape):
    """
    Get two 2D-arrays with the wave numbers in the x and y directions.
    """
    nx, ny = shape
    dx = (x.max() - x.min())/(nx - 1)
    fx = 2*np.pi*np.fft.fftfreq(padshape[0], dx)
    dy = (y.max() - y.min())/(ny - 1)
    fy = 2*np.pi*np.fft.fftfreq(padshape[1], dy)
    return np.meshgrid(fy, fx)[::-1]

class prismfreq():
    pass
## 3. wave2num

if __name__ == '__main__':
    
    print('hello freqinv!')
    shape = (156, 156)
    gu, gulist, xp, yp, zp = gen_grav(shape)
    titles = ['potential', 'gx', 'gy', 'gz',
              'gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
    plt.figure(figsize=(8, 8))
    plt.axis('scaled')
    plt.title(titles[0])
    levels = giplt.contourf(yp * 0.001, xp * 0.001, gu, shape, 15)
    cb = plt.colorbar()
    giplt.contour(yp * 0.001, xp * 0.001, gu, shape,
                levels, clabel=False, linewidth=0.1)
    plt.show()
    
    plt.figure(figsize=(8, 8))
    plt.axis('scaled')
    plt.title('gz by freq')
    levels = giplt.contourf(yp * 0.001, xp * 0.001, gulist, shape, 15)
    cb = plt.colorbar()
    giplt.contour(yp * 0.001, xp * 0.001, gulist, shape,
                levels, clabel=False, linewidth=0.1)
    # plt.subplots_adjust(left=0.03, right=0.95, bottom=0.05, top=0.92, hspace=0.3)
    # plt.suptitle("Potential fields produced by a 3 prism model")
    # for i, field in enumerate(gulist):
    #     plt.subplot(3, 3, i + 1)
    #     plt.axis('scaled')
    #     plt.title(titles[i+1])
    #     levels = giplt.contourf(yp * 0.001, xp * 0.001, field, shape, 15)
    #     cb = plt.colorbar()
    #     giplt.contour(yp * 0.001, xp * 0.001, field, shape,
    #                 levels, clabel=False, linewidth=0.1)
    plt.show()    
    df1 = pd.DataFrame(columns=['x','y','z','g1','g2'])
    df1['x'] = yp
    df1['y'] = xp
    df1['z'] = zp
    df1['g1'] = gu
    df1['g2'] = gulist
    df1.to_csv('D:\\MyResearch\\abic_inversion\\freqinv\\freq.txt')
    
    # img = image.imread('D:\\testhli\\realdata-forward.png') 
    # gray_img = np.dot(img[:,:,:3], [.21, .72, .07]) 
    # gray_img.shape 
    # plt.imshow(gray_img, cmap = plt.get_cmap('gray')) 
    # plt.show()
    # # fft2 是二维数组的傅里叶变换
    # # 将空域转换为频域
    # fft = np.fft.fft2(gray_img) 
    # amp_spectrum = np.abs(fft) 
    # plt.imshow(np.log(amp_spectrum)) 
    # plt.show()    
    # angle_spectrum = np.angle(fft) 
    # plt.imshow(angle_spectrum) 
    # plt.show()   
    # fft_shift = np.fft.fftshift(fft) 
    # plt.imshow(np.log(np.abs(fft_shift))) 
    # plt.show() 
    # m, n = fft_shift.shape
    # b = np.zeros((int(m / 2), n)) 
    # c = np.zeros((2 * m, int(n / 2))) 
    # fft_shift = np.concatenate((b, fft_shift, b), axis = 0) 
    # fft_shift = np.concatenate((c, fft_shift, c), axis = 1) 
    
    # # 然后再转换回去
    # ifft = np.fft.ifft2(np.fft.ifftshift(fft_shift)) 
    # ifft.shape 
    # # (633L, 1321L) 
    # ifft = np.real(ifft) #np.imag(ifft)
    # plt.imshow(ifft, cmap = plt.get_cmap('gray')) 
    # plt.show() 
