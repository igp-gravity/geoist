# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:08:26 2019

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
    Adjmethod 

"""
import numpy as np
import h5py, json
import scipy.optimize as opt
import scipy.linalg as slg

class Adjustment(object):
    """Base Class of Gravity Adjustment Methods
    
    Properties:
                    
    functions:
    - extract_subset
                
    """       
    method = []
    obs_mat = []
    result = []
    error = []
    residual = []

    def __init__(self, name):
        """
        """
        self._name = name

    def __str__(self):
        """Override the built-in method 'print' when applied to such object
        
        """
        print('Adjustment name is %s'%(self._name))       
        return "====Adjustment Print by BACGS===="
    
    def load_json_mat(self, filename):
        """load matrix from JSON file
        """
        data_dict = {}
        try:
            f = open(filename, mode = 'r')
            data_dict = json.load(f)
            f.close
        except IOError:
            print('No file : %s' %(filename))            
        except ValueError:
            print('check raw data file')
        except IndexError:
            print('check raw data file: possibly last line?')           
        return data_dict
    
    def load_hdf5_mat(self, filename):
        """load matrix from HDF5 file
        """
        data_dict = {}
        h = h5py.File(filename,'r') 
        for key in h.keys():
            #print(key)
            print(h[key].name)
            print(h[key].shape)
            data_dict[key] = h[key].value
            #print(h[key].value)
        return data_dict
    @staticmethod
    def forward(W, A, b):
        """Compute X and residual using the optimized weigths
        Input Matrix: 
           W : diagnoal weights matrix 
           A : core matrix
           b : observed vector
        Output:
           X   : the adjustment result
           err : the estimated error of result
           res : the residual noise  
        """
        X = np.linalg.inv(A.T*W*A)*A.T*W*b
        err = np.sqrt(np.diag(np.linalg.inv(A.T*W*A)))
        res = b- A*X
        return X, err, res


    @staticmethod
    def optimization(likelihood, x0, args, method = 1, maxiter = 1000, xmin = [], xmax = [] ):
        """There are 4 method available. 
        1: Nelder-Mead simplex ; 
        2: Netwon; 
        3: Dual Annealing, can set bounds
        4: Basinhopping, can set bounds
        """
        Nfeval = 1
        #print(Nfeval, id(Nfeval))
        def callback(xk):
            nonlocal Nfeval
            xx = np.sqrt(np.exp(xk[:]))*1000
            #print(id(Nfeval))
            if Nfeval%10 == 0:
                if (len(xk) == 2):
                    print('{0:4d} {1: 3.6f} {2: 3.6f}'.format(Nfeval, *xx))
                elif (len(xk) == 4):
                    print('{0:4d} {1: 3.6f} {2: 3.6f} {3: 3.6f} {4: 3.6f}'.format(Nfeval, *xx))
                else:
                    print('{0:4d} {1: 3.6f} {2: 3.6f} {3: 3.6f} {4: 3.6f} ...'.format(Nfeval, *xx))

            Nfeval += 1

        if (method == 1):
            print('Nelder-Mead simplex method used for optimization')
            if (len(x0) == 2):
                print('{0:4s}   {1:9s}   {2:9s} '.format('Iter', ' X1', ' X2'))
            elif (len(x0) == 4):
                print('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s} '.format('Iter', ' X1', ' X2', ' X3', ' X4'))
            else:
                print('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s} ...'.format('Iter', ' X1', ' X2', ' X3', ' X4'))
            xopt = opt.minimize(likelihood, x0, args, method='nelder-mead',
                                options={'maxiter':maxiter,'disp': True}, callback = callback)
        elif (method == 2):
            #test = lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2
            #xopt = opt.fmin(func = likelihood, x0)
            print('BFGS method used for optimization')
            xopt = opt.minimize(likelihood, x0, args, method='BFGS',
                                options={'maxiter':maxiter, 'disp': True})
        elif (method == 3):
            print('L-BFGS-B method used for optimization')
            xopt = opt.minimize(likelihood, x0, args, method='L-BFGS-B',
                                options={'disp': None})
        else:
            #x0 = [10., 10.] # the starting point
            #xmin = [1., 1.] # the bounds
            #xmax = [11., 11.]
            # rewrite the bounds in the way required by L-BFGS-B
            bounds = [(low, high) for low, high in zip(xmin, xmax)]
            # use method L-BFGS-B because the problem is smooth and bounded
            minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)
            xopt = opt.basinhopping(likelihood, x0, minimizer_kwargs=minimizer_kwargs)        
        
        return xopt
    
class Clsadj(Adjustment):
    """Classical Adjustment using Linear drift model
    A derived class of Adjustment
    the staticmethod and classmethod have been designed.
    the instance of object has not be needed.
    Properties:
       NONE             
    functions:
    - likelihood: the user likelihood with the defined problem
    - goadj : to optimize the weights or hyperparameters
    - result : to forward the gravity values with optimized weights                
    """   
    @staticmethod
    def likelihood(x, *args):
        """给定超参数后，计算AIC数值, Max(ln(L))
        W       S         W Y       
        |Wrg   ||Udd  Urr| =   |Dbb|
        |   Wag||     Uag|     |dag|
        W~1/sigma^2  
        Input params(S,Y) must be np.matrix
        """
        S, Y, wag, gravlen = args
        lens = len(gravlen)
        
        S = np.matrix(S)
        Y = np.matrix(Y).T
        k = 0 
        w = np.ones(gravlen[0,0])/np.exp(x[0])
        W = np.diag(w)
        for k in range(1, lens):
            w = np.ones(gravlen[k,0])/np.exp(x[k])
            W = slg.block_diag(W, np.diag(w))
       
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))       
        W = np.matrix(W)
        ff = np.diag(np.linalg.inv(W))
        #print(sum(np.log(ff)))
        #print(Y.shape)
        #print(S.shape)
        #print(W.shape)
        f1 = sum(np.log(ff))
        #print(min(ff))
        X = np.linalg.inv(S.T*W*S)*S.T*W*Y
        f2 = np.sum((S*X-Y).T*W*(S*X-Y))  #min AT P A
        return f1 + f2
    @classmethod
    def result(cls, xopt, mat_list, glen):
        """forward result with the optimated weights
        """
        udd = mat_list[0]
        ubb = mat_list[2]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]

        m, n = udd.shape
        aa = np.zeros([m, len(glen)], dtype = float)
        k = 0
        for ii in range(len(glen)):
            kk = glen[ii, 0]
            for jj in range(kk):
                aa[k, ii] = sum(udd[k,:])
                k = k + 1
        
        aa = np.hstack([aa ,mat_list[3]])
        aa = np.vstack([aa ,uag[:, n-len(glen):]])
        bb = np.hstack([ubb, dag])        
        
        A = np.matrix(aa)
        b = np.matrix(bb).T

        k = 0
        w = np.ones(glen[0,0])/np.exp(xopt[0])
        W = np.diag(w)
        for k in range(1, len(glen)):
            w = np.ones(glen[k,0])/np.exp(xopt[k])
            W = slg.block_diag(W, np.diag(w))
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))       
        W = np.matrix(W)
        
        return cls.forward(W, A, b)
    @classmethod
    def goadj(cls, mat_list, glen, xinit):
        """Go adjustment using the matrix list
        the inversed weights will be estimated.
        """
        udd = mat_list[0]
        ubb = mat_list[2]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]
        m, n = udd.shape
        aa = np.zeros([m, len(glen)], dtype = float)
        k = 0
        for ii in range(len(glen)):
            kk = glen[ii, 0]
            for jj in range(kk):
                aa[k, ii] = sum(udd[k,:])
                k = k + 1
        
        aa = np.hstack([aa ,mat_list[3]])
        aa = np.vstack([aa ,uag[:, n-len(glen):]])
        bb = np.hstack([ubb, dag])        
        
        x0 = np.log(np.ones(len(glen))*xinit**2)
        #print(x0.shape)
        ua = np.matrix(aa)
        ub = np.matrix(bb)
        #uw = np.matrix(uw)
        args = (ua, ub, wag, glen)
        return cls.optimization(cls.likelihood, x0, args)

class Bayadj(Adjustment):
    """Bayesian Adjustment using non-Linear drift model
    A derived class of Adjustment
    the staticmethod and classmethod have been designed.
    the instance of object has not be needed.
    Properties:
    NONE             
    functions:
    - likelihood: the user likelihood with the defined problem
    - goadj : to optimize the weights or hyperparameters
    - result : to forward the gravity values with optimized weights                
    """     
    bsm = 2 #the order of smooth of meter drift

    @staticmethod
    def likelihood(x, *args): #x,Uall,Uss,Y,gravlen,wag
        """给定超参数后，计算ABIC数值
        W       S         W Y       
        |Wrg       ||Urr  Udd| =   |Dbb|
        |   Wag    ||Uag     |     |dag|
        |       Wb ||     Uss|     |0  | 
        Input params(S,Y) must be np.matrix
        """
        uall, uss, Y, wag, gravlen, bsm= args
        lens = len(gravlen)

        uall = np.matrix(uall)
        Y = np.matrix(Y).T
        k = 0 
        w = np.ones(gravlen[0,0])/np.exp(x[0])
        W = np.diag(w)
        for k in range(1, lens):
            w = np.ones(gravlen[k,0])/np.exp(x[k])
            W = slg.block_diag(W, np.diag(w))
       
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))
        #len2 = int(len(uss)/lens) #多台重力仪观测不一致
        #bsm = 2
        len2 = gravlen[0,2] - bsm
        wb = np.ones(len2)/np.exp(x[lens])
        Wb = np.diag(wb)
        for k in range(1, lens):
            len2 = gravlen[k,2] - bsm #20191014
            wb = np.ones(len2)/np.exp(x[lens + k])
            Wb = slg.block_diag(Wb, np.diag(wb))

        Wall = slg.block_diag(W, Wb)

        Wall = np.matrix(Wall)
        Wb = np.matrix(Wb)

        u11 = uall.T*Wall*uall
        u22 = uss.T*Wb*uss

        l0,u0,v0 = np.linalg.svd(u11) #SVD分解
        l20,u20,v1 = np.linalg.svd(u22)
        u21 = np.abs(u20)
        #print(u21)
        u210 = u21[:-2*lens]

        f1 = sum(np.log(np.abs(u0)))
        f2 = sum(np.log(u210))
        X = v0.T*(np.mat(np.diag(1./u0)))*l0.T*uall.T*Wall*Y
        f3 = (uall*X-Y).T*Wall*(uall*X-Y) #min AT P A
        f0 = sum(np.log(1./np.diag(W)))
        f3 = f3[0,0]

        return f0 + f1 - f2 + f3
    @classmethod
    def result(cls, xopt, mat_list, glen):
        """forward result with the optimated weights
        """
        udd = mat_list[0]
        uss = mat_list[1]
        ubb = mat_list[2]
        urr = mat_list[3]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]
        #print(glen)

        m, n = uss.shape
        m1, n1 = urr.shape
        m2, n2 = uag.shape

        aa = np.hstack([urr ,udd])
        zb0 = np.zeros([m2, n])
        ab = np.hstack([uag[:,n:] ,zb0])
        zc0 = np.zeros([m, n1], dtype = float)
        ac = np.hstack([zc0 ,uss])

        aa = np.vstack([aa, ab])
        aa = np.vstack([aa, ac])

        bb = np.hstack([ubb, dag])
        bb0 = np.zeros(m, dtype = float)
        bb = np.hstack([bb, bb0])
        
        A = np.matrix(aa)
        b = np.matrix(bb).T
        lens = len(glen)
        k = 0 
        w = np.ones(glen[0,0])/np.exp(xopt[0])
        W = np.diag(w)
        for k in range(1, lens):
            w = np.ones(glen[k,0])/np.exp(xopt[k])
            W = slg.block_diag(W, np.diag(w))
       
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))
        #len2 = int(len(uss)/lens)
        #cls.bsm = 2
        len2 = glen[0,2] - cls.bsm #20191014
        wb = np.ones(len2)/np.exp(xopt[lens])
        Wb = np.diag(wb)
        for k in range(1, lens):
            len2 = glen[k,2] - cls.bsm
            wb = np.ones(len2)/np.exp(xopt[lens + k])
            Wb = slg.block_diag(Wb, np.diag(wb))

        Wall = slg.block_diag(W, Wb)
        Wall = np.matrix(Wall)
        #print(np.diag(Wall))
        #print(np.linalg.cond(A))
        #print(np.linalg.cond(Wall))

        return cls.forward(Wall, A, b)
    @classmethod
    def goadj(cls, mat_list, glen, xinit, dinit, method = 1, maxiter = 1000):
        """Go adjustment using the matrix list
           the inversed weights will be estimated.  
        """
        udd = mat_list[0]
        uss = mat_list[1]
        ubb = mat_list[2]
        urr = mat_list[3]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]

        m, n = uss.shape
        m1, n1 = urr.shape
        m2, n2 = uag.shape

        aa = np.hstack([urr ,udd])
        zb0 = np.zeros([m2, n])
        ab = np.hstack([uag[:,n:] ,zb0])
        zc0 = np.zeros([m, n1], dtype = float)
        ac = np.hstack([zc0 ,uss])
        #print(aa.shape, ac.shape)
        #print(glen)
        #print(len(glen))
        aa = np.vstack([aa, ab])
        aa = np.vstack([aa, ac])

        bb = np.hstack([ubb, dag])
        bb0 = np.zeros(m, dtype = float)
        bb = np.hstack([bb, bb0])

        x0 = np.ones(len(glen))*xinit**2
        d0 = np.ones(len(glen))*dinit**2
        x0 = np.log(np.hstack([x0, d0]))

        ua = np.matrix(aa)
        ub = np.matrix(bb)
        uss = np.matrix(uss)

        args = (ua, uss, ub, wag, glen, cls.bsm)
        return cls.optimization(cls.likelihood, x0, args, method, maxiter)


class Bayadj1(Bayadj):
    """Bayesian Adjustment using non-Linear drift model
    for estimated scale factor of relative gravimeter 
    A derived class of Bayadj
    the staticmethod and classmethod have been designed.
    the instance of object has not be needed.
    Properties:
    NONE             
    functions:
    - likelihood1: the user likelihood with the defined problem
    - goadj1 : to optimize the weights or hyperparameters
    - result1 : to forward the gravity values with optimized weights                
    """    
    @staticmethod
    def likelihood1(x, *args): #x,Uall,Uss,Y,gravlen,wag
        """给定超参数后，计算ABIC数值
        W       S         W Y       
        |Wrg       ||Urr  Udd| =   |Dbb|
        |   Wag    ||Uag     |     |dag|
        |       Wb ||     Uss|     |0  | 
        Input params(S,Y) must be np.matrix
        kstart is the index of fix meter scale factor
        """
        uall, uss, Y, wag, gravlen, Y2, bsm, kstart = args
        lens = len(gravlen)

        uall = np.matrix(uall)

        kk = 0
        for i in range(lens):   #scale factor of meter
            len3 = int(gravlen[i,0])
            for j in range(kstart,len3):
                sf = np.exp(x[-i])
                Y[:,kk] = Y[:,kk]*sf + Y2[:,kk]
                kk += 1

        Y = np.matrix(Y).T
        k = 0
        w = np.ones(gravlen[0,0])/np.exp(x[0])
        W = np.diag(w)
        for k in range(1, lens):
            w = np.ones(gravlen[k,0])/np.exp(x[k])
            W = slg.block_diag(W, np.diag(w))
       
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))
        #len2 = int(len(uss)/lens)
        #bsm = 2
        len2 = gravlen[0,2] - bsm #20191014
        wb = np.ones(len2)/np.exp(x[lens])
        Wb = np.diag(wb)
        for k in range(1, lens):
            len2 = gravlen[k,2] - bsm
            wb = np.ones(len2)/np.exp(x[lens + k])
            Wb = slg.block_diag(Wb, np.diag(wb))

        Wall = slg.block_diag(W, Wb)

        Wall = np.matrix(Wall)
        Wb = np.matrix(Wb)

        u11 = uall.T*Wall*uall
        u22 = uss.T*Wb*uss

        l0,u0,v0 = np.linalg.svd(u11) #SVD分解
        l20,u20,v1 = np.linalg.svd(u22)
        u21 = np.abs(u20)
        #print(u21)
        u210 = u21[:-2*lens]

        f1 = sum(np.log(np.abs(u0)))
        f2 = sum(np.log(u210))
        X = v0.T*(np.mat(np.diag(1./u0)))*l0.T*uall.T*Wall*Y
        f3 = (uall*X-Y).T*Wall*(uall*X-Y) #min AT P A
        f0 = sum(np.log(1./np.diag(W)))
        f3 = f3[0,0]

        return f0 + f1 - f2 + f3
    @classmethod
    def result1(cls, xopt, mat_list, glen, kstart):
        """forward result with the optimated weights
        """
        udd = mat_list[0]
        uss = mat_list[1]
        ubb = mat_list[2]
        urr = mat_list[3]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]
        ubb2 = mat_list[7]

        m, n = uss.shape
        m1, n1 = urr.shape
        m2, n2 = uag.shape

        aa = np.hstack([urr ,udd])
        zb0 = np.zeros([m2, n])
        ab = np.hstack([uag[:,n:] ,zb0])
        zc0 = np.zeros([m, n1], dtype = float)
        ac = np.hstack([zc0 ,uss])

        aa = np.vstack([aa, ab])
        aa = np.vstack([aa, ac])

        bb = np.hstack([ubb, dag])
        bb0 = np.zeros(m, dtype = float)
        bb = np.hstack([bb, bb0])
        lens = len(glen)
        kk = 0
        for i in range(lens):   #scale factor of meter
            len3 = int(glen[i,0])
            for j in range(kstart, len3):
                sf = np.exp(xopt[-i])
                bb[:,kk] = bb[:,kk]*sf + ubb2[:,kk]
                kk += 1

        A = np.matrix(aa)
        b = np.matrix(bb).T

        k = 0
        w = np.ones(glen[0,0])/np.exp(xopt[0])
        W = np.diag(w)
        for k in range(1, lens):
            w = np.ones(glen[k,0])/np.exp(xopt[k])
            W = slg.block_diag(W, np.diag(w))
       
        W = slg.block_diag(W, np.linalg.inv(np.diag(wag)))
        #len2 = int(len(uss)/lens)
        len2 = glen[0,2] - cls.bsm #20191014
        wb = np.ones(len2)/np.exp(xopt[lens])
        Wb = np.diag(wb)
        for k in range(1, lens):
            len2 = glen[k,2] - cls.bsm
            wb = np.ones(len2)/np.exp(xopt[lens + k])
            Wb = slg.block_diag(Wb, np.diag(wb))

        Wall = slg.block_diag(W, Wb)

        Wall = np.matrix(Wall)
        
        return cls.forward(Wall, A, b)
    @classmethod
    def goadj1(cls, mat_list, glen, xinit, dinit, sfinit, kstart = 0):
        """Go adjustment using the matrix list
           the inversed weights will be estimated.  
        """
        udd = mat_list[0]
        uss = mat_list[1]
        ubb = mat_list[2]
        urr = mat_list[3]
        uag = mat_list[4]
        dag = mat_list[5]
        wag = mat_list[6]
        ubb2 = mat_list[7]

        m, n = uss.shape
        m1, n1 = urr.shape
        m2, n2 = uag.shape

        aa = np.hstack([urr ,udd])
        zb0 = np.zeros([m2, n])
        ab = np.hstack([uag[:,n:] ,zb0])
        zc0 = np.zeros([m, n1], dtype = float)
        ac = np.hstack([zc0 ,uss])

        aa = np.vstack([aa, ab])
        aa = np.vstack([aa, ac])

        bb = np.hstack([ubb, dag])
        bb0 = np.zeros(m, dtype = float)
        bb = np.hstack([bb, bb0])

        x0 = np.ones(len(glen))*xinit**2
        d0 = np.ones(len(glen))*dinit**2
        sf0 = np.ones(len(glen))*sfinit**2

        x0 = np.log(np.hstack([x0, d0, sf0]))

        ua = np.matrix(aa)
        ub = np.matrix(bb)
        uss = np.matrix(uss)
        ub2 = np.matrix(ubb2)

        args = (ua, uss, ub, wag, glen, ub2, cls.bsm, kstart)
        return cls.optimization(cls.likelihood1, x0, args)

class Bayadj2(Bayadj):
    """Bayesian Adjustment using non-Linear drift model
       for multi-repetated gravity observations
    
    A derived class of Bayadj
    the staticmethod and classmethod have been designed.
    the instance of object has not be needed.
    Properties:
       NONE             
    functions:
    - likelihood2: the user likelihood with the defined problem
    - goadj2 : to optimize the weights or hyperparameters
    - result2 : to forward the gravity values with optimized weights                
    """    
    pass
 
if __name__ == '__main__':
    
    adj1 = Adjustment('adj test')
    print(adj1)
    matdata = adj1.load_json_mat('./data/gravmat.json')
    for key,values in matdata.items():
        print(key)
        np.asarray(matdata[key], dtype=np.float32)
