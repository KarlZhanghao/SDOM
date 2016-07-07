#!/opt/local/bin/python
#****************************************************************************
#*                                                                          *
#*                            OM_SPoD - Algorithm                           *
#*                                                                          *          
#****************************************************************************
# This program is modified from 'Supplementary Software' of [1]

#Referrence:
#[1]Hafi, N. et al. Fluorescence nanoscopy by polarization modulation and polarization angle narrowing. Nature Methods 11, 579-584 (2014).
#[2]Beck A, Teboulle M. A fast iterative shrinkage-thresholding algorithm for linear inverse problems[J]. SIAM journal on imaging sciences, 2009, 2(1): 183-202.

import numpy as np
import scipy as sp
import scipy.fftpack
import scipy.misc
import scipy.signal
import scipy.io as sio
import math
import os
import sys
import matplotlib.pyplot as plt
#*****************************************************************************************
#                                      Function
#*****************************************************************************************
#Discrete Cosine Transform
def ctf(a): 
    return sp.fftpack.dct(sp.fftpack.dct(a, axis=1, norm='ortho'), axis=0, norm='ortho')

#Inverse Discrete Cosine Transform
def ictf(a): 
    return sp.fftpack.idct(sp.fftpack.idct(a, axis=1, norm='ortho'), axis=0, norm='ortho')

#The operator p_L() in [2]
def step(L, lamb, y, grad): 
    c = (y[0] - grad[0]/L, y[1] - grad[1]/L)
    tmp1 = np.maximum(c[0] - lamb[0]/L, 0)
    tmp2 = np.minimum(c[0] + lamb[0]/L, 0)
    tmp3 = np.maximum(c[1] - lamb[1]/L, 0)
    tmp4 = np.minimum(c[1] + lamb[1]/L, 0)
    return (np.maximum(tmp1 + tmp2, 0), tmp3 + tmp4)

#Possion parameter \mu
def forward(tup, psf, innerRect, mod): 
    g = tup[0]
    period = g.shape[2]
    tmp1 = np.zeros_like(g)
    for j in range(period):
        tmp1[:,:,j] = mod[j]*sp.signal.fftconvolve( g[:,:,j], psf, 'same' )
    return np.maximum(tmp1[innerRect] + ictf(tup[1])[:,:,np.newaxis]*mod, 1e-6)

def maximumLikelihood(forwardProjection, I):
    return (forwardProjection - I*np.log(forwardProjection)).sum()

def quadraticApprox(gx, gy, grady, L):
    delta0 = gx[0] - gy[0]
    delta1 = gx[1] - gy[1]
    return (delta0*grady[0] + L/2*delta0*delta0).sum() + (delta1*grady[1] + L/2*delta1*delta1).sum()

def backward(h, psf, mod):
    period = h.shape[2]
    tmp0 = h
    tmp1 = np.zeros_like(g)
    tmp2 = np.zeros(h.shape[0:2])
    for j in range(period):
        tmp1[:,:,j] = mod[j]*sp.signal.fftconvolve(tmp0[:,:,j], psf, 'full' )
        tmp2 = tmp2 + mod[j]*ctf(h[:,:,j])
    return (tmp1, tmp2)

def gradient(forwardProjection, I, psf, mod): 
    tmp = 1 - I/forwardProjection
    return backward(tmp, psf, mod)

# Solution algorithm for optimization model : Figure the following Notations by reference to 'FISTA with backtracking' in [2] 
def oneIteration(i):
    global xkm1, xk, yk, ykp1, tk, tkp1, L, lamb, innerRect, minimumLikelihood, oldFileList, psf, back, a, modulation
    xkm1 = xk
    yk = ykp1
    tk = tkp1
    forwardProjectionY = forward(yk, psf, innerRect, modulation) # \mu induced by yk
    maxLikelihoodY = maximumLikelihood(forwardProjectionY, a)  # f(yk) 
    grad = gradient(forwardProjectionY, a, psf, modulation)   # grade of f
    for j in range(1000):  # search L
        Ltest = L*math.pow(1.1, j)  
        xtest = step(Ltest, lamb, yk, grad) # xtest is p_L(y_k) and step is the operator p_L()
        newForwardProjection = forward(xtest, psf, innerRect, modulation) # new \mu by xtest = p_L(y_k) 
        newMaxLikelihood = maximumLikelihood(newForwardProjection, a)  # f(p_L(y_k))
        quadratic = maxLikelihoodY + quadraticApprox(xtest, yk, grad, Ltest) # The first three items of QL() 
        if j > 0:
            print "Ltest =", Ltest, "newMaxLikelihood =", newMaxLikelihood, "quadratic =", quadratic,  "difference =", newMaxLikelihood - quadratic
        if newMaxLikelihood <= quadratic:   # inequality condition F(pL(yk))<QL(pL(yk),yk)
            xk = xtest
            L = Ltest
            maxLikelihood = newMaxLikelihood
            break
    else:
        print "Lipshitz factor still too small after 1000 tries."
        exit()
    tkp1 = (1 + np.sqrt(1 + 4*tk*tk)/2)
    over = (tk - 1)/tkp1
    ykp1 = (xk[0] + over*(xk[0] - xkm1[0]), xk[1] + over*(xk[1] - xkm1[1]))    
#    print 'iter=', i
#    if i == 1200:
#        sio.savemat( dataDir+'tmp/' + 'ReconResult1.mat', {'data': xk[0][innerRect]})
#    if i%10==0 and i>0:
#        newFileList = []
#        fileName = dataDir+'tmp/' + 'output_'+str(i)+'_bg'
#        np.savetxt(fileName, ictf(xk[1]))
#        newFileList.append(fileName)
#        for j in range(period):
#            fileName = dataDir+'tmp/' + 'output_'+str(i)+'_'+str(j)
#            np.savetxt(fileName,  xk[0][:,:,j])
#            newFileList.append(fileName)
#        for j in range(len(oldFileList)):
#            os.remove(oldFileList[j])
## the output_?_? files can be opened by  text
#        oldFileList = newFileList
#        sp.misc.imsave(dataDir+'tmp/' + 'SR.png', xk[0].sum(2)[innerRect])
#        for k in range(period):
#             sp.misc.imsave(dataDir+'tmp/' + str(k)+'.png', xk[0][:,:,k][innerRect])
#        sp.misc.imsave(dataDir+'tmp/' + 'bg.png',ictf(xk[1]))
#        sio.savemat( dataDir+'tmp/' + 'g.mat', {'g': xk[0][innerRect]})
#        sio.savemat( dataDir+'tmp/' + 'b.mat', {'b': ictf(xk[1])})
#*****************************************************************************************
#                                    Main Function
#***************************************************************************************** 

#************************************setting parameters***************************************  
dataDir = sys.argv[1]
#dataDir = 'K:/_PaperSubmit/Correspondence_nmeth/manuscript_submit_v3/SupplementarySoftware/Data/NeuronalSpine2/'
data = sio.loadmat( dataDir + 'tmp/ReconData.mat')
data = data['data']
reconPara = data['ReconPara'][0][0]
lamb = (reconPara[0][0], reconPara[1][0])
#lamb = (0.075, 10)
iterations = np.int(reconPara[2][0])
#iterations = 1201
L = 0.001
#----------------------------read the observed images I(r,t)---------------------------------
a = data['image'][0][0]
a = np.float64( a)
period = a.shape[2]
avg = a.sum((0,1)) / a.shape[0] / a.shape[1]
imgShape = a.shape[0:2]
#info = os.listdir( dataDir+'WideField/')
#period = len(info)
#imgTmp = plt.imread( dataDir+'WideField/'+info[0])
#imgShape = imgTmp.shape
#a = np.empty(imgTmp.shape+(period,), dtype="float64") # a is I(r,t)
#avg = np.empty((period,), dtype="float64")
#for i in range(period):
#    a[:,:,i] = np.rint(plt.imread(dataDir+'WideField/'+info[i]).astype("float64"))
#    sp.misc.imsave( dataDir +  'tmp/original_'+ str(i)+'.png', a[:,:,i])
#    avg[i] = np.sum(a[:,:,i])/a[:,:,i].size
#sp.misc.imsave( dataDir  +  'tmp/original.png', a.sum(2))
# modulation is the ratio of average intensity,namely,I(t)
modulation = avg/avg[0]
#--------------------------read psf(point spread function) files and reshape------------------    
#psfShape = ( 129, 129)
#psf = np.fromfile(dataDir + 'tmp/figure_4_psf.raw', dtype="float32").astype("float64").reshape(psfShape)
#psf /= psf.sum()
psf = data['psf'][0][0];
psf = psf/psf.sum()
psfShape = psf.shape
offset=(psfShape[0]/2, psfShape[1]/2)
upImgShape = (imgShape[0], imgShape[1])
innerRect=(slice(offset[0], offset[0]+upImgShape[0]), slice(offset[1], offset[1]+upImgShape[1]))
extendedShape = (upImgShape[0]+psfShape[0]-1, upImgShape[1]+psfShape[1]-1)
#initial value: g is in spatial domain and 0 initial value. b is in frequency domain and just not zero in (0,0)
average = np.sum(a)/a.size
b = np.zeros(imgShape, dtype="float64")
b[0, 0] = average*np.sqrt(imgShape[0]*imgShape[1])
g = np.zeros(extendedShape + (period, ))
xk = (g, b)
ykp1 = (g, b)
tkp1 = 1.0
oldFileList = []
#----------------------------FISTA with backtracking algorithm in [2]---------------------------------
for i in range(iterations):
    oneIteration(i)
#----------------------------save results---------------------------------
sio.savemat( dataDir+'tmp/' + 'ReconResult.mat', {'data': xk[0][innerRect]})


