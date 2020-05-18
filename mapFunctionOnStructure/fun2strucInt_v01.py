import numpy as np
import pylab as pl
import itertools as it
from sklearn.covariance import LedoitWolf
#from ips import prop_scaling
from StrucSup import StrucSup
from my_cs_amd import cs_amdW
from run_randomisedLasso import run_randomisedLasso
from prepChol import prepChol
from visual_res import visual_res
#import pdb

#declare which functional connections to visualise
roisF = [2,0]
thresV = 0.2;

#percentage of connections to include
perCent = 60./100.

#subject ids
folder_names = [ '3188','3227','3295','3323','3333','3367','3371','3382',              '3388N','3397N','3409','3415','3422','3455','4504','4523','4524','4525','4526','4546','4563','4662',              '5185','5313','5587','5599','5020','5031','5067','5110','5142','5156','5193','5217','5301','5735','5821','5934','6125','6152','6160' ]
#folder_names = [ '3188','3227','3295' ]
subjNum = len(folder_names)

#vector with indeces to the rois selected
rois = range(0,18)
lenR = len(rois)

prefS=str(lenR)

#output main directory: intermediate results and final are written there
pathMainDir='../preprocessing/INRIA/subj41_RSNsG'

#path to structural connectivity matrices
pathDTI = pathMainDir + '/' + 'StrucCons' + '/'

#path to fmri timeseries
pathfMRI = pathMainDir + '/' + 'timeseries_clean' + '/'

#path to library cs_amd library
myCSparceLibP = "./libctest.so.1.0"

#load structural connectivity matrices
print "Loading structural data ..."
scAll = np.empty([lenR,lenR,subjNum])
for ii,subj in enumerate(folder_names):
  inputP =  pathDTI + 'con_' + subj + '.txt'
  sc = np.loadtxt(inputP)
  sc = sc[rois,]
  sc = sc[:,rois]
  sc += sc.T  		#connectivity matrix need to be symmetrised: here I use the average
  sc = sc/2.0
  scAll[:,:,ii] = sc 

#run the script that estimates structural support 
#given the percentage of desirable connections included
print "Estimasting structural support ..."
sp,thres = StrucSup(perCent,scAll)

#provide the ordering that would result in the sparser cholesky
order = 1		#0:natural ordering, 1:Chol, 2:LU, 3:QR
print "Estimasting AMD ordering ..."
spOrder = cs_amdW(sp, order, myCSparceLibP)

#estimate precision matrix from either Gaussian graphical model or
#the LedoitWolf regularisation (here we are using the second)
#note that with the typical correlation coefficients the estimation
#of the covariance matrix is not well regularised especially with a large number of rois.
#This problem is exacerbated with the inversion of the correlation matrix 
print "Loading fMRI averaged signal ..."
print "Estimate the precision matrix ..."
precAll = np.empty([lenR,lenR,subjNum])
for ii,subj in enumerate(folder_names):
  inputP = pathfMRI + 'timeSeriesR' + prefS + '_' + subj + '.txt'
  signal = np.loadtxt(inputP)
  LedW = LedoitWolf(store_precision=True, assume_centered=False)
  LedW.fit(signal.T)
  prec = LedW.get_precision()   #estimate the precision matrix based on the LedoitWolf regularisation
  precAll[:,:,ii] = prec
    
print "Preparing Data for randomised Lasso: Cholesky decomposition, scaling diagonal and so on..."    
#prepare data for ransomised Lasso
Y,X,ri0,ci0 = prepChol(scAll,precAll,sp,spOrder)

print "Randomised Lasso is running ..."
prob_all = run_randomisedLasso(X.T,Y)

print "Reorder data and visualise ..."
adj2,adj = visual_res(prob_all,lenR,roisF,ri0,ci0,thres)

pl.matshow(adj2)
pl.colorbar()
pl.show()

print "Good bye!"

