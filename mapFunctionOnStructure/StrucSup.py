import numpy as np
from scipy import stats as stats 
#import pylab as pl


def StrucSup(perCent,strucCon):
  """ Estimate the structural support.
  
  Parameters
  ----------
  perCent: a scalar number of percentage of connections to be included (from 0 to 1)
  strucCon: A 3D ndarray shape (n,n, p), n is the number of rois and p is the number of subjects
  containing the structural connectivity matrix for each subject
  
  Returns
  -------
  struc_sup: 2D ndarray
  the structural support with ones where a connection is supported and zero for the rest
  thresh: float
  a p-value that corresponds to the threshold of the statistical significance test that would 
  include this percentage of connections
  
  Notes
  -----
  The 
  """
  
  n, m, p = strucCon.shape
  totConA = n*m 			#total number of connections
  numConS = round(perCent*totConA)	#number of connections to be selected
  t1,p1 = stats.ttest_1samp(strucCon,popmean=0,axis=2)
  p2 = p1.reshape([totConA])
  p2 = p1.copy()
  p2 = p2.reshape([totConA])
  p2.sort()
  thres = p2[numConS]
  ri,ci = np.where(p1>=thres)
  struc_sup = np.ones([n,m])
  struc_sup[ri,ci] = 0.
  for i in range(0,n):  		#put one at diagonal
    struc_sup[i,i] = 1.
  
  return struc_sup, thres