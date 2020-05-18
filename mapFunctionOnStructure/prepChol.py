import numpy as np
#import pylab as pl


def prepChol(scAll,precAll,sp,spOrder):
  """ Estimate the cholesky decomposition and prepare data for randomised lasso.
  
  Parameters
  ----------
  scAll: A 3D ndarray shape (n,n,p) containing all the structural connectivity matries.
	  n is the number of rois and p is the number of subjects
  precAll: A 3D ndarray shape (n,n, p), cotaining all the precision matrices
	  n is the number of rois and p is the number of subjects
  sp: A 2D ndarray with shape( n,n) containing the structural support over all subjects
  spOrder: A 1D ndarray with shape n containing the ordering of the support matrix
  rois: labels of rois
  
  Returns
  -------
  precAll: A 2D ndarray shape (N,p), where N is the number of connections and p is the number of subjects.
  It contains the functional connections  across all subjects.
  scAll: A 2D ndarray shape (N,p), where N is the number of connections and p is the number of subjects.
  It contains the structural connections  across all subjects.
  ri0: row indices that correspond to N element-by-element 
  ci0: column indices that correspond to N element-by-element
      Therefore, N[i] corresponds to the connection between regions: ri0[i] and ci0[i]
  
  Notes
  Cholesky decomposition according to AMD ordering.
  -----
  """
  
  n, m, p = scAll.shape
  n2,m2,p2 = precAll.shape
  n3,m3 = sp.shape
  n4 = len(spOrder)
  assert n == m, ('structural matrix must be square')
  assert n2 == m2, ('precision matrix must be square')
  assert n3 == m3, ('support matrix must be square')
  assert n4 ==n, ('vector with order must have same number rois as structural matrix')
  assert n == n2 and p == p2, ('structural matrix and precision matrix must have same dimensions')
  assert n == n3, ('support matrix should have same dimensions with the structural matrix')

  lenR = n
  
 
    
  totCons = lenR*lenR
  sp = sp[spOrder,:]
  sp = sp[:,spOrder]
  scAll = scAll[spOrder,:,:]
  scAll = scAll[:,spOrder,:]
  precAll = precAll[spOrder,:,:]
  precAll = precAll[:,spOrder,:]


  precDiag = np.empty([p,lenR])
  for ii in range(0,p):
    precDiag[ii,:] = np.diag(precAll[:,:,ii])
    InvdiagM = np.zeros([lenR,lenR])
    for j in range (0,lenR):
      InvdiagM[j,j] = 1.0/precDiag[ii,j]  
    precAll[:,:,ii] = np.linalg.cholesky(precAll[:,:,ii]) 	 	#apply cholesky
    precAll[:,:,ii] = np.dot( InvdiagM,  precAll[:,:,ii] ) 	 #use diagonal to scale cholesky decomposition

  sp = np.tril(sp)			#elements above the main diagonal will be zeroed 
  sp = sp.reshape([totCons])
  scAll = scAll.reshape([totCons,p])
  precAll = precAll.reshape([totCons,p])



  ci = np.where(sp>0) 
  ci = ci[0]
  scAll = scAll[ci,:]			#apply support to structural data
  precAll = precAll[ci,:]			#apply support tp functional data

  #recover the correct indices
  rois = range(0,lenR)
  rois2 = np.asarray(rois)
  rois2 = rois2.reshape(lenR,1)
  tt = np.ones([1,lenR])
  rii = np.dot(rois2,tt)			#rows' indices
  cii = rii.T				#columns' indices

  rii = rii[spOrder,:]
  cii = cii[:,spOrder]
  rii = rii.reshape([totCons])
  cii = cii.reshape([totCons])
  ri0 = rii[ci]     	#this must be the row indices that correspond to prob_all element by element
  ci0 = cii[ci]	    	#this must be the column indices that correspond to prob_all element by element	
 
  return precAll, scAll, ri0,ci0
