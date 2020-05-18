import numpy as np
import itertools as it
#import pylab as pl


def visual_res(prob_all,lenR,roisF,ri0,ci0,thres):
  """ Estimate the cholesky decomposition and prepare data for randomised lasso.
  
  Parameters
  ----------
  prob_all: A list with the probabilities that each connection is selected
  lenR: Number of total regions
  roisF: A list of regions. The results associated with each pair of regions will be examined.
  ri0: Rows' indices that correspond to N element-by-element 
  ci0: Columns' indices that correspond to N element-by-element
      Therefore, prob_all[i] corresponds to the connection between regions: ri0[i] and ci0[i]
  thres: A total of thres*total nonzero connections will be displayed    

  Returns
  -------
  adj2: Thresholded adjacency matrix
  adj: Total adjacency matrix for the functional connection specified in roisF
  
  Notes
  Cholesky decomposition according to AMD ordering.
  -----
  """
  
  CombCons=list(it.combinations(roisF,2))

  consI = np.empty([len(ri0),2])
  consI[:,0] = ri0
  consI[:,1] = ci0
  consI = np.sort(consI,1)     

  lenComb = len(CombCons)
  index0 = []
  prob0 = []
  adj = np.zeros([lenR,lenR])
  countA = np.zeros([lenR,lenR])
  for combi in CombCons:
    combs = np.sort(combi)
    print(combs)
    index = np.where( (consI[:,0]==combs[0])&(consI[:,1]==combs[1]) )
    index = index[0]
    prob = prob_all[index]		#structural connections that correspond to the functional: 'combi'
    index2 = np.where(prob!=0)
    index2 = index2[0]
    sCon = [ ri0[index2], ci0[index2] ]   #rows and columns that correspond to
    for jj in index2:
      adj[ ri0[jj],ci0[jj] ] = adj[ ri0[jj],ci0[jj] ] + prob[jj]
      adj[ ci0[jj],ri0[jj] ] = adj[ ri0[jj],ci0[jj] ] 
      countA[ri0[jj],ci0[jj] ] = countA[ri0[jj],ci0[jj] ] + 1.0
      countA[ci0[jj],ri0[jj] ] =  countA[ri0[jj],ci0[jj] ]
      
  in0 = np.nonzero(countA.flat)
  adj.flat[in0] = adj.flat[in0]/countA.flat[in0]

  nzLen = len(np.nonzero(adj.flat)[0])
  numC = np.rint(nzLen*thres)

  tt = np.sort(adj.flat)[-numC]
  ind2 = np.where(adj.flat<tt)
  ind2 = ind2[0]
  adj2 = adj
  adj2.flat[ind2] = 0
 
  return adj2, adj
