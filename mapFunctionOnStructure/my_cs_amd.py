import numpy as np
from ctypes import *


def cs_amdW ( sp, order, myCSparceLibP ):  #/* order 0:natural, 1:Chol, 2:LU, 3:QR */
  """ Estimate the structural support.
  
  Parameters
  ----------
  sp: A 2D ndarray shape (n,n), n is the number of rois. This must be a symmetric matrix
  containing the structural connectivity matrix across subjects
  order: 0:natural ordering, 1:Chol, 2:LU, 3:QR
  myCSparceLibP: Path to the dynamically linked cs_amd library 
  
  Returns
  -------
  spOrder: 1D ndarray vector of shape n that contains the permutation of the matrix that 
  would provide a sparser cholesky decomposition
  
  Notes
  -----
  This function use ctypes to access the CSparse library:
  http://www.cise.ufl.edu/research/sparse/CSparse/
  Direct Methods for Sparse Linear Systems, T. A. Davis, SIAM, Philadelphia, Sept. 2006. 
  Part of the SIAM Book Series on the Fundamentals of Algorithms. 

  Note that the library has been build dynamically (instead of statically which is the default)
  !IMPORTANT NOTE! This code works for 64-bit machines only. 
  For 32-bit machines 'c_int64' has to be replaced with 'c_int32
  
  This my attempt to translate the cs_amd function to python so I can directly retrieve the ordering:  
  #define csi ptrdiff_t
  csi *cs_amd (csi order, const cs *A)
  typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
  {
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
  } cs ;  
  """
  
  n,m = sp.shape
  assert n == m, ('structural matrix must be square')
  lenR = n
  
  myCSparceLib = CDLL(myCSparceLibP)

  ri,ci = np.where(sp!=0)
  cmeth = c_int(order) #choose the AMD for the cholesky: 1
  #construct a structure to accommodate the support matrix appropriate to pass it to the function
  class cs(Structure): #define structure
    _fields_ = [("nzmax", c_int64),("m", c_int64),("n",c_int64),("p",POINTER(c_int64)),("i",POINTER(c_int64)),("x",POINTER(c_double)),("nz",c_int64)]

  #initialise structure: Note that myCind, myRind, mysp define a sparse array the same way that 
  # mxGetJc, mxGetIr, and mxGetPr functions work, respectively (mex functions in matlab).
  numIt = lenR*lenR
  numNZ = len(ri)
  myCindType = c_int64*(lenR+1)     #Need to be contiguous blocks of memory (the C way)
  myCind = myCindType()
  myRindType = c_int64*numNZ
  myRind = myRindType()
  myspType = c_double*numNZ
  mysp = myspType()

  count = 0
  for cc in range(0,lenR):
    countCE = 0
    for rr in range(0,lenR):
      if sp[rr,cc]:
	myRind[count] = rr
	mysp[count] = sp[rr,cc]
	count = count+1
	countCE = countCE+1
    if(cc==0):
      myCind[cc] = 0
    elif count==0:
      myCind[cc] = 0
    else:
      myCind[cc] = count-countCE
  
  myCind[lenR] = numNZ

  my_csStruc = cs( c_int64(numNZ),c_int64(lenR),c_int64(lenR),(myCind),(myRind),(mysp), c_int64(-1) ) #insert values to structure
  myCSparceLib.cs_amd.restype = POINTER(c_int64)		#specify the type of the return pointer so results are readable
  ordObj = myCSparceLib.cs_amd(cmeth,byref(my_csStruc) )  	#call function cs_amd via ctypes
  spOrder = np.empty([lenR],dtype=int)				#copy results in an np array
  for i in range(0,lenR):
    spOrder[i] = ordObj[i]
  
  #suspect that ordObj need to be freed here
  myCSparceLib.cs_free(ordObj)
  
  return spOrder
    
    