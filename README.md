# mapFunctionOnStructure
Python package for multi-modal connectivity analysis

Sample code in python for model identification is provided. Model identification is the process we use to identify structural connections related to functional connections assigned with a probability. The precision matrix is estimated based on the Ledoit-Wolf estimator. This is a regularised version of estimating covariance matrices compared to simple correlation coefficient. (Main script: fun2strucInt_v01.py)

I also provide the scripts used to implement the bootstrap procedure in our IEEE TMI publication.
The Gaussian graphical model based on the Iterative Proportional Scaling algorithm is written by Dr. Gael Varoquaux. (Main script: fun2strucBootstrap_v01.py)

## Description of the Software
**Inputs:**
- The structural connectivity matrices. 
  n-by-n arrays contained in text files within the folder specified. Where n is the number of total regions examined.
- The average time-series of fMRI data. 
  n-by-vols array contained in text file within the folder specified. Where vols is the number of fMRI volumes.
  Data should be preprocessed. Typically, we also remove the average signal of white matter, csf and motion parameters.
- The percentage of connections to be included in the model (controls sparsity). 
  Typically, we use model selection to choice this parameter. Model selection involves cross-validation to select the model with the least prediction error. In our experiments, a value of around 60% to 70% provided the best results. This parameter may differ according to the description of structural and functional connectivity as well as the noise in the data.
- For visualisation of the results:
  - A subset of regions should be provided. 
    For a comprehensive visualisation of the results, we recommend a small sub-set of regions.
  - A threshold, thres. 
    The structural connection with probability within the highest thres (as percentage) will be displayed.

**Output of fun2strucInt_v01.py:**
- An adjacency matrix that displays the structural connections related to the functional connections estimated from the subset of regions provided.

**Output of fun2strucBootstrap_v01.py:**
- For each bootstrap iteration the probabilities for all the connections are stored along with the row and column indices to allow reconstructing the n-by-n-by-n*n matrix. This matrix provides the probabilities of each structural connection to be relate with each functional connection.

**Functions:**
The software includes the following functions:
- *fun2strucInt_v01.py:*
  Main script to run. All the parameters and paths are set here.
  - Load structural connectivity matrices and symmetrizes them.
  - Load functional time series and estimate the precision matrix based on the LedoitWolf estimator.
  - It calls: /StrucSup/, /cs_amdW/, /prepChol/, /visual_res/.  
- fun2strucBootstrap_v01.py:
  Main script to run the bootstrap algorithm to estimate confidence intervals. All the parameters and paths are set here.
  - Load structural connectivity matrices and symmetrizes them.
  - Load functional time series and estimate the precision matrix based on the LedoitWolf estimator.
  - It calls: /StrucSup/, /cs_amdW/, /prepChol/, /prop_scaling/. 
- StrucSup: 
  Estimate the structural support. This is a square array /n/-by-/n/ with with ones where a connection is supported and zero for the rest. 
- *cs_amdW:* 
  Estimate the /Minimum degree ordering/ that provide a sparser cholesky decomposition. Note that this is a function that wraps [http://www.cise.ufl.edu/research/sparse/CSparse/CSparse/Source/cs_amd.c cs_amd.c] function of the [http://www.cise.ufl.edu/research/sparse/CSparse/ CSparse library], directly.
  -- *libctest.so.1.0:* This is a compiled version of [http://www.cise.ufl.edu/research/sparse/CSparse/ CSparse library] on /Ubuntu 10.04.1 LTS/ (64-bit machine), linked dynamically. If you have a different system you need to rebuild the library. In addition, if your system is a 32-bit machine, you will need to modify /my_cs_amd.py/ to use /c_int32/.
- *prop_scaling:* 
  Estimate the precision matrix based on the Iterative Proportional Scaling (IPS) and the structural support.  
- *prepChol:* 
  Estimate the cholesky decomposition based on the ordering provided. It also reshapes the data to be /N/-by-/p/ arrays. Where /N/ is the total non-zero number of connections and /p/ is the number of subjects.
- *run_randomisedLasso:* 
  It runs randomised Lasso and it returns a list of elements. Each element corresponds to a functional connection and it is a list with * *N* * values. Each value is the probability the underlying structural connection to be selected. 
- *visual_res:*
  Visualise the results for a number of functional connections picked. Note that if a structural connection is picked in association with more than one functional connection, its probability will be averaged. Once randomised Lasso has finished. You can call this function independently to visualise the results for any functional connection. 

