# mapFunctionOnStructure
Python package for multi-modal connectivity analysis

Sample code in python for model identification is provided. Model identification is the process we use to identify structural connections related to functional connections assigned with a probability. The precision matrix is estimated based on the Ledoit-Wolf estimator. This is a regularised version of estimating covariance matrices compared to simple correlation coefficient. (Main script: fun2strucInt_v01.py)

I also provide the scripts used to implement the bootstrap procedure in our IEEE TMI publication.
The Gaussian graphical model based on the Iterative Proportional Scaling algorithm is written by Dr. Gael Varoquaux. (Main script: fun2strucBootstrap_v01.py)

## Description of the Software
**Inputs**
- The structural connectivity matrices. \n
  /n/-by-/n/ arrays contained in text files within the folder specified. Where /n/ is the number of total regions examined.
- The average time-series of fMRI data. \n
  /n/-by-/vols/ array contained in text file within the folder specified. Where /vols/ is the number of fMRI volumes.
  Data should be preprocessed. Typically, we also remove the average signal of white matter, csf and motion parameters.
- The percentage of connections to be included in the model (controls sparsity). \n
  Typically, we use /model selection/ to choice this parameter. Model selection involves cross-validation to select the model with the least prediction error. In our experiments, a value of around 60\% to 70\% provided the best results. This parameter may differ according to the description of structural and functional connectivity as well as the noise in the data.
- For visualisation of the results:
  -- A subset of regions should be provided. \n
    For a comprehensive visualisation of the results, we recommend a small sub-set of regions.
  -- A threshold, /thres/. \n
    The structural connection with probability within the highest /thres/ (as percentage) will be displayed.
