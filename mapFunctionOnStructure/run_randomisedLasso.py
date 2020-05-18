"""
===================================================
Model Identification 
===================================================

"""
print __doc__

import numpy as np
import my_randomized_lasso as ms

def run_randomisedLasso(X,Y):
  """ Estimate the cholesky decomposition and prepare data for randomised lasso.
  
  Parameters
  ----------
  X: A 2D ndarray shape (N,p), where N is the number of connections and p is the number of subjects.
  It contains the structural connections  across all subjects.
  Y: A 2D ndarray shape (N,p), where N is the number of connections and p is the number of subjects.
  It contains the functional connections  across all subjects.
  
  Returns
  -------
  prob_all: A list with the probabilities that each connection is selected
  
  Notes
  -----
  It uses my_randomised_lasso, written by Gael Varoquaux: http://gael-varoquaux.info/
  which requires sklearn: http://scikit-learn.org/stable/
  """

  prob_all = list()
  alphas = list()
  for i_y, y in enumerate(Y):
    print(i_y)
    randLasso = ms.MyRandomizedLasso(alpha='cv',n_jobs=1)
    randLasso.fit(X,y)

    prob = randLasso.scores_
    prob_all.append(prob)
  
  return prob_all
