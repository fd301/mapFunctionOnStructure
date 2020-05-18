"""
Use randomised lasso to select non-zero coefficients
Written by Gael Varoquaux: http://gael-varoquaux.info/
It requires sklearn: http://scikit-learn.org/stable/

"""

import numpy as np

from sklearn.linear_model import LassoLarsCV, RandomizedLasso
from sklearn.linear_model.randomized_l1 import _randomized_lasso

class MyRandomizedLasso(RandomizedLasso):

    def _make_estimator_and_params(self, X, y):
        alpha = self.alpha
        if alpha in ('auto', 'cv'):
            alpha_max = np.max(np.abs(np.dot(X.T, y))) / X.shape[0]
            if alpha == 'auto':
                alpha_min = .1 * alpha_max
            elif alpha == 'cv':
                alpha_min = LassoLarsCV(max_iter=self.max_iter,
                                        eps=self.eps,
                                        n_jobs=self.n_jobs,
                                        precompute=self.precompute).fit(
                                            X, y).alpha
            alpha = np.linspace(alpha_max, alpha_min, 10)
            return _randomized_lasso, dict(alpha=alpha,
                    max_iter=self.max_iter, eps=self.eps,
                    precompute=self.precompute)
        else:
            return RandomizedLasso._make_estimator_and_params(self, X, y)


if __name__ == '__main__':
    from sklearn import datasets
    diabetes = datasets.load_diabetes()
    X = diabetes.data
    y = diabetes.target
    model = MyRandomizedLasso(alpha='cv', n_jobs=-1)
    model.fit(X, y)
