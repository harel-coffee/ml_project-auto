#! /usr/bin/env python

"""
Copyright of the program: Andrea Agazzi, UNIGE
Project started the 19.02.2016

Module containing crossvalidation routines for the main program
"""

import numpy as np

def leave_x_out(nlst,x,nsamples=300):
	""" Returns a random sample of size "samplesize" of a leave-x-out crossvalidation sampling"""

	# enumerate the predictors
	preds = range(len(nlst))
	outlst = []

	# sample nsamples times a new cv sample of size 112-x
	while len(outlst) < nsamples:
		bsample = np.random.choice(preds,size=len(nlst)-x,replace=False)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))

	return outlst

def l20o(nlst):
	""" Returns a random sample of size "samplesize" of a leave-x-out crossvalidation sampling"""

	# enumerate the predictors
	preds = range(len(nlst))
	outlst = []

	# sample nsamples times a new cv sample of size 112-x
	while len(outlst) < 300:
		bsample = np.random.choice(preds,size=92,replace=False)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))

	return outlst

def looo(nlst):
	""" Returns a random sample of size "samplesize" of a leave-x-out crossvalidation sampling"""

	# enumerate the predictors
	preds = range(len(nlst))
	outlst = []

	# sample nsamples times a new cv sample of size 112-x
	while len(outlst) < 100:
		bsample = np.random.choice(preds,size=52,replace=False)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))

	return outlst



"""
from sklearn import linear_model, decomposition, datasets
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV
"""

"""
logreg = sk.linear_model.LinearDiscriminantAnalysis()

pca = decomposition.PCA()
pipe = Pipeline(steps=[('pca', pca), ('logistic', logistic)])

digits = datasets.load_digits()
X_digits = digits.data
y_digits = digits.target

###############################################################################
# Plot the PCA spectrum
pca.fit(X_digits)

plt.figure(1, figsize=(4, 3))
plt.clf()
plt.axes([.2, .2, .7, .7])
plt.plot(pca.explained_variance_, linewidth=2)
plt.axis('tight')
plt.xlabel('n_components')
plt.ylabel('explained_variance_')

###############################################################################



# Prediction

n_components = [20, 40, 64]
Cs = np.logspace(-4, 4, 3)
estimator = GridSearchCV(pipe,dict(pca__n_components=n_components,logistic__C=Cs))
estimator.fit(X_digits, y_digits)

plt.axvline(estimator.best_estimator_.named_steps['pca'].n_components,linestyle=':', label='n_components chosen')
plt.legend(prop=dict(size=12))
plt.show()
"""
