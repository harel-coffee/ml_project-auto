#! /usr/bin/env python

"""
Copyright of the program: Andrea Agazzi, UNIGE


"""

import numpy as np
import importd
import cv_sampling as cv
import matplotlib.pyplot as plt
import sklearn as sk

class Pipe(object):
    """
    Class for the successive application of different preprocessing/Learning
    methods for the classification and biomarker extraction for type II
    diabetic mice metabolomic data.
    ------
    Attributes:

        scoring

        steps

        estimator

        fit_params

        refit

        cv

        cv_scores

        verbose

    """

    def __init__(self, estimator, scoring=None,
                 fit_params=None, refit=True, cv=None, verbose=0 ):

        self.steps =
        self.estimator = estimator
        self.scoring = scoring
        self.fit_params = fit_params if fit_params is not None else {}
        self.refit = refit
        self.cv = cv
        self.cv_scores = []
        self.verbose = verbose
        self.bestrank = []

    def cv_grid(self, ):
        """
        Computes the cross-validation scores for the ranges of input parameters
        indicated in the inputs.

        The results are saved in the corresponding object attributes.
        """

    def _fit(self):
        """
        Fits the pipeline with the parameters given as input.
        """
        Xt, fit_params = self._pre_transform(X, y, **fit_params)
        self.steps[-1][-1].fit(Xt, y, **fit_params)
        return self

    def predict(self, X):
        """Applies transforms to the data, and the predict method of the
        final estimator. Valid only if the final estimator implements
        predict.
        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].predict(Xt)

    def score(self, X, y=None):
        """Applies transforms to the data, and the score method of the
        final estimator. Valid only if the final estimator implements
        score.
        Parameters
        ----------
        X : iterable
            Data to score. Must fulfill input requirements of first step of the
            pipeline.
        y : iterable, default=None
            Targets used for scoring. Must fulfill label requirements for all steps of
            the pipeline.
        """
        Xt = X
        for name, transform in self.steps[:-1]:
            Xt = transform.transform(Xt)
        return self.steps[-1][-1].score(Xt, y)

    def
