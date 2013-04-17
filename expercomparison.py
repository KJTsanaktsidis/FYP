"""Module for shifting and comparing model + experimental results
"""

import numpy as np
from calcsim import CalcSimWrapper


class ComparisonEngine():

    def __init__(self, cs):
        #need a calc sim wrapper object
        assert isinstance(cs, CalcSimWrapper)
        self.cs = cs

    def calibrate(self, model, experiment):
        """
        This method will attempt to align model to experiment by left-padding with ones and right padding with zeros
        and then computing a cross-correlation
        The optimum shift is saved in this object, and (least-squares, shift) is returned

        :param model: The vector of model data
        :param experiment: the vector of experimental data
        """

        if len(model) != len(experiment):
            raise ValueError("Model and experiment must be the same length")

        self.model = model
        self.shift = self.cs.fast_pad_shift(experiment, model)
        return self.lsq(model, experiment), self.shift

    def shifted_lsq(self, model, experiment):
        """
        This method will shift model by the pre-saved shift factor and then return the least squares error betwen
        the two.

        :param model: The vector of model data
        :param experiment: the vector of experimental data
        """
        if self.shift > 0:
            #right shift
            keep = model[:-self.shift]
            add = np.ones(self.shift)
            shifted = np.concatenate((add, keep))
        elif self.shift < 0:
            keep = model[-self.shift:]
            add = np.zeros(-self.shift)
            shifted = np.concatenate((keep, add))
        else:
            shifted = model

        return self.lsq(shifted, experiment)

    def lsq(self, y1, y2):
        """
        This method computes the sum of the squares of the differences between y1 and y2
        """
        return np.sum((y2 - y1) ** 2)