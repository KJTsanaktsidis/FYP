"""Class for interacting with simulation input data
"""

import numpy as np
import os.path
import scipy.interpolate as interp


class InputDatastore():

    def __init__(self, data_dir, prefix, T):
        """
        Initialiser for this datastore

        :param data_dir: Directory where the csv files are found
        :param prefix: Prefix for the filenames
        :param T: Temperature in kelvin
        """
        #we store 3 kinds of data
        diff_name = os.path.join(data_dir, prefix + '_Diffusivity_' + str(T) + 'K.csv')
        resist_name = os.path.join(data_dir, prefix + '_Resistivity_' + str(T) + 'K.csv')
        exper_name = os.path.join(data_dir, prefix + '_Experimental_' + str(T) + 'K.csv')

        self.diffusivity_raw = np.genfromtxt(diff_name, skip_header=1)
        self.resistivity_raw = np.genfromtxt(resist_name, skip_header=1)
        self.experimental_raw = np.genfromtxt(exper_name, skip_header=1)

        nSets = np.shape(self.experimental_raw)[1] - 1
        self.experimental_dict = dict()
        self.experimental_C = self.experimental_raw[1:, 0]
        for col_index in range(1, nSets + 1):
            current = self.experimental_raw[0, col_index]
            self.experimental_dict[current] = self.experimental_raw[1:, col_index]

    def interpolated_vector(self, d_type, size):
        if d_type == 'diffusivity':
            operative_vec = self.diffusivity_raw
        elif d_type == 'resistivity':
            operative_vec = self.resistivity_raw
        else:
            raise ValueError('Unknown d_type ' + d_type)

        spline = interp.InterpolatedUnivariateSpline(operative_vec[:, 0], operative_vec[:, 1])
        x = np.linspace(0, 1, num=size)
        return spline(x)

    def interpolated_diffusivity(self, size):
        return self.interpolated_vector('diffusivity', size)

    def interpolated_resistivity(self, size):
        return self.interpolated_vector('resistivity', size)