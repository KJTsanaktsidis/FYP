"""Class for interacting with simulation input data
"""

import numpy as np
import os.path
import scipy.interpolate as interp


class InputDatastore():

    def __init__(self, data_dir, prefix):
        """
        Initialiser for this datastore

        :param data_dir: Directory where the csv files are found
        :param prefix: Prefix for the filenames
        :param Tfloat: Temperature in kelvin
        """

        #TODO: Make these do a temperature extrapolation instead of assume 973K

        #we store 3 kinds of data
        diff_name = os.path.join(data_dir, prefix + '_Diffusivity_973K.csv')
        resist_name = os.path.join(data_dir, prefix + '_Resistivity_973K.csv')
        fw_exper_name = os.path.join(data_dir, str.format('{}_Experimental_{}_{}K.csv', prefix, 'forward', '973'))
        rv_exper_name = os.path.join(data_dir, str.format('{}_Experimental_{}_{}K.csv', prefix, 'reverse', '973'))

        self.diffusivity_raw = np.genfromtxt(diff_name, skip_header=1)
        self.resistivity_raw = np.genfromtxt(resist_name, skip_header=1)
        self.experimental_fw_raw = np.genfromtxt(fw_exper_name, skip_header=0, delimiter=',')
        self.experimental_rv_raw = np.genfromtxt(rv_exper_name, skip_header=0, delimiter=',')

        #xaxis assumed to be the same between fw and rv
        self.experimental_x = self.experimental_rv_raw[1:, 0]

        self.experimental_fw_dict = {}
        self.experimental_rv_dict = {}

        for raw, edict in ((self.experimental_fw_raw, self.experimental_fw_dict),
                          (self.experimental_rv_raw, self.experimental_rv_dict)):
            nSets = np.shape(raw)[1] - 1
            for col_index in range(1, nSets + 1):
                current = int(np.abs(raw[0, col_index]))
                edict[current] = raw[1:, col_index]

    def interpolated_vector(self, operative_vec, size):

        spline = interp.InterpolatedUnivariateSpline(operative_vec[:, 0], operative_vec[:, 1])
        x = np.linspace(0, 1, num=size)
        return spline(x)

    def interpolated_diffusivity(self, size, T):
        return self.interpolated_vector(self.diffusivity_raw, size)

    def interpolated_resistivity(self, size, T):
        return self.interpolated_vector(self.resistivity_raw, size)

    def interpolated_experiment(self, current, x, edict):
        edict = self.edict_for_direction(edict)
        data = edict[current]
        spline = interp.InterpolatedUnivariateSpline(self.experimental_x, data)
        return spline(x)

    def interpolated_experiment_dict(self, x, edict):
        edict = self.edict_for_direction(edict)
        newdict = dict()
        for key in edict.keys():
            newdict[key] = self.interpolated_experiment(key, x, edict)
        return newdict

    def interpolated_experiment_fw_dict(self, x):
        return self.interpolated_experiment_dict(x, self.experimental_fw_dict)

    def interpolated_experiment_rv_dict(self, x):
        return self.interpolated_experiment_dict(x, self.experimental_rv_dict)

    def edict_for_direction(self, edict):
        if isinstance(edict, str):
            if edict == 'forward':
                return self.experimental_fw_dict
            elif edict == 'reverse':
                return self.experimental_rv_dict
            else:
                raise ValueError('Unknown direction ' + edict)
        elif isinstance(edict, dict):
            return edict
        else:
            raise ValueError('Unknown type of edict')