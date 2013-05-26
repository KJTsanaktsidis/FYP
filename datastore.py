"""Class for interacting with simulation input data
"""

import numpy as np
import os.path
import scipy.interpolate as interp
import functools


class InputDatastore():

    def __init__(self, data_dir, prefix):
        """
        Initialiser for this datastore

        :param data_dir: Directory where the csv files are found
        :param prefix: Prefix for the filenames
        :param Tfloat: Temperature in kelvin
        """

        #we store 3 kinds of data
        precise_diff_name = os.path.join(data_dir, prefix + '_Diffusivity_973K.csv')
        resist_name = os.path.join(data_dir, prefix + '_Resistivity_973K.csv')
        diff_arrhenius_name = os.path.join(data_dir, prefix + '_Diffusivity_Arrhenius_Constants.csv')
        fw_exper_name = os.path.join(data_dir, str.format('{}_Experimental_{}_{}K.csv', prefix, 'forward', '973'))
        rv_exper_name = os.path.join(data_dir, str.format('{}_Experimental_{}_{}K.csv', prefix, 'reverse', '973'))

        self.precise_diff_raw = np.genfromtxt(precise_diff_name, skip_header=1)

        #self.resistivity_raw = np.genfromtxt(resist_name, skip_header=1)
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
                current = int(np.round(np.abs(raw[0, col_index])))
                edict[current] = raw[1:, col_index]

        #the 3 resistivity functions
        self.resistivity_raw = {}
        self.resistivity_splines = {}
        for T in (900, 1000, 1100):
            fname = os.path.join(data_dir, "{}_Resistivity_{}K.csv".format(prefix, T))
            self.resistivity_raw[T] = np.genfromtxt(fname, skip_header=1, delimiter=',')
            self.resistivity_raw[T][:, 0] = 1 - self.resistivity_raw[T][:, 0]
            self.resistivity_splines[T] = interp.InterpolatedUnivariateSpline(self.resistivity_raw[T][:, 0],
                                                                              self.resistivity_raw[T][:, 1])

        #load up arrhenius constants
        self.diffusivity_arrhenius_constants = np.genfromtxt(diff_arrhenius_name, delimiter=',')

    def interpolated_vector(self, operative_vec, size):

        spline = interp.InterpolatedUnivariateSpline(operative_vec[:, 0], operative_vec[:, 1])
        x = np.linspace(0, 1, num=size)
        return spline(x)

    @functools.lru_cache(maxsize=512)
    def interpolated_diffusivity(self, size, T, precise=False):
        if precise:
            return self.interpolated_vector(self.precise_diff_raw, size)
        x = np.linspace(0, 1, size)
        Cmid = self.diffusivity_arrhenius_constants[:, 0]
        Dmid = self.diffusivity_arrhenius_constants[:, 1] * 1e-4 * np.exp(-self.diffusivity_arrhenius_constants[:, 2]
                                                                          * 1e3/(8.31 * T))
        #straight line on log axes...
        lnDmid = np.log(Dmid)
        cfs = np.polyfit(Cmid, lnDmid, 1)
        Dout = np.exp(cfs[0] * x + cfs[1])

        #but make it constant at the edges
        startfilter = np.where(x < Cmid[0])
        endfilter = np.where(x > Cmid[-1])
        startcopy = Dout[len(startfilter[0])]
        endcopy = Dout[-len(endfilter[0])]
        Dout[startfilter] = startcopy
        Dout[endfilter] = endcopy

        dmean = Dout.mean()
        #return np.ones(len(Dout)) * dmean
        return Dout

    @functools.lru_cache(maxsize=512)
    def interpolated_resistivity(self, size, T):
        x = np.linspace(0, 1, size)
        a2 = 0.0461 * T - 144.7
        a1 = -0.0275 * T + 161.1
        a0 = 0.00888 * T - 2.605
        #oops, forgot to switch around to Cu solute, full Cu c == 1
        #hack at it by doing 1-c
        f = lambda c: (a2 * (1-c) ** 2 + a1 * (1-c) + a0) * 10e-8
        return f(x)

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