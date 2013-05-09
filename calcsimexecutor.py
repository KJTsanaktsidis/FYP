"""This file contains a class that simplifies setting up a simulation
"""
from calcsim import CalcSimWrapper
from datastore import InputDatastore
import defaults
import numpy as np


class CalcSimExecutor():

    def __init__(self, dstore, T):
        assert isinstance(dstore, InputDatastore)

        self.T = T

        #first, get defaults
        self.dx = defaults.simulation_dx
        self.dt = defaults.simulation_dt
        self.ndx = defaults.simulation_xsteps
        self.ndt = defaults.simulation_tsteps

        #now we can set up the initial conditions
        #x in micron
        self.x = np.arange(0, self.ndx * self.dx, self.dx) * 1e6
        self.init_cond = np.ones(self.ndx)
        self.init_cond[(self.ndx // 2):] = 0

        #D/R vectors
        self.Dvector = dstore.interpolated_diffusivity(10001, T)
        self.Rvector = dstore.interpolated_resistivity(10001, T)

        self.cs = CalcSimWrapper()
        #now we're ready to fire on demand

    def compute(self, z, cvf, I, direction):
        """
        Actually runs the simulation that's been set up, with the supplied
        parameters. I is supplied in A/cm^2
        Returns sim results as a 2 column array of x, y
        """

        if direction == 'forward':
            pass
        elif direction == 'reverse':
            I = -I
        else:
            raise ValueError('Unknown direction ' + str(direction))

        r = self.cs.emigration_factor(z, I * 100 * 100, self.T)
        outy = self.cs.calc_simulation(self.Dvector, self.Rvector, self.init_cond,
                                       self.ndt, self.dt, self.dx, r, cvf)
        return np.column_stack((self.x, outy))