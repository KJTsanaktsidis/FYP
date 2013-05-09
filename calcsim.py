"""Python ctypes wrapper for calcsim.c
"""
import logging
import numpy as np
import ctypes


class CalcSimWrapper:

    def __init__(self):
        try:
            self.libcalcsim = ctypes.cdll.LoadLibrary('libcalcsim.so')
        except OSError:
            self.libcalcsim = ctypes.cdll.LoadLibrary('./libcalcsim.so')

        self.libcalcsim.calc_simulation.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.c_int32, ctypes.POINTER(ctypes.c_double),
            ctypes.c_int32, ctypes.c_int32, ctypes.c_double, ctypes.c_double,
            ctypes.c_double, ctypes.POINTER(ctypes.c_double)
        ]
        self.libcalcsim.fast_pad_shift.argtypes = [
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int32
        ]

    def emigration_factor(self, z, I_density, T):
        """
        Computes the electromigration factor r

        :param z: Effective valence
        :param I_density: Current density in A/m^2
        :param T: Temperature in K
        """

        return (1.602e-19 * I_density * z) / (1.381e-23 * T)

    def num_sim_steps(self, dt, tmax):
        """
        Computes the number of timesteps needed to get the simulation past tmax

        :param dt: The size of a timestep
        :param tmax: The desitred simulation time
        """

        return int(np.ceil(tmax / dt))

    def optimum_dt(self, dx, D_vector, cv_factor):
        """
        Computes a dt that should make the simulation stable

        :param dx: Size of one timestep, in metres
        :param D_vector: Vector of diffusivity values in SI units
        :param cv_factor: multiplicative factor for diffusivity
        """

        max_diff = np.max(D_vector) * cv_factor
        return 0.1 * dx ** 2 / max_diff

    def calc_simulation(self,
                        D_vector, R_vector, init_cond,
                        ndt, dt, dx, r, cv_factor):
        """
        Wrapper for the calc_simulation function in calcsim.c
        This function will compute the forward-difference diffusion equation involving electromigration and accelerated
        diffusion (i.e. Cv/Cve > 1).

        :param D_vector: identically sized vectors of diffusivity and resistivity as a function of concentration from
            0 to 1. DO NOT MODIFY THESE FROM ANOTHER THREAD
        :param R_vector: identically sized vectors of diffusivity and resistivity as a function of concentration from
            0 to 1. DO NOT MODIFY THESE FROM ANOTHER THREAD
        :param init_cond: initial concentration profile. DO NOT MODIFY FROM ANOTHER THREAD
        :param ndt: number of timesteps to run
        :param dt: size of one timestep in seconds
        :param dx: size of one spacestep in metres
        :param r: multiplicative factor involving z*, T and Idensity, as well as diverse other constants
        :param cv_factor: multiplicative factor for diffusivity
        """

        if len(D_vector) != len(R_vector):
            raise ValueError("D and R must be the same length")
        if len(D_vector) < 1:
            raise ValueError("There needs to be at least one element in D and R")

        #easy scalars
        nIV = len(D_vector)
        ndx = len(init_cond)

        #multiply diffusivity by the relevant factor
        D_revised = D_vector * cv_factor

        #make an output container
        out_array = np.zeros(ndx)

        #need things to be contiguous
        D_contig = np.ascontiguousarray(D_revised, dtype=np.float64)
        R_contig = np.ascontiguousarray(R_vector, dtype=np.float64)
        init_cond_contig = np.ascontiguousarray(init_cond, dtype=np.float64)
        out_contig = np.ascontiguousarray(out_array, dtype=np.float64)

        #vars
        D_ptr = D_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        R_ptr = R_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        init_cond_ptr = init_cond_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        out_ptr = out_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        #and shell out
        logging.debug(str.format('About to move to native code (r = {})', r))
        res = self.libcalcsim.calc_simulation(D_ptr, R_ptr, nIV, init_cond_ptr,
                                              ndt, ndx, dt, dx, r, out_ptr)
        if res == 1:
            raise SimulationUnstableError("Simulation was unstable, aborting...")
        return out_contig

    def fast_pad_shift(self, y1, y2):
        """
        Wrapper around fast_pad_shift() in fastshift.c

        Works out the optimum number of indicies by which y2 needs to be shifted to have the minimum
        least squares error between the two

        :param y1: The array that will not be shifted
        :param y2: The array that will be shifted
        """

        if len(y1) != len(y2):
            raise ValueError("Input sizes must be the same")

        y1_contig = np.ascontiguousarray(y1, dtype=np.float64)
        y2_contig = np.ascontiguousarray(y2, dtype=np.float64)
        y1_ptr = y1_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        y2_ptr = y2_contig.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        return self.libcalcsim.fast_pad_shift(y1_ptr, y2_ptr, len(y1))


class SimulationUnstableError(Exception):
    pass
