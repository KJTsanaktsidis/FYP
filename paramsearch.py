"""Module for performing parameter search over z* and Cv/Cve
"""
import logging

import numpy as np
import itertools
from concurrent import futures

from concurrent.futures import ThreadPoolExecutor
from threading import RLock

from calcsim import CalcSimWrapper, SimulationUnstableError
from datastore import InputDatastore
from expercomparison import ComparisonEngine


class ParamSearchEngine():

    def __init__(self, calcsim_wrapper, input_datastore):
        """
        Constructor for injecting dependancies into ParamSearchEngine

        :param calcsim_wrapper: instance of CalcSimWrapper
        :param input_datastore: instance of InputDatastore
        """
        assert isinstance(calcsim_wrapper, CalcSimWrapper)
        assert isinstance(input_datastore, InputDatastore)

        self.calcsim_wrapper = calcsim_wrapper
        self.input_datastore = input_datastore

    def search(self, z_list, dmult_list, emigration_T, progress_cb):
        """
        This method will compute calcsim_wrapper.calc_simulation() for every value of z* and Cv in
        z_range and dmult_range and for every current in input_datastore.
        The error of the simulation against experimental data is compared for each current for which there is data
        The cell of the array representing z, Cv is set to the sum of all of the least squared errors

        """

        #set up the work

        self.x = np.linspace(0, 25, num=200)
        self.exper_data = self.input_datastore.interpolated_experiment_dict(self.x)
        self.diffusivity = self.input_datastore.interpolated_diffusivity(10001)
        self.resistivity = self.input_datastore.interpolated_resistivity(10001)

        self.init_cond = np.ones(200)
        self.init_cond[100:] = 0

        self.emigration_T = emigration_T
        self.dt = 0.01
        self.ndt = int(2 * 60 * 60 / 0.01)
        self.dx = 25e-6 / 200

        #now create a queue
        work_queue = itertools.product(z_list, dmult_list)

        #for q in work_queue:
        #    print('(' + str(q[0]) + ', ' + str(q[1]) + ')\n')

        #create storage
        ret_storage = np.zeros((len(z_list), len(dmult_list)))

        unstable_count = 0

        logging.info('Testing logging')

        with ThreadPoolExecutor(max_workers=4) as executor:
            future_dict = dict(
                ((executor.submit(self.do_work, qit[0], qit[1]), qit)
                 for qit in work_queue)
            )
            pcount = 0
            for res_future in futures.as_completed(future_dict):
                qit = future_dict[res_future]

                z_index = np.nonzero(z_list == qit[0])[0][0]
                d_index = np.nonzero(dmult_list == qit[1])[0][0]
                try:
                    simres = res_future.result()
                    logging.info(str.format(
                        'Writing ({}, {}) (val = {}) into position ({}, {})',
                        qit[0], qit[1], simres, z_index, d_index
                    ))
                    ret_storage[z_index, d_index] = simres
                except SimulationUnstableError:
                    logging.warning(str.format(
                        'Simulation unstable for ({}, {}) for position ({}, {})',
                        qit[0], qit[1], z_index, d_index
                    ))
                    unstable_count += 1

                pcount += 1
                progress_cb(pcount)

        logging.info('Simulation unstable count: ' + str(unstable_count))
        return ret_storage

    def do_work(self, z, dmult):
        """
        Actually does one item of work from the search() method

        :param z: The effective valence
        :param dmult: The premultiplication factor for D
        :return: A measure of how good the model fits the data for these parameters
        """

        logging.info(str.format('Executing workload ({}, {}))', z, dmult))
        lsq_results = list()

        for I in self.exper_data.keys():
            if I == 0:
                dv = 1
            else:
                dv = dmult
            r = self.calcsim_wrapper.emigration_factor(z, I, self.emigration_T)
            simresults = self.calcsim_wrapper.calc_simulation(self.diffusivity, self.resistivity,
                                                              self.init_cond, self.ndt, self.dt, self.dx,
                                                              r, dv)
            ce = ComparisonEngine(self.calcsim_wrapper)
            lsq, shift = ce.calibrate(simresults, self.exper_data[I])
            lsq_results.append(lsq)

        return sum(lsq_results)