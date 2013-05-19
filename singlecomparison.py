"""Contains the code for computing and plotting a single comparison between
simdata and experimental data
"""
from argparse import ArgumentParser
from calcsim import CalcSimWrapper, SimulationUnstableError
from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from diffsimtask import DiffSimTask
from expercomparison import ComparisonEngine
from outstore import OutputDatastore
import numpy as np
import dmplots


class SingleComparisonTask(DiffSimTask):

    def __init__(self, dstore, ostore, args):
        assert isinstance(dstore, InputDatastore)
        assert isinstance(ostore, OutputDatastore)

        self.dstore = dstore
        self.ostore = ostore
        self.args = args

    @staticmethod
    def add_arg_options(subparsercol):

        argparser = subparsercol.add_parser('singlecomp')

        argparser.add_argument('--z', type=float, required=True,
                               help='Effective valence of an atom in the simulation')
        argparser.add_argument('--cvf', type=float, required=True,
                               help='Vacancy concentration boost factor')
        argparser.add_argument('--current', type=int, required=True,
                               help='Current that should be used for simulation, and to select experimental '
                                    'data to compare with')
        argparser.add_argument('--temperature', type=int, default=973,
                               help='Temperature to conduct simulation at, in K')
        argparser.add_argument('--direction', type=str, default='forward',
                               help='Direction of the current being applied')

    def do_stage1(self):
        """Our stage 1 is to do the requested simulation and spit out a csv file with (x, sim, exper)
        """
        simexec = CalcSimExecutor(self.dstore, self.args.temperature)
        try:
            simdata = simexec.compute(self.args.z, self.args.cvf, self.args.current, self.args.direction)
        except SimulationUnstableError:
            print('Bang! Dmax = {}, Dmin = {}'.format(simexec.Dvector.max(), simexec.Dvector.min()))
            raise
        experdata = self.dstore.interpolated_experiment(self.args.current, simdata[:, 0], self.args.direction)

        ce = ComparisonEngine(CalcSimWrapper())
        ce.calibrate(simdata[:, 1], experdata)
        simdata[:, 1] = ce.shift_data(simdata[:, 1])

        self.fulldata = np.column_stack((simdata, experdata))

        self.ostore.add_matrix(self.fulldata, '.csv')

    def do_stage2(self):
        """
        Our stage 2 is to do the required plotting
        """

        if self.args.resume:
            self.fulldata = self.ostore.retrieve_matrix('.csv')

        simdata = self.fulldata[:, 0:2]
        experdata = np.column_stack((self.fulldata[:, 0], self.fulldata[:, 2]))
        fig = dmplots.plot_sim_fit(simdata, experdata, self.args.current, self.args.z,
                                   self.args.cvf, self.args.direction)
        self.ostore.add_plot(fig, '.png')