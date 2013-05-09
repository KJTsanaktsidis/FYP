"""This class is a place to stash files against; they'll be saved in a directory
"""
from os import path
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import numpy as np


class OutputDatastore():

    def __init__(self, directory, prefix=''):
        """Constructor should specifiy what directory output data will be saved in
        """

        self.directory = directory
        self.prefix = prefix

    def add_plot(self, fig, fname):
        assert isinstance(fig, Figure)

        fullfname = path.join(self.directory, self.prefix + fname)
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(fullfname)

    def add_matrix(self, mat, fname):
        fullfname = path.join(self.directory, self.prefix + fname)
        np.savetxt(fullfname, mat, delimiter=',')

    def retrieve_matrix(self, fname):
        fullfname = path.join(self.directory, self.prefix + fname)
        return np.genfromtxt(fullfname, delimiter=',')

    def complete(self):
        pass
