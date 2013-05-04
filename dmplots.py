"""Plotting utilities for the vacancy search
"""
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import cm


def plot_search_map(datamap, zlim, cvflim, I, direction, fname):
    fig = Figure()
    ax = fig.add_subplot(111)
    extent = cvflim[0], cvflim[1], zlim[0], zlim[1]
    im = ax.imshow(datamap, cmap=cm.jet, interpolation='nearest', extent=extent)
    fig.colorbar(im)
    ax.set_aspect('auto')
    ax.set_xlabel('Vacancy concentration multiplier')
    ax.set_ylabel('z* (Effective valence)')
    ax.set_title(str.format('Least Squares Error for I = {}, {} bias', I, direction))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fname)


def plot_cvf_function(data, direction, fname):
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:, 0], data[:, 1], 'bx')
    ax.set_xlim((0, data[:, 0].max() * 1.2))
    ax.set_ylim((0, data[:, 1].max() * 1.2))
    ax.set_xlabel('Current Density (A/cm^2)')
    ax.set_ylabel('Vacancy concentration multiplier')
    ax.set_title(str.format('Vacancy concentration for {} bias', direction))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fname)


def plot_z_function(data, direction, fname):
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.plot(data[:, 0], data[:, 1], 'bx')
    ax.set_xlim((0, data[:, 0].max() * 1.2))
    ax.set_ylim((0, data[:, 1].max() * 1.2))
    ax.set_xlabel('Current Density (A/cm^2)')
    ax.set_ylabel('Effective valence (z*)')
    ax.set_title(str.format('Effective valence for {} bias', direction))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fname)


def plot_sim_fit(simdata, experdata, I, z, cvf, direction, fname):
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.plot(simdata[:, 0], simdata[:, 1], 'b', label='Simulation')
    ax.plot(experdata[:, 0], experdata[:, 1], 'r', label='Experimental')
    ax.set_xlim((0, 25))
    ax.set_ylim((0, 1))
    ax.set_xlabel('Position (micron)')
    ax.set_ylabel('Concentration fraction of Cu')
    ax.legend()
    ax.set_title(str.format('Simulation for I = {} A/cm^2, {} bias, z* = {}, Cvf = {}', I, direction, z, cvf))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fname)


def plot_cvf_function_ebars(data, ebars, direction, fname):
    fig = Figure()
    ax = fig.add_subplot(111)
    ax.errorbar(data[:, 0], data[:, 1], xerr=None, yerr=ebars, elinewidth=1, fmt='bx')
    ax.set_xlim((0, data[:, 0].max() * 1.2))
    ax.set_ylim((0, (data[:, 1].max() + ebars.max().max()) * 1.2))
    ax.set_xlabel('Current Density (A/cm^2)')
    ax.set_ylabel('Vacancy concentration multiplier')
    ax.set_title(str.format('Vacancy concentration for {} bias', direction))
    canvas = FigureCanvas(fig)
    canvas.print_figure(fname)