from calcsim import CalcSimWrapper
from datastore import InputDatastore
from expercomparison import ComparisonEngine
from paramsearch import ParamSearchEngine
from progressbar import ProgressBar, Percentage, Bar
from io import StringIO

import logging

from pylab import *


def do_param_search():
    z_list = arange(500, 1550, 50)
    diff_list = arange(1, 5.5, 0.5)
    cs = CalcSimWrapper()
    ds = InputDatastore('../InputData', 'NiCu', 973)
    ps = ParamSearchEngine(cs, ds)

    pbar = ProgressBar(widgets=['Parameter search: ', Percentage(), ' ', Bar()],
                       maxval=len(z_list) * len(diff_list))
    pbar.start()

    def progress_update(curi):
        pbar.update(curi)

    #logstream = StringIO()
    #logging.basicConfig(stream=logstream, level=logging.DEBUG)

    logging.basicConfig(level=logging.INFO)

    retim = ps.search(z_list, diff_list, 973, progress_update)

    #print(logstream.getvalue())

    extent = 1, 50, 0, 3000
    immap = imshow(retim, cmap=cm.jet, interpolation='nearest', extent=extent)
    colorbar(immap)
    ylabel('z* (arb. units)')
    xlabel('Cv/Cve (arb. units)')
    gca().set_aspect('auto')
    title('Search for minimum fit error')
    savetxt('paramsearch.csv', retim)
    savefig('paramsearch.png')


def plot_for(z, d):
    cs = CalcSimWrapper()
    ds = InputDatastore('../InputData', 'NiCu', 973)
    ce = ComparisonEngine(cs)
    D = ds.interpolated_diffusivity(10001)
    R = ds.interpolated_resistivity(10001)

    dx = 0.5 * 35e-8
    ndx = 200
    dt = 0.01
    ndt = int(2 * 60 * 60 / 0.01)

    init = ones(ndx)
    init[ndx/2:] = 0

    x = linspace(0, 25, 200)

    ddict = ds.interpolated_experiment_dict(x)

    for I in ddict.keys():
        if I == 0:
            dv = 1
        else:
            dv = d
        r = cs.emigration_factor(z, I, 973)
        mdl = cs.calc_simulation(D, R, init, ndt, dt, dx, r, dv)
        ce = ComparisonEngine(cs)
        lsq, shfit = ce.calibrate(mdl, ddict[I])
        smdl = ce.shift_data(mdl)
        plot(x, ddict[I], label=str.format('Exper. (I={} A/cm^2)', I/100/100))
        plot(x, smdl, label=str.format('Model. (I={} A/cm^2)', I/100/100))

    legend(loc=3)
    show()


def basic_test():
    cs = CalcSimWrapper()
    ds = InputDatastore('../InputData', 'NiCu', 973)
    ce = ComparisonEngine(cs)

    D = ds.interpolated_diffusivity(10001)
    R = ds.interpolated_resistivity(10001)

    dx = 0.5*35e-8
    ndx = 200
    dt = cs.optimum_dt(dx, D, 1)
    ndt = cs.num_sim_steps(dt, 2 * 60 * 60)

    init = np.ones(ndx)
    init[ndx/2:] = 0

    res = cs.calc_simulation(D, R, init, ndt, dt, dx, 0, 1)
    res = cs.calc_simulation(D, R, init, ndt, dt, dx, 0, 1)
    x = np.linspace(0, dx * ndx, num=ndx)
    plot(x, res)
    show()


def plot_for_I400(z, d):
    cs = CalcSimWrapper()
    ds = InputDatastore('../InputData', 'NiCu', 973)
    ce = ComparisonEngine(cs)
    D = ds.interpolated_diffusivity(10001)
    R = ds.interpolated_resistivity(10001)

    dx = 0.5 * 35e-8
    ndx = 200
    dt = 0.01
    ndt = int(2 * 60 * 60 / 0.01)

    init = ones(ndx)
    init[ndx/2:] = 0

    x = linspace(0, 25, 200)

    ddict = ds.interpolated_experiment_dict(x)

    for I in ddict.keys():
        if I != 400*100*100:
            continue
        if I == 0:
            dv = 1
        else:
            dv = d
        r = cs.emigration_factor(z, I, 973)
        mdl = cs.calc_simulation(D, R, init, ndt, dt, dx, r, dv)
        ce = ComparisonEngine(cs)
        lsq, shfit = ce.calibrate(mdl, ddict[I])
        smdl = ce.shift_data(mdl)
        plot(x, ddict[I], label=str.format('Exper. (I={} A/cm^2)', I/100/100))
        plot(x, smdl, label=str.format('Model. (I={} A/cm^2)', I/100/100))

    legend(loc=3)
    show()

if __name__ == '__main__':
    do_param_search()
