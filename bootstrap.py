from calcsim import CalcSimWrapper
from datastore import InputDatastore
from expercomparison import ComparisonEngine
from paramsearch import ParamSearchEngine
from progressbar import ProgressBar, Percentage, Bar
from io import StringIO

import logging

from pylab import *


def do_param_search():
    z_list = arange(0, 210, 10)
    diff_list = arange(0, 210, 10)
    cs = CalcSimWrapper()
    ds = InputDatastore('../InputData', 'NiCu', 973)
    ps = ParamSearchEngine(cs, ds)

    pbar = ProgressBar(widgets=['Parameter search: ', Percentage(), ' ', Bar()],
                       maxval= len(z_list) * len(diff_list))
    pbar.start()

    def progress_update(curi):
        pbar.update(curi)

    logstream = StringIO()
    logging.basicConfig(stream=logstream, level=logging.DEBUG)

    retim = ps.search(z_list, diff_list, 973, progress_update)

    print(logstream.getvalue())

    extent = 0, 200, 0, 200
    immap = imshow(retim, cmap=cm.jet, interpolation='nearest', extent=extent)
    colorbar(immap)
    xlabel('z* (arb. units)')
    ylabel('Cv/Cve (arb. units)')
    title('Search for minimum fit error')
    savetxt('paramsearch.csv', retim)
    savefig('paramsearch.png')


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