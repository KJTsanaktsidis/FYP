from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import RLock
from matplotlib import cm
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from calcsim import SimulationUnstableError
from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
import defaults
from diffsimtask import DiffSimTask
from expercomparison import ComparisonEngine
from outstore import OutputDatastore
import numpy as np
import itertools
import progressbar

dstore = InputDatastore('../InputData', 'NiCu')
dstorelock = RLock()

z = 160
cvfunc = lambda ID: 0.53e-3 * abs(ID) + 1
direction = 'forward'

Idensities = np.arange(0, 5020, 100)
Temps = np.arange(1030, 1070, 2)
omap = np.zeros((len(Temps), len(Idensities)))

taskqueue = itertools.product(Temps, Idensities)

print('no. iters: ' + str(len(Idensities) * len(Temps)))
pbar = progressbar.ProgressBar(maxval=len(Idensities) * len(Temps), widgets=[progressbar.Percentage(), progressbar.Bar()])
pbar.start()
i = 0


def do_work(T, I):
    with dstorelock:
        Davg = dstore.interpolated_diffusivity(1001, T).mean()
        n_secs_sim = 1e-12 / Davg
        dt = n_secs_sim / (defaults.simulation_tsteps * 1)
        cse = CalcSimExecutor(dstore, T, dt=dt, ndt=(defaults.simulation_tsteps * 1))

    try:
        current = cse.compute(z, cvfunc(I), I, direction)[:, 1]
        nocurrent = cse.compute(0, 1, 0, direction)[:, 1]
    except SimulationUnstableError:
        print('Fucked up on T = {}, I = {}'.format(T, I))
        print('Max D = {}'.format(cse.Dvector.max()))
        raise
    shiftengine = ComparisonEngine(cse.cs)
    lsq, _ = shiftengine.calibrate(current, nocurrent)
    return lsq

with ThreadPoolExecutor(max_workers=8) as executor:
    future_dict = dict(
        (executor.submit(do_work, qit[0], qit[1]), qit) for qit in taskqueue
    )
    for res_future in as_completed(future_dict):
        qit = future_dict[res_future]

        Tindex = np.where(Temps == qit[0])[0][0]
        Iindex = np.where(Idensities == qit[1])[0][0]

        omap[Tindex, Iindex] = res_future.result()
        i += 1
        pbar.update(i)

np.savetxt('../regionmap.csv', omap)


fig = Figure()
ax = fig.add_subplot(111)
extent = Idensities[0], Idensities[-1], Temps[0], Temps[-1]
im = ax.imshow(omap, cmap=cm.jet, interpolation='nearest', extent=extent, origin='lower')
fig.colorbar(im)
ax.set_aspect('auto')
ax.set_xlabel('Current Density (A/cm^2)')
ax.set_ylabel('Temperature (K)')
ax.set_title(str.format('Error relative to zero-current ({} bias)', direction))
canvas = FigureCanvasAgg(fig)
canvas.print_figure('../regionmap.png')
