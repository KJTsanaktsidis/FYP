from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import RLock
from matplotlib import cm
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
import sys
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
direction = sys.argv[1]

Idensities = np.arange(0, 5100, 100)
Temps = np.arange(850, 1210, 10)
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

# with ThreadPoolExecutor(max_workers=8) as executor:
#     future_dict = dict(
#         (executor.submit(do_work, qit[0], qit[1]), qit) for qit in taskqueue
#     )
#     for res_future in as_completed(future_dict):
#         qit = future_dict[res_future]
#
#         Tindex = np.where(Temps == qit[0])[0][0]
#         Iindex = np.where(Idensities == qit[1])[0][0]
#
#         omap[Tindex, Iindex] = res_future.result()
#         i += 1
#         pbar.update(i)
#
# np.savetxt('../regionmap.csv', omap)
omap = np.genfromtxt('../pg_regionmap_{}.csv'.format(direction))

matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['legend.fontsize'] = 22
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20

fig = Figure()
ax = fig.add_subplot(111)
extent = Idensities[0], Idensities[-1], Temps[0], Temps[-1]
im = ax.imshow(omap, cmap=cm.jet, interpolation='nearest', extent=extent, origin='lower')
cb = fig.colorbar(im)
cb.set_ticks(np.arange(0, np.max(np.max(omap)), 1))
cb.set_ticklabels(['0.00', '1.00', '2.00', '3.00'])
ax.set_aspect('auto')
ax.set_xlabel(r'Current Density ($A/cm^2$)')
ax.set_ylabel('Temperature (K)')
#ax.set_title(str.format('Error relative to zero-current ({} bias)', direction))
fig.tight_layout()
canvas = FigureCanvasAgg(fig)
canvas.print_figure('../pg_regionmap_{}.png'.format(direction))
canvas.print_figure('../pg_regionmap_{}.svg'.format(direction))
