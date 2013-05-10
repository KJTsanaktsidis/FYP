from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import RLock
from matplotlib import cm
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from calcsim import SimulationUnstableError
from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from diffsimtask import DiffSimTask
from outstore import OutputDatastore
import numpy as np
import itertools
import progressbar

dstore = InputDatastore('../InputData', 'NiCu')
dstorelock = RLock()

z = 1875
cvfunc = lambda ID: 0.0014 * ID

Idensities = np.arange(0, 10200, 200)
Temps = np.arange(850, 1425, 25)
omap = np.zeros((len(Temps), len(Idensities)))

taskqueue = itertools.product(Temps, Idensities)

print('no. iters: ' + str(len(Idensities) * len(Temps)))
pbar = progressbar.ProgressBar(maxval=len(Idensities) * len(Temps), widgets=[progressbar.Percentage(), progressbar.Bar()])
pbar.start()
i = 0


def do_work(T, I):
    with dstorelock:
        cse = CalcSimExecutor(dstore, T)
    try:
        current = cse.compute(z, cvfunc(I), I, 'forward')[:, 1]
        nocurrent = cse.compute(0, 1, 0, 'forward')[:, 1]
    except SimulationUnstableError:
        print('Fucked up on T = {}, I = {}'.format(T, I))
        print('Max D = {}'.format(cse.Dvector.max()))
        raise
    lsq_arr = (current - nocurrent) ** 2
    return np.sum(lsq_arr)

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
im = ax.imshow(omap, cmap=cm.jet, interpolation='nearest', extent=extent)
fig.colorbar(im)
ax.set_aspect('auto')
ax.set_xlabel('Current Density (A/cm^2)')
ax.set_ylabel('Temperature (K)')
ax.set_title(str.format('Error relative to zero-current'))
canvas = FigureCanvasAgg(fig)
canvas.print_figure('../regionmap.png')
