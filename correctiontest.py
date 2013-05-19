from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine

dstore = InputDatastore('../InputData', 'NiCu')

z = 1875
I = 3000
cvfunc = lambda ID: 0.0014 * ID

ioff()
figure()
for T in (1000, 1100, 1200):
    Davg = dstore.interpolated_diffusivity(1001, T).mean()
    n_secs_sim = 1e-12 / Davg
    ndt = int(n_secs_sim / defaults.simulation_dt)
    cse = CalcSimExecutor(dstore, T)

    current_whole = cse.compute(z, cvfunc(I), I, 'forward')
    nocurrent = cse.compute(0, 1, 0, 'forward')[:, 1]

    x = current_whole[:, 0]
    current = current_whole[:, 1]

    shiftengine = ComparisonEngine(cse.cs)
    lsq, _ = shiftengine.calibrate(current, nocurrent)
    s_current = shiftengine.shift_data(current)

    plot(x, s_current, label='{}K, current'.format(T))
    plot(x, nocurrent, label='{}K, no current'.format(T))

xlabel('x coordinate (micron)')
ylabel('Cu concentration (at. fraction)')
title('Comparison of different Temps with correction')
legend(loc='lower left')

show()