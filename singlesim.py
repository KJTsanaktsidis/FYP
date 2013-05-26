from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine


T = 1052
z = 160
cvf = 0.53e-3 * 0 + 1
current = 0
T2 = 1040
direction = 'forward'
ofile = 'lel.png'

ds = InputDatastore('../InputData', 'NiCu')

#--CHANGEME
#simdata_I = cse.compute(z, cvf, current, direction)
ndt = int(10 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simdata_10s = cse.compute(z, cvf, 0, 'forward')

ndt = int(60 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simdata_1m = cse.compute(z, cvf, 0, 'forward')

ndt = int(600 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simdata_10m = cse.compute(z, cvf, 0, 'forward')

ndt = int(30 * 60 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simdata_30m = cse.compute(z, cvf, 0, 'forward')

ndt = int(60 * 60 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simdata_1h = cse.compute(z, cvf, 0, 'forward')

#shiftengine = ComparisonEngine(cse.cs)
#shiftengine.calibrate(simdata_I[:, 1], simdata_noI[:, 1])
#simdata_I[:, 1] = shiftengine.shift_data(simdata_I[:, 1])

figure()

#--CHANGEME
plot(simdata_10s[:, 0], simdata_10s[:, 1], label='t = 10s')
plot(simdata_1m[:, 0], simdata_1m[:, 1], label='t = 1m')
plot(simdata_10m[:, 0], simdata_10m[:, 1], label='t = 10m')
plot(simdata_30m[:, 0], simdata_30m[:, 1], label='t = 30m')
plot(simdata_1h[:, 0], simdata_1h[:, 1], label='t = 1h')

xlabel('Position (micron)')
ylabel('Cu Composition (at. fraction)')

#--CHANGEME
title('Simulation run, T = {}K'.format(T))

legend(loc='best')
savefig(ofile)