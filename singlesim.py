from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine


T = 1200
z = 160
cvf = 0.53e-3 * 0 + 1
current = 0
T2 = 1040
direction = 'forward'
ofile = 'lel.png'

ds = InputDatastore('../InputData', 'NiCu')

#--CHANGEME
#simdata_I = cse.compute(z, cvf, current, direction)
ndt = int(600 // defaults.simulation_dt)
cse = CalcSimExecutor(ds, T, ndt=ndt)
simd = cse.compute(z, cvf, 0, 'forward')

#shiftengine = ComparisonEngine(cse.cs)
#shiftengine.calibrate(simdata_I[:, 1], simdata_noI[:, 1])
#simdata_I[:, 1] = shiftengine.shift_data(simdata_I[:, 1])

figure()

ndxmult = 1
#plot(simd[(105*ndxmult):(120*ndxmult), 0], simd[(105*ndxmult):(120*ndxmult), 1], label='T = {}K'.format(T))
plot(simd[:, 0], simd[:, 1], label='T = {}K'.format(T))

xlabel('Position (micron)')
ylabel('Cu Composition (at. fraction)')

#--CHANGEME
title('Simulation run, T = {}K'.format(T))

legend(loc='best')
savefig(ofile)