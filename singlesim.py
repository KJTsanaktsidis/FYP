from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
from expercomparison import ComparisonEngine


T = 1052
z = 160
cvf = 0.53e-3 * 3500 + 1
current = 3500
T2 = 1040
direction = 'forward'
ofile = 'lel.png'

ds = InputDatastore('../InputData', 'NiCu')

cse = CalcSimExecutor(ds, T)
simdata_I = cse.compute(z, cvf, current, direction)
simdata_noI = cse.compute(z, 1, 0, direction)

shiftengine = ComparisonEngine(cse.cs)
shiftengine.calibrate(simdata_I[:, 1], simdata_noI[:, 1])
simdata_I[:, 1] = shiftengine.shift_data(simdata_I[:, 1])

figure()
plot(simdata_I[:, 0], simdata_I[:, 1], label='I = 3500A/cm^2')
plot(simdata_noI[:, 0], simdata_noI[:, 1], label='I = 0A/cm^2')
xlabel('Position (micron)')
ylabel('Cu Composition (at. fraction)')
title('Simulation run, Extrapolated diffusivity, T = {}K'.format(T))
legend(loc='best')
savefig(ofile)