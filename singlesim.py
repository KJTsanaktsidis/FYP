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
ndt = int(600 // defaults.simulation_dt)
Tks = [1049, 1049.5, 1050, 1050.5, 1051, 1051.5, 1052]
simds = []

for Tk in Tks:
    cse = CalcSimExecutor(ds, Tk, ndt=ndt)
    sd = cse.compute(z, cvf, 0, 'forward')
    simds.append(sd)

#shiftengine = ComparisonEngine(cse.cs)
#shiftengine.calibrate(simdata_I[:, 1], simdata_noI[:, 1])
#simdata_I[:, 1] = shiftengine.shift_data(simdata_I[:, 1])

figure()

ndxmult = 1
for Tk, simd in zip(Tks, simds):
    plot(simd[(105*ndxmult):(120*ndxmult), 0], simd[(105*ndxmult):(120*ndxmult), 1], label='T = {}K'.format(Tk))

xlabel('Position (micron)')
ylabel('Cu Composition (at. fraction)')

#--CHANGEME
title('Simulation run, T = {}K'.format(T))

legend(loc='best')
savefig(ofile)