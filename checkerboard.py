from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine

__author__ = 'kj'
z = 160
cvfunc = lambda ID: 0.53e-3 * abs(ID) + 1
dstore = InputDatastore('../InputData', 'NiCu')

figure(figsize=(11, 7.7))
subplot(211)
ylabel('Cu. Composition (at. fraction)')
title('Investigation of checkerboard')
subplot(212)
xlabel('Position (micron)')
ylabel(r'$\Delta$Cu Composition (at. fraction)')

for T in (1046, 1048, 1050, 1052):
    Davg = dstore.interpolated_diffusivity(1001, T).mean()
    n_secs_sim = 1e-12 / Davg
    dt = n_secs_sim / (defaults.simulation_tsteps * 1)
    cse = CalcSimExecutor(dstore, T, dt=dt, ndt=(defaults.simulation_tsteps * 1))

    simd_lc = cse.compute(z, cvfunc(2600), 2600, 'forward')
    simd_hc = cse.compute(z, cvfunc(2700), 2700, 'forward')

    shiftengine = ComparisonEngine(cse.cs)
    shiftengine.calibrate(simd_lc[:, 1], simd_hc[:, 1])
    simd_lc[:, 1] = shiftengine.shift_data(simd_lc[:, 1])

    subplot(211)
    plot(simd_lc[:, 0], simd_lc[:, 1], label='T = {}K, I = {}A/cm^2'.format(T, 2600))
    plot(simd_hc[:, 0], simd_hc[:, 1], label='T = {}K, I = {}A/cm^2'.format(T, 2700))
    subplot(212)
    plot(simd_lc[:, 0], simd_hc[:, 1] - simd_lc[:, 1], label='T = {}K'.format(T))

subplot(211)
legend(loc='best')
subplot(212)
legend(loc='best')
savefig('../checkerboard_investigation.png')