from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine


T = 700+273
z = 160
I = 1000
cvf = 0.53e-3 * I + 1
#cvf = 1
ofile = 'lel.png'
ds = InputDatastore('../InputData', 'NiCu')

cse = CalcSimExecutor(ds, T)
simd_forward = cse.compute(z, cvf, I, 'forward')
simd_reverse = cse.compute(z, cvf, I, 'reverse')
x = simd_forward[:, 0]
fedict = ds.interpolated_experiment_fw_dict(x)
redict = ds.interpolated_experiment_rv_dict(x)
expr_forward = fedict[I]
expr_reverse = redict[I]

shifter = ComparisonEngine(cse.cs)
shifter.calibrate(simd_forward[:, 1], expr_forward)
simd_forward[:, 1] = shifter.shift_data(simd_forward[:, 1])
shifter.calibrate(simd_reverse[:, 1], expr_reverse)
simd_reverse[:, 1] = shifter.shift_data(simd_reverse[:, 1])

figure()

plot(x, expr_forward, 'b', label=r'Ni->Cu $e^-$ flow (Zhao et. al.)')
plot(x, simd_forward[:, 1], 'r--', label=r'Ni->Cu $e^-$ flow (simulation)')
plot(x, expr_reverse, 'g', label=r'Cu->Ni $e^-$ flow (Zhao et. al.)')
plot(x, simd_reverse[:, 1], color='orange', linestyle='--', label='Cu->Ni $e^-$ flow (simulation)')

xlabel('Position (micron)')
ylabel('Cu Composition (at. fraction)')
xlim([5, 20])

legend(loc='best', prop={'size': 10})
savefig(ofile)