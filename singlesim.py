from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from pylab import *
import defaults
from expercomparison import ComparisonEngine


T = 700+273
z = 700
#I = 1000
cvf = 0.53e-3 * 0 + 1
#cvf = 1
ofile = 'lel'
ds = InputDatastore('../InputData', 'NiCu')
direction = 'forward'

#matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['axes.labelsize'] = 24
matplotlib.rcParams['legend.fontsize'] = 22
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20

def get_edict(x):
    return ds.interpolated_experiment_fw_dict(x) if direction == 'forward' else ds.interpolated_experiment_rv_dict(x)

figure()
cse = CalcSimExecutor(ds, T)
colours = ['#0000FF', '#FF0000', '#003300', '#DD8500']
colours.reverse()
for I, sh in [[0, 0.7], [400, 0.5], [800, 0.06], [1000, 0]]:
    simd = cse.compute(z, cvf, I, direction)
    x = simd[:, 0]
    edict = get_edict(x)
    expr = edict[I]
    shifter = ComparisonEngine(cse.cs)
    shifter.calibrate(simd[:, 1], expr)
    #simd[:, 1] = shifter.shift_data(simd[:, 1])

    x += sh
    cl = colours.pop()
    #plot(x, expr, linestyle='-', label='I = {}A/cm^2'.format(I))
    plot(x, simd[:, 1], linestyle='-', color=cl, linewidth=2, label=r'$\mathrm{'+str(I)+'A/cm}^2$')

xlim([5, 20])
xticks(arange(5, 21, 5))
xlabel('Position (micron)')
ylabel('Cu Composition (at. frac.)')
lg = legend(loc='center left')
lg.get_frame().set_alpha(0)
lg.get_frame().set_edgecolor('white')
savefig(ofile + '.svg')
savefig(ofile + '.png')