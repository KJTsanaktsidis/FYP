from argparse import ArgumentParser
from calcsim import CalcSimWrapper
from datastore import InputDatastore
from expercomparison import ComparisonEngine
import numpy as np
import dmplots

aparser = ArgumentParser(description='Plots a comparison of experiment and simulation')
aparser.add_argument('--current', required=True, type=float,
                     help='Current density to simulate at')
aparser.add_argument('--inputdata', metavar='DIR', type=str, required=True,
                     help='Directory containing input data')
aparser.add_argument('--dataprefix', metavar='PREFIX', type=str, default='NiCu',
                     help='Prefix to data files')
aparser.add_argument('--output', metavar='FILE', type=str, required=True,
                     help='File to output to')
aparser.add_argument('--z', type=float, required=True,
                     help='Effective valence')
aparser.add_argument('--cvf', type=float, required=True,
                     help='Vacancy concentration factor')
aparser.add_argument('--direction', type=str, default='forward',
                     help='Direction of application of current')

args = aparser.parse_args()

accelcs = CalcSimWrapper()
dstore = InputDatastore(args.inputdata, args.dataprefix, 973, args.direction)
ce = ComparisonEngine(accelcs)

x = np.linspace(0, 25, num=100)
exper_data = dstore.interpolated_experiment_dict(x)[args.current]
diffusivity = dstore.interpolated_diffusivity(10001)
resistivity = dstore.interpolated_resistivity(10001)
init_cond = np.ones(100)
init_cond[50:] = 0

emigration_T = 973
dt = 0.05
ndt = int(2 * 60 * 60 / 0.05)
dx = 25e-6 / 100

r = accelcs.emigration_factor(args.z, args.current * 100 * 100, emigration_T)
simd = accelcs.calc_simulation(diffusivity, resistivity, init_cond, ndt, dt, dx, r, args.cvf)
lsq, shift = ce.calibrate(simd, exper_data)
shifted_simd = ce.shift_data(simd)
full_simd = np.column_stack((x, shifted_simd))
full_exper = np.column_stack((x, exper_data))
dmplots.plot_sim_fit(full_simd, full_exper, args.current, args.z, args.cvf, args.direction, args.output)
print('lsq = ' + str(lsq))