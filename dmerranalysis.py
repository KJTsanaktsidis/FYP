from argparse import ArgumentParser
from os import path
import os
import re
from matplotlib.backends.backend_agg import FigureCanvasAgg
from calcsim import CalcSimWrapper
from datastore import InputDatastore
from expercomparison import ComparisonEngine
import numpy as np
import dmplots


aparser = ArgumentParser(description='Searches for the optimum vacancy concentration multiplier')
aparser.add_argument('--inputdata', metavar='DIR', type=str, required=True,
                     help='Directory containing input data')
aparser.add_argument('--outputdir', metavar='DIR', type=str, default='.',
                     help='Directory to stash the output results (plots, tables etc.)')
aparser.add_argument('--dataprefix', metavar='PREFIX', type=str, default='NiCu',
                     help='Prefix to data files')
aparser.add_argument('--temperature', type=float, default=973,
                     help='Temperature in K')

args = aparser.parse_args()

accelcs = CalcSimWrapper()

dstore = InputDatastore(args.inputdata, args.dataprefix)
x = np.linspace(0, 25, num=100)
init_cond = np.ones(100)
init_cond[50:] = 0
dt = 0.05
ndt = int(2 * 60 * 60 / 0.05)
dx = 25e-6 / 100


def quicksim(z, cvf, I, exper, direction):
    ce = ComparisonEngine(accelcs)
    if direction == 'reverse':
        I = -np.abs(I)
    else:
        I = np.abs(I)
    diffusivity = dstore.interpolated_diffusivity(10001, args.temperature, precise=True)
    resistivity = dstore.interpolated_resistivity(10001, args.temperature)
    r = accelcs.emigration_factor(z, I * 100 * 100, args.temperature)
    simd = accelcs.calc_simulation(diffusivity, resistivity, init_cond, ndt, dt, dx, r, cvf)
    lsq, shift = ce.calibrate(simd, exper)
    shifted_simd = ce.shift_data(simd)
    return shifted_simd, lsq


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


#http://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

#OK, get the data out of outputdir
fresults = {}
rresults = {}
results = {'forward': fresults, 'reverse': rresults}

zrange = np.genfromtxt(path.join(args.outputdir, 'zrange.csv'), delimiter=',')
cvfrange = np.genfromtxt(path.join(args.outputdir, 'cvfrange.csv'), delimiter=',')

files = [path.join(args.outputdir, f) for f in os.listdir(args.outputdir)
         if path.isfile(path.join(args.outputdir, f))]
forwardfiles = [f for f in files if f.endswith('forward.csv')]
reversefiles = [f for f in files if f.endswith('reverse.csv')]
rx = r'.*Searchmap_I([\-0-9\.]+)_.*.csv'
for f in forwardfiles:
    print(str.format('Processing {}', f))
    I = int(float(re.match(rx, f).groups()[0]))
    fresults[I] = np.genfromtxt(f, skip_header=0, delimiter=',')
for f in reversefiles:
    print(str.format('Processing {}', f))
    I = np.abs(int(float(re.match(rx, f).groups()[0])))
    rresults[I] = np.genfromtxt(f, skip_header=0, delimiter=',')

#now, find average z* for forward case
fzaverage = 0
ict = 0
zvalues = []
cresults = {}

#multiply forward/reverse maps to get current-min maps
for I in fresults.keys():
    combined_map = fresults[I] * rresults[I]
    cresults[I] = combined_map
    if I == 0:
        continue
    min_indicies = np.unravel_index(combined_map.argmin(), combined_map.shape)
    min_z = zrange[min_indicies[0]]
    print(str.format('Combined z* for I = {}: {}', I, min_z))
    zvalues.append(min_z)

fzaverage = np.array(zvalues).mean()
print(str.format('Average z*: {}', fzaverage))
zaverage_rounded, zaverage_index = find_nearest(zrange, fzaverage)
print(str.format('Rounding to {}', zaverage_rounded))

#find (min, best, max) cvf at zaverage_rounded for each I
cbounds = {}
for I in cresults.keys():
    bounds_column = cresults[I][zaverage_index, :]
    bounds_truth = bounds_column < 0.005
    contig_idxs = contiguous_regions(bounds_truth)
    contig_sizes = contig_idxs[:, 1] - contig_idxs[:, 0]
    contig_largest_idx = contig_sizes.argmax()
    contig_bounds = bounds_column[contig_idxs[contig_largest_idx, 0]: contig_idxs[contig_largest_idx, 1]]

    bounds_min = cvfrange[contig_idxs[contig_largest_idx][0]]
    bounds_max = cvfrange[contig_idxs[contig_largest_idx][1] - 1]
    bounds_best = cvfrange[bounds_column.argmin()]

    assert(bounds_column.argmin() >= contig_idxs[contig_largest_idx][0])
    assert(bounds_column.argmin() <= contig_idxs[contig_largest_idx][1])

    print(str.format('I = {}, min, max, best = {}, {}, {}', I, bounds_min, bounds_max, bounds_best))
    cbounds[I] = bounds_best - bounds_min, bounds_best, bounds_max - bounds_best

#now plot with error bars
I_plot = np.array([b for b in cbounds.keys()])
cvf_plot_list = []
for I in I_plot:
    cvf_plot_list.append(cbounds[I][1])
cvf_plot = np.array(cvf_plot_list)
ebars_list_lower = []
ebars_list_upper = []
for I in I_plot:
    ebars_list_lower.append(cbounds[I][0])
    ebars_list_upper.append(cbounds[I][2])
ebars = np.array([ebars_list_lower, ebars_list_upper])
plot_data = np.column_stack((I_plot, cvf_plot))
outfile = path.join(args.outputdir, 'cvfplot_combined.png')
print(I_plot)
print(ebars)
dmplots.plot_cvf_function_ebars(plot_data, ebars, 'combined', outfile)

#and politely output simulations
for direction in ('forward', 'reverse'):
    for I in cresults.keys():
        for idx, edge in ((0, 'lower'), (1, 'best'), (2, 'upper')):
            outfile = path.join(args.outputdir, str.format('Comparison_I{}_{}_{}bound.png', I, direction, edge))
            cvf = cbounds[I][1]
            if idx == 0:
                cvf -= cbounds[I][idx]
            elif idx == 2:
                cvf += cbounds[I][idx]
            print('using cvf = ' + str(cvf))
            edict = dstore.edict_for_direction(direction)
            exper = dstore.interpolated_experiment_dict(x, edict)[I]
            simd, lsq = quicksim(zaverage_rounded, cvf, I, exper, direction)
            simd = np.column_stack((x, simd))
            exper = np.column_stack((x, exper))
            f = dmplots.plot_sim_fit(simd, exper, I, zaverage_rounded, cvf, direction)
            c = FigureCanvasAgg(f)
            c.print_figure(outfile)


