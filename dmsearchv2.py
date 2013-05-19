from matplotlib.backends.backend_agg import FigureCanvasAgg
from progressbar import ProgressBar, Percentage, Bar
from calcsim import CalcSimWrapper
from calcsimexecutor import CalcSimExecutor
from datastore import InputDatastore
from dmplots import plot_search_map
import dmplots
from expercomparison import ComparisonEngine
from paramsearch import ParamSearchEngine
from os import path

from pylab import *

T = 973
zrange = arange(0, 2010, 10)
cvfrange = arange(1, 3.1, 0.1)

cs = CalcSimWrapper()
ds = InputDatastore('../InputData', 'NiCu')
psearch = ParamSearchEngine(cs, ds)

pbar = ProgressBar(widgets=['Parameter search: ', Percentage(), ' ', Bar()],
                   maxval=len(zrange) * len(cvfrange) * 8)
pbar.start()


def update_pbar(val):
    pbar.update(pbar.currval + 1)


def get_fname(name):
    return path.join('../dmsearch_v2_output', name)


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
    idx.shape = (-1, 2)
    return idx

results = {}

for direct in ('forward', 'reverse'):
    edict = ds.edict_for_direction(direct)
    results[direct] = {}
    for I in edict.keys():
        #sres = psearch.search(zrange, cvfrange, I, T, update_pbar, direct)
        #savetxt(get_fname('{}_I{}.csv'.format(direct, I)), sres, delimiter=',')
        sres = genfromtxt(get_fname('{}_I{}.csv'.format(direct, I)), delimiter=',')

        results[direct][I] = sres
        plot_search_map(sres, (zrange[0], zrange[-1]), (cvfrange[0], cvfrange[-1]), I, direct,
                        get_fname('{}_I{}.png'.format(direct, I)))

results['combined'] = {}
zs = []
for I in results['forward'].keys():
    results['combined'][I] = results['forward'][I] * results['reverse'][I]
    savetxt(get_fname('{}_I{}.csv'.format('combined', I)), results['combined'][I], delimiter=',')
    plot_search_map(results['combined'][I], (zrange[0], zrange[-1]), (cvfrange[0], cvfrange[-1]), I, 'combined',
                    get_fname('{}_I{}.png'.format('combined', I)))

    r = results['combined'][I]
    min_indicies = unravel_index(r.argmin(), r.shape)
    min_val = (zrange[min_indicies[0]], cvfrange[min_indicies[1]])
    print('I = {}: z* = {}, cvf = {}'.format(I, min_val[0], min_val[1]))
    if I != 0:
        zs.append(min_val[0])


zaverage = array(zs).mean()
print('z* average: {}'.format(zaverage))
zaverage_rounded, zaverage_index = find_nearest(zrange, zaverage)

cresults = results['combined']
#find (min, best, max) cvf at zaverage_rounded for each I
cbounds = {}
for I in cresults.keys():
    bounds_column = cresults[I][zaverage_index, :]
    bounds_truth = bounds_column < 0.0015
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
outfile = get_fname('cvfplot_combined.png')
dmplots.plot_cvf_function_ebars(plot_data, ebars, 'combined', outfile)

#and politely output simulations
cse = CalcSimExecutor(ds, T)
ce = ComparisonEngine(cse.cs)
x = linspace(0, 25, 100)
for direction in ('forward', 'reverse'):
    for I in cresults.keys():
        for idx, edge in ((0, 'lower'), (1, 'best'), (2, 'upper')):
            outfile = get_fname(str.format('Comparison_I{}_{}_{}bound.png', I, direction, edge))
            cvf = cbounds[I][1]
            if idx == 0:
                cvf -= cbounds[I][idx]
            elif idx == 2:
                cvf += cbounds[I][idx]
            print('using cvf = ' + str(cvf))
            edict = ds.edict_for_direction(direction)
            exper = ds.interpolated_experiment_dict(x, edict)[I]
            simd = cse.compute(zaverage_rounded, cvf, I, direction)
            exper = np.column_stack((x, exper))
            ce.calibrate(simd[:, 1], exper[:, 1])
            simd[:, 1] = ce.shift_data(simd[:, 1])
            f = dmplots.plot_sim_fit(simd, exper, I, zaverage_rounded, cvf, direction)
            c = FigureCanvasAgg(f)
            c.print_figure(outfile)