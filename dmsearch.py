"""Entry point for optimum vacancy concentration search
"""
from matplotlib.backends.backend_agg import FigureCanvasAgg

import numpy as np
import dmplots
import os
import re

from os import path
from argparse import ArgumentParser
from progressbar import ProgressBar, Percentage, Bar

from datastore import InputDatastore
from calcsim import CalcSimWrapper
from expercomparison import ComparisonEngine
from paramsearch import ParamSearchEngine


aparser = ArgumentParser(description='Searches for the optimum vacancy concentration multiplier')
aparser.add_argument('--inputdata', metavar='DIR', type=str, required=True,
                     help='Directory containing input data')
aparser.add_argument('--outputdir', metavar='DIR', type=str, default='.',
                     help='Directory to stash the output results (plots, tables etc.)')
aparser.add_argument('--dataprefix', metavar='PREFIX', type=str, default='NiCu',
                     help='Prefix to data files')
aparser.add_argument('--temperature', type=float, default=973,
                     help='Temperature in K')
aparser.add_argument('--zlim', metavar=('ZMIN', 'ZMAX', 'ZSTEP'), nargs=3, type=float, required=True,
                     help='Limits of the effective valence search (up to but not including ZMAX)')
aparser.add_argument('--cvflim', metavar=('CVFMIN', 'CVFMAX', 'CVFSTEP'), nargs=3, type=float, required=True,
                     help='Limits of the vacancy concentration multiplier search (up to but not including CVFSTEP)')
aparser.add_argument('--resume', action='store_true', default=False,
                     help='Resume using .csv files in outputdir')

args = aparser.parse_args()

dstore = InputDatastore(args.inputdata, args.dataprefix)
accelcs = CalcSimWrapper()

zrange = np.arange(args.zlim[0], args.zlim[1], args.zlim[2])
cvfrange = np.arange(args.cvflim[0], args.cvflim[1], args.cvflim[2])

fresults = dict()
rresults = dict()

np.savetxt(path.join(args.outputdir, 'zrange.csv'), zrange, delimiter=',')
np.savetxt(path.join(args.outputdir, 'cvfrange.csv'), cvfrange, delimiter=',')

if args.resume:
    files = [path.join(args.outputdir, f) for f in os.listdir(args.outputdir)
             if path.isfile(path.join(args.outputdir, f))]
    forwardfiles = [f for f in files if f.endswith('forward.csv')]
    reversefiles = [f for f in files if f.endswith('reverse.csv')]
    rx = r'.*Searchmap_I([\-0-9\.]+)_.*.csv'
    for f in forwardfiles:
        print(str.format('Processing {}', f))
        I = float(re.match(rx, f).groups()[0])
        fresults[I] = np.genfromtxt(f, skip_header=0, delimiter=',')
    for f in reversefiles:
        print(str.format('Processing {}', f))
        I = -float(re.match(rx, f).groups()[0])
        rresults[I] = np.genfromtxt(f, skip_header=0, delimiter=',')
else:
    #do search for forward/back
    for direction, result_stash in [('forward', fresults), ('reverse', rresults)]:
        search_engine = ParamSearchEngine(accelcs, dstore)
        edict = dstore.edict_for_direction(direction)
        for I in edict.keys():
            print(str.format('{} bias, I = {}', direction, I))
            pbar = ProgressBar(widgets=['Parameter search: ', Percentage(), ' ', Bar()],
                               maxval=len(zrange) * len(cvfrange))
            pbar.start()
            search_map = search_engine.search(zrange, cvfrange, I, args.temperature, pbar.update, direction)
            print('\n')
            result_stash[I] = search_map
            outpath = path.join(args.outputdir, str.format('Searchmap_I{}_{}.png', I, direction))
            dmplots.plot_search_map(search_map, (zrange[0], zrange[-1]), (cvfrange[0], cvfrange[-1]),
                                    I, direction, outpath)
            outpath = path.join(args.outputdir, str.format('Searchmap_I{}_{}.csv', I, direction))
            np.savetxt(outpath, search_map, delimiter=',')

fzaverage = 0
ict = 0
zvalues = []
cresults = {}

#multiply forward/reverse maps to get current-min maps
for I in fresults.keys():
    combined_map = fresults[I] * rresults[I]
    outpath = path.join(args.outputdir, str.format('Searchmap_I{}_combined.png', I))
    dmplots.plot_search_map(combined_map, (zrange[0], zrange[-1]), (cvfrange[0], cvfrange[-1]), I,
                            'combined', outpath)
    cresults[I] = combined_map
    if I == 0:
        continue
    min_indicies = np.unravel_index(combined_map.argmin(), combined_map.shape)
    min_z = zrange[min_indicies[0]]
    print(str.format('Combined z* for I = {}: {}', I, min_z))
    zvalues.append(min_z)

fzaverage = np.array(zvalues).mean()
print(str.format('Average z*: {}', fzaverage))


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

zaverage_rounded, zaverage_index = find_nearest(zrange, fzaverage)
print(str.format('Rounding to {}', zaverage_rounded))

#now find our optimum positions
for direction, result_stash in [('forward', fresults), ('reverse', rresults), ('combined', cresults)]:
    for I in result_stash.keys():
        min_indicies = np.unravel_index(result_stash[I].argmin(), result_stash[I].shape)
        min_val = (zrange[min_indicies[0]], cvfrange[min_indicies[1]])
        print(str.format('{} bias, I = {}: best (z*, cvf) = ({}, {})', direction, I,
              min_val[0], min_val[1]))

    cvf_plotlist = []
    I_plotlist = []
    for I in result_stash.keys():
        cvf_best = cvfrange[result_stash[I][zaverage_index, :].argmin()]
        cvf_plotlist.append(cvf_best)
        I_plotlist.append(np.abs(I))

    cvf_plotarr = np.array(cvf_plotlist)
    I_plotarr = np.array(I_plotlist)
    print(str.format('Plotting I = [{}], Cvf = [{}]', I_plotarr, cvf_plotarr))
    plotarr = np.column_stack((I_plotarr, cvf_plotarr))
    outfname = os.path.join(args.outputdir, str.format('cvfplot_{}.png', direction))
    dmplots.plot_cvf_function(plotarr, direction, outfname)

    #let's plot z* as well
    z_plotlist = []
    for I in I_plotlist:
        z_best_idx = np.unravel_index(result_stash[I].argmin(), result_stash[I].shape)[1]
        z_best = zrange[z_best_idx]
        z_plotlist.append(z_best)

    z_plotarr = np.array(z_plotlist)
    plotarr = np.column_stack((I_plotarr, z_plotarr))
    outfname = os.path.join(args.outputdir, str.format('zplot_{}.png', direction))
    dmplots.plot_z_function(plotarr, direction, outfname)

    #we should be nice and print comparison plots
    if direction == 'forward' or direction == 'reverse':
        x = np.linspace(0, 25, num=100)
        dstore = InputDatastore(args.inputdata, args.dataprefix)
        edict = dstore.edict_for_direction(direction)
        accelcs = CalcSimWrapper()
        ce = ComparisonEngine(accelcs)
        diffusivity = dstore.interpolated_diffusivity(10001, args.temperature, precise=True)
        resistivity = dstore.interpolated_resistivity(10001, args.temperature)
        init_cond = np.ones(100)
        init_cond[50:] = 0
        emigration_T = args.temperature
        dt = 0.05
        ndt = int(2 * 60 * 60 / 0.05)
        dx = 25e-6 / 100
        for I in result_stash.keys():
            cvf_best = cvfrange[result_stash[I][zaverage_index, :].argmin()]
            if direction == 'forward':
                exper_data = dstore.interpolated_experiment_dict(x, edict)[I]
                r = accelcs.emigration_factor(zaverage_rounded, I * 100 * 100, emigration_T)
            else:
                exper_data = dstore.interpolated_experiment_dict(x, edict)[I]
                r = accelcs.emigration_factor(zaverage_rounded, -I * 100 * 100, emigration_T)
            simd = accelcs.calc_simulation(diffusivity, resistivity, init_cond, ndt, dt, dx, r, cvf_best)
            lsq, shift = ce.calibrate(simd, exper_data)
            shifted_simd = ce.shift_data(simd)
            full_simd = np.column_stack((x, shifted_simd))
            full_exper = np.column_stack((x, exper_data))
            outfname = os.path.join(args.outputdir, str.format('SimExperComp_I{}_{}.png', I, direction))
            f = dmplots.plot_sim_fit(full_simd, full_exper, I, zaverage_rounded, cvf_best, direction)
            c = FigureCanvasAgg(f)
            c.print_figure(outfname)

