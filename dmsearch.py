"""Entry point for optimum vacancy concentration search
"""

import numpy as np
import dmplots
import os
import re

from os import path
from argparse import ArgumentParser
from progressbar import ProgressBar, Percentage, Bar

from datastore import InputDatastore
from calcsim import CalcSimWrapper
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

fdata = InputDatastore(args.inputdata, args.dataprefix, args.temperature, 'forward')
rdata = InputDatastore(args.inputdata, args.dataprefix, args.temperature, 'reverse')
accelcs = CalcSimWrapper()

zrange = np.arange(args.zlim[0], args.zlim[1], args.zlim[2])
cvfrange = np.arange(args.cvflim[0], args.cvflim[1], args.cvflim[2])

fresults = dict()
rresults = dict()

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
        I = float(re.match(rx, f).groups()[0])
        rresults[I] = np.genfromtxt(f, skip_header=0, delimiter=',')
else:
    #do search for forward/back
    for direction, dstore, result_stash in [('forward', fdata, fresults), ('reverse', rdata, rresults)]:
        search_engine = ParamSearchEngine(accelcs, dstore)
        for I in dstore.experimental_dict.keys():
            print(str.format('{} bias, I = {}', direction, I))
            pbar = ProgressBar(widgets=['Parameter search: ', Percentage(), ' ', Bar()],
                               maxval=len(zrange) * len(cvfrange))
            pbar.start()
            search_map = search_engine.search(zrange, cvfrange, I, args.temperature, pbar.update)
            print('\n')
            result_stash[I] = search_map
            outpath = path.join(args.outputdir, str.format('Searchmap_I{}_{}.png', I, direction))
            dmplots.plot_search_map(search_map, (zrange[0], zrange[-1]), (cvfrange[0], cvfrange[-1]),
                                    I, direction, outpath)
            outpath = path.join(args.outputdir, str.format('Searchmap_I{}_{}.csv', I, direction))
            np.savetxt(outpath, search_map, delimiter=',')

fzaverage = 0
ict = 0

#find average for forward
for I in fresults.keys():
    if I == 0:
        continue
    min_indicies = np.unravel_index(fresults[I].argmin(), fresults[I].shape)
    min_val = (min_indicies[0] * args.zlim[2], min_indicies[1] * args.cvflim[2])
    fzaverage += min_val[0]
    ict += 1

fzaverage /= ict
print(str.format('Average z* for {} bias: {}', 'forward', fzaverage))


def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

zaverage_rounded, zaverage_index = find_nearest(zrange, fzaverage)
print(str.format('Rounding to {}', zaverage_rounded))

#now find our optimum positions
for direction, result_stash in [('forward', fresults), ('reverse', rresults)]:
    zaverage = 0
    ict = 0
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