#!/bin/bash

from datastore import InputDatastore
from calcsim import CalcSimWrapper 
from numpy import *
from scipy.optimize import *
from expercomparison import ComparisonEngine
from itertools import product

ds = InputDatastore('../InputData', 'NiCu')
x = linspace(0, 25, num=100)
fedict = ds.edict_for_direction('forward')
redict = ds.edict_for_direction('reverse')
fexpr = ds.interpolated_experiment_dict(x, fedict)
rexpr = ds.interpolated_experiment_dict(x, redict)
diffusivity = ds.interpolated_diffusivity(1001, 973)
resistivity = ds.interpolated_resistivity(1001, 973)

cs = CalcSimWrapper()
ce = ComparisonEngine(cs)
initcond = ones(100)
initcond[50:] = 0


dt = 0.05
ndt = int(2 * 60 * 60 / 0.05)
dx = 25e-6 / 100

def make_objective(I, direction):
    if direction == 'forward':
        exprd = fexpr[I]
    else:
        exprd = rexpr[I]
        I = - I
    
    def rfunc(z, m):
        cvf = 1 + m * I
        r = cs.emigration_factor(z, I*100*100, 973)
        simd = cs.calc_simulation(diffusivity, resistivity, initcond, ndt, dt, dx, r, cvf)
        lsq, shift = ce.calibrate(simd, exprd)
        return lsq
    return rfunc

objs = [make_objective(I, d) for (I, d) in product([0, 400, 800, 1000], ['forward', 'reverse'])]

def master_objective(tup):
    z = tup[0]
    m = tup[1]
    lsqs = [f(z, m) for f in objs]
    return sum(lsqs)

optr = minimize(master_objective, (160, 0.3e-3), method='TNC', bounds=[(10, 800), (0, 0.7e-3)], options={'maxiter' : 100})

print('Result: z* = {}, m = {}'.format(optr.x[0], optr.x[1]))
print('Termination reason: {}'.format(optr.message))

    
