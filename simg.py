#!/bin/bash 
from datastore import InputDatastore 
from calcsim import CalcSimWrapper 
from numpy import *
from expercomparison import ComparisonEngine
from itertools import product, chain

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

def make_datafunc(I, direction):
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
        simd = ce.shift_data(simd) 
        return (simd, exprd)
    return rfunc

pd = product(['forward', 'reverse'], [0, 400, 800, 1000])
dfs = [make_datafunc(I, d) for (d, I) in pd]

pairs = [(421.6, 0.000131), (160, 0.53e-3), (160, 0.352e-3), (300, 0.7e-3), (200, 0.7e-3), (250, 0.7e-3), (160, 0.6e-3), (0, 0.3e-3), (0, 0.4e-3), (0, 0.5e-3), (0, 0.6e-3), (0, 0.7e-3)]

for j,(z, m) in enumerate(pairs):
    simrs = [f(z, m) for f in dfs]
    cols = [x for y in simrs for x in y] 
    rmat = zeros((100, 17))
    rmat[:,0] = x
    for i,c in enumerate(cols, start=1):
        rmat[:,i] = c
    savetxt('Pair_{}.csv'.format(j), rmat, delimiter=',')  
