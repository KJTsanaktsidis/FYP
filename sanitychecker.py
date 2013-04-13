from calcsim import CalcSimWrapper
from datastore import InputDatastore
import numpy as np
from pylab import *

ds = InputDatastore('../InputData', 'NiCu', 973)
cs = CalcSimWrapper()

D = ds.interpolated_diffusivity(10001)
R = ds.interpolated_resistivity(10001)

dx = 0.5*35e-8
ndx = 200
dt = cs.optimum_dt(dx, D, 1)
ndt = cs.num_sim_steps(dt, 2 * 60 * 60)

init = np.ones(ndx)
init[ndx/2:] = 0

res = cs.calc_simulation(D, R, init, ndt, dt, dx, 0, 1)
x = np.linspace(0, dx * ndx, num=ndx)
plot(x, res)
show()