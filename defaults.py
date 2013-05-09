"""This file contains default values for some simulation-related parameters that I don't
anticipate needing much control over
"""

simulation_xsteps = 101
simulation_dx = 25e-6 / 100
simulation_tsteps = int(2 * 60 * 60 / 0.05)
simulation_dt = 0.05