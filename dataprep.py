"""Prepares a complete experimental .csv file from parts
"""

from argparse import ArgumentParser
import numpy as np
from os import path
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt


aparser = ArgumentParser(description='Prepares a complete experimental .csv file from parts')
aparser.add_argument('directory', metavar='DIR', type=str,
                     help='Directory containing constituent .csv files')
aparser.add_argument('--currents', dest='currents', type=int, nargs='+', default=[0],
                     help='Current values to attempt to read')
aparser.add_argument('--plot', dest='plot', action='store_true', default=False,
                     help='Whether or not data should be plotted')


args = aparser.parse_args()
x = np.linspace(0, 25, 100)
f_smoothed = np.zeros((101, len(args.currents) + 1))
r_smoothed = np.zeros((101, len(args.currents) + 1))


def process_file(direction, I):
    filename = path.join(args.directory, str.format('{}_{}.csv', direction, I))
    data = np.genfromtxt(filename, skip_header=1, delimiter=',')
    data[np.where(data[:, 1] > 1), 1] = 1
    data[np.where(data[:, 1] < 0), 1] = 0
    front_x = np.arange(0, 5, 0.5)
    back_x = np.arange(20.5, 25.5, 0.5)
    tmp_x = np.concatenate((front_x, data[:, 0] + 5, back_x))
    tmp_y = np.concatenate((np.ones(10), data[:, 1], np.zeros(10)))
    sp = InterpolatedUnivariateSpline(tmp_x, tmp_y)
    spx = sp(x)
    spx[np.where(spx < 0)] = 0
    spx[np.where(spx > 1)] = 1
    return spx

f_smoothed[1:, 0] = x
r_smoothed[1:, 0] = x
#try and open each forward/reverse file
for idx, I in enumerate(args.currents):
    f_smoothed[1:, idx + 1] = process_file('forward', I)
    f_smoothed[0, idx + 1] = I
    r_smoothed[1:, idx + 1] = process_file('reverse', I)
    r_smoothed[0, idx + 1] = I

    if args.plot:
        f = plt.figure()
        plt.plot(f_smoothed[1:, 0], f_smoothed[1:, idx + 1])
        plt.xlabel('x (micron)')
        plt.ylabel('Concentration')
        plt.xlim((0, 25))
        plt.ylim((0, 1))
        plt.title('I = ' + str(I) + ' forward')
        plt.show()

        f = plt.figure()
        plt.plot(f_smoothed[1:, 0], r_smoothed[1:, idx + 1])
        plt.xlabel('x (micron)')
        plt.ylabel('Concentration')
        plt.xlim((0, 25))
        plt.ylim((0, 1))
        plt.title('I = ' + str(I) + ' reverse')
        plt.show()

np.savetxt(path.join(args.directory, 'forward.csv'), f_smoothed, delimiter=',')
np.savetxt(path.join(args.directory, 'reverse.csv'), r_smoothed, delimiter=',')