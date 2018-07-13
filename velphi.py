
import numpy as np
import os
import sys
import h5py
import scidata
import scidata.carpet.hdf5 as h5
import matplotlib.pyplot as plt
from math import sqrt
from scidata.utils import locate
import argparse

parser = argparse.ArgumentParser(description=
        "Plot the angular velocity of Carpet data \
         around the grid's center.")
parser.add_argument("plane", type=str, help=
        "Plane centered at the origin. 'all' plots all three planes.")

args = parser.parse_args()

if len(args) == 0:
    parser.error("No plane specified")

plane = args.plane
planes = ["xy", "xz", "yz"]
if plane=='all':
    sys.stderr.write("Plotting for all planes ... \n")
    for p in planes:
        plot(p)
else:
    assert plane in planes
    plot(plane)
sys.stderr.write("All done!")

def plot(plane):
    # Find the velocity files
    os.chdir('../')
    xlist = locate('vel[0].{}.h5'.format(plane))
    print(xlist)
    ylist = locate('vel[1].{}.h5'.format(plane))
    print(ylist)
    os.chdir('analysis')

    reflevels = (0,1,2,3,4,5,6)

    # Read the metadata
    xset = h5.dataset(xlist)
    yset = h5.dataset(ylist)

    sys.stderr.write("velphi {}-plane\n".format(plane))
    for reflevel in reflevels:
        sys.stderr.write("Plotting reflevel {} ... ".format(reflevel))
        # Find the max and min velocities
        vamax = []
        vamin = []
        for it in xset.iterations:
            va, = get_velphi(xset, yset, it, reflevel)
            vamax.append(np.max(va))
            vamin.append(np.min(va))
        vamax = max(vamax)
        vamin = min(vamin)

        # Plot each iteration
        for it in xset.iterations:
            va, dx, dy = get_velphi(xset, yset, it, reflevel)

            # Set the tick marks for the axes
            tics = get_ticks(reflevel)
            X = center[0] + (tics / dx)
            Y = center[1] + (tics / dy)

            # Plot
            plt.figure(figsize=(6,6))
            plt.imshow(va, interpolation=None,  origin='lower',
                        vmin=vamin, vmax=vamax)
            cbar = plt.colorbar()
            cbar.set_label(r'$\omega$')
            plt.clim=(0.0, np.max(va))
            if plane=='xy':
                plt.xlabel(r'x [M$_{\odot}$]')
                plt.ylabel(r'y [M$_{\odot}$]')
            elif plane=="xz":
                plt.xlabel(r'x [M$_{\odot}$]')
                plt.ylabel(r'z [M$_{\odot}$]')
            else:
                plt.xlabel(r'y [M$_{\odot}$]')
                plt.ylabel(r'z [M$_{\odot}$]')
            plt.xticks(X, tics)
            plt.yticks(Y, tics)

        sys.stderr.write("done!\n")

def get_velphi(xset, yset, it, reflevel):
    xgrid = xset.get_reflevel(iteration=it, reflevel=reflevel)
    ygrid = yset.get_reflevel(iteration=it, reflevel=reflevel)

    vx = xset.get_reflevel_data(xgrid, iteration=it, timelevel=0)
    vy = yset.get_reflevel_data(ygrid, iteration=it, timelevel=0)

    assert vx.shape==vy.shape, 'Shapes are not equal'
    row = np.size(vx,0)
    col = np.size(vx,1)
    center = np.array([row, col])*0.5 # Grid center

    assert np.all([xgrid.delta,ygrid.delta]), 'Grid spacing is not equal'
    dx = xgrid.delta[0]
    dy = ygrid.delta[1]

    va = np.full((row,col), np.NAN, dtype=np.float32)

    # Calculate the angular velocity at each point
    for i in range(row):
        for j in range(col):
            r = sqrt((dx*(i - center[0]))**2 + (dy*(j - center[1]))**2)
            if r==0:
                r += 1.0e-10 # Avoid division by zero
            va[i][j] = sqrt(vx[i][j]**2 + vy[i][j]**2) / r

    # Eliminate the 'infinity' at the center
    va[int(center[0])][int(center[1])] = 0.0

    return va, dx, dy

def get_ticks(reflevel):
    if reflevel==0:
        tic = np.array([250, 500, 750, 1000])
    elif reflevel==1:
        tic = np.array([75, 150, 225, 300, 375])
    elif reflevel==2:
        tic = np.array([50, 100, 150, 200])
    elif reflevel==3:
        tic = np.array([25, 50, 75, 100])
    elif reflevel==4:
        tic = np.array([10, 20, 30, 40, 50])
    elif reflevel==5:
        tic = np.array([10, 20])
    elif reflevel==6:
        tic = np.array([4, 8, 12, 14])
    else:
        raise Exception('This is not a ticked reflevel')

    return np.int_(np.concatenate((-tics[::-1], [0], tics)))
