import matplotlib
matplotlib.use("Agg")
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
import re
import tempfile

parser = argparse.ArgumentParser(description=
        "Plot the angular velocity of Carpet data \
         around the grid's center.")
parser.add_argument("plane", type=str, help=
        "Plane centered at the origin. 'all' plots all three planes.")
parser.add_argument("--reflevel", default='all',
        help="Choose reflevel to plot. Default are 0-6.")

args = parser.parse_args()

def plot(plane, ref):
    # Create plane directory if non-existing
    if not os.path.exists(plane):
        os.mkdir(plane)

    # Find the velocity files
    ret = os.getcwd()
    os.chdir('../../')
    xlist = locate('vel[0].{}.h5'.format(plane))
    #print(xlist)
    ylist = locate('vel[1].{}.h5'.format(plane))
    #print(ylist)
    os.chdir(ret)

    if ref=='all':
        reflevels = (0,1,2,3,4,5,6)
    else:
        reflevels = [int(ref)]

    # Read the metadata
    xset = h5.dataset(xlist)
    yset = h5.dataset(ylist)

    sys.stderr.write("\nvelphi {}-plane\n".format(plane))

    # Find the max and min velocities of the finest reflevel
    sys.stderr.write("Finding the min and max ... ")
    vamax = []
    vamin = []
    for it in xset.iterations:
        va = get_velphi(xset, yset, it, 6)[0]
        vamax.append(np.max(va))
        vamin.append(np.min(va))
    vamax = max(vamax)
    vamin = min(vamin)
    sys.stderr.write("done!\n")

    for reflevel in reflevels:
        sys.stderr.write("Reflevel {}: ".format(reflevel))
        pdir = "{0}/r{1}".format(plane,reflevel)
        # Make directory for reflevel
        if not os.path.exists(pdir):
            os.mkdir(pdir)

        # Plot each iteration
        sys.stderr.write("Plotting ... ")
        for it in xset.iterations:
            va, center, dx, dy = get_velphi(xset, yset, it, reflevel)

            # Set the tick marks for the axes
            tics = get_ticks(reflevel)
            X = center[0] + (tics / dx)
            Y = center[1] + (tics / dy)

            # Plot
            plt.figure(figsize=(6,6))
            plt.imshow(va.T, interpolation=None,  origin='lower',
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

            t = xset.get_dataset(iteration=it).attrs["time"]
            plt.title("time = %10.2f" %t)

            plt.savefig("{0}/{1}.png".format(pdir, it))

        # Make movie
        sys.stderr.write("making movie ... ")

        os.chdir(pdir)
        ldir   = tempfile.mkdtemp()
        fnames = os.listdir('./')
        fnames = [f for f in fnames if re.match(r".*\.png$", f) is not None]
        fnames = sorted(fnames)

        for i in range(len(fnames)):
            os.symlink(os.getcwd() + "/" + fnames[i],
                    ldir + "/%05d.png" % i)
        os.system("ffmpeg -y -i {}/%05d.png -vcodec mpeg4 -qscale 1 movie.mp4"
            " &> /dev/null".format(ldir))
        for i in range(len(fnames)):
            os.unlink(ldir + "/%05d.png" % i)
        os.rmdir(ldir)
        os.chdir(ret)

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

    return va, center, dx, dy

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

    return np.int_(np.concatenate((-tic[::-1], [0], tic)))

plane = args.plane
ref = args.reflevel
planes = ["xy", "xz", "yz"]
if plane=='all':
    sys.stderr.write("Plotting for all planes ... \n")
    for p in planes:
        plot(p,ref)
else:
    assert plane in planes
    plot(plane, ref)
sys.stderr.write("All done!\n")
