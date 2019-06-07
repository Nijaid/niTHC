import argparse
parser = argparse.ArgumentParser(description=
         "Plot rho for all segments.")
parser.add_argument("--skip", "-s", nargs="*", type=str,
        help="Skip specific segments (format: ####).")
args = parser.parse_args()

import sys
import os
if all(os.getcwd() != s for s in sys.path):
    sys.path = [os.getcwd()] + sys.path

from shutil import copy

import config.base
import config.rho

import scivis.config
scivis.config.validate()

from scivis.vis.field import ColorMap2D
from scivis.vis.vector import VectorColorMap2D
from scivis.data.carpet2d import ScalarData2D, VectorData2D
from scivis.fields.hydro import Density, MassFlux
from scivis.utils import mkdir, printer
from scivis import movie

# Turn off making movies in each segment
movie_bool = scivis.config.options['base.makemovies']
scivis.config.options['base.makemovies'] = False

# Get output directories
from glob import glob
segments = [x[-4:] for x in sorted(glob("../output-????"))]
segments = [x for x in segments if x not in args.skip]

# Plot each segment
for seg in segments:
    print("*** %s ***\n"%seg)

    # Make output directory structure for segment
    # mkdir(data_dir + "uber_analysis/")
    # mkdir(data_dir + "uber_analysis/output")
    # mkdir(data_dir + "uber_analysis/log")

    scivis.config.options['base.datapath'] = "../output-%s/data/"%seg

    if scivis.config.options['rho.massflux']:
        VectorColorMap2D(
                ScalarData2D(Density()),
                VectorData2D(MassFlux()),
                "rho_%s"%seg
        ).run()
    else:
        ColorMap2D(ScalarData2D(Density()), "rho_%s"%seg).run()

# Setup directories for all images
planes = scivis.config.options['plot.planes']
for p in planes:
    mkdir("output/rho/%s"%p)

root = 'rho'
if scivis.config.options['plot.scale'] == 'log':
    root = 'log_' + root
elif scivis.config.options['plot.scale'] == 'log_abs':
    root = 'log_abs_' + root

# Copy each segment's images into a single directory
# for each reflevel that is included across all segments.
for p in planes:
    print("For %s:"%p)
    level_lis = []
    for seg in segments:
        levels = glob("output/{0}_{1}/{2}/r?".format(root, seg, p))
        level_lis.append(set([os.path.basename(x) for x in levels]))

    # Find the levels covered by all segments
    uber_levels = sorted(list(set.intersection(*level_lis)))
    print("Copying levels: ")
    print(uber_levels)
    print("\n")

    for l in uber_levels:
        level_dir = "output/%s/%s/%s"%(root, p, l)
        mkdir(level_dir)
        for seg in segments:
            copy("output/{0}_{1}/{2}/{3}/*.png".format(root, seg, p, l),
                        level_dir)
        if movie_bool:
            printer.progress("---> {0} {1} : {2}/movie.mp4".format(root, p, level_dir))
            movie(level_dir)
            printer.final("---> {0} {1} : movies done!".format(root, p))
