#!/usr/bin/env pandda.python

###############################################################################
# Packages
##############################################################################
# System tasks
import os, sys, copy, glob
# Panda Dataframes
import pandas

import numpy
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh = logging.FileHandler('ringer_script.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


def normalise_and_sort_ringer_results(current_dataset_results,params):
    """Sorts ringer results by angle"""

    # Extract from current residue in dataset

    residue = current_dataset_results.index[0]
    start_ang  = current_dataset_results.values[0,2]
    ang_rel = params.settings.angle_sampling * current_dataset_results.columns.values[3:]-3
    map_values = current_dataset_results.values[0,3:]

    logger.info('Showing data for {}'.format(residue))

    # Set angles
    ang = (start_ang+ang_rel)%360

    ###################################################
    # Sort Angles
    ##################################################
    sorted_idx = sorted(range(len(ang)), key=lambda i: ang[i])
    sorted_angles = [ang[i] for i in sorted_idx]
    sorted_map_values = [map_values[i] for i in sorted_idx]

    return (sorted_angles, sorted_map_values)


def linear_interpolate_ringer_results(sorted_angles,sorted_map_values,angle_sampling):
    """ Interpolate ringer results to run across same angle range """

    # Extend the map values, and angles to include the first element at end 
    # (over 360 deg),and the last elment at the start (below 0 deg)
    sorted_map_values.insert(0,sorted_map_values[-1])
    # now need to append 2nd value to end of list, as 1st values is the appended value
    sorted_map_values.append(sorted_map_values[1])
    sorted_angles.insert(0,sorted_angles[0]-angle_sampling)
    sorted_angles.append(sorted_angles[-1]+angle_sampling)

    # Generate a set of angles to interpolate to based on the angle sampling
    interpolated_angles = numpy.arange(1, 360, angle_sampling)

    # interpolate
    interpolated_map_values = numpy.interp(interpolated_angles,sorted_angles,sorted_map_values)

    # Offset to set peaks [60,180,300]
    offset_map_values_end = interpolated_map_values[150:]
    offset_map_values_start = interpolated_map_values[0: 150]
    interpolated_map_values = numpy.concatenate((offset_map_values_end,
                                               offset_map_values_start))

    return (interpolated_angles,interpolated_map_values)

