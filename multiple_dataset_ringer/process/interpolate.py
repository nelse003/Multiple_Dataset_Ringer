#!/usr/bin/env pandda.python

###############################################################################
# Packages
##############################################################################
# System tasks
import os

# Panda Dataframes
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)


def interpolate_all_ringer_results(ref_set, all_results, params):

    """
    Interpolate all ringer results to be on the same angle range

    Parameters
    ----------
    ref_set:
    all_results:
    params:

    Returns
    -------

    """
    datasets = all_results.keys()
    # Name for storage of interpolated results (without residue)
    interpolate_base_csv = "_{}_Datasets_{}_{}-ringer.csv".format(
        len(datasets), params.settings.map_type, params.settings.angle_type
    )

    # Iterate through the residues
    interpolated_datsets = {}
    for residue, data in ref_set.iterrows():
        residue_data_list = []

        logger.info("Interpolating ringer results for resiude {}".format(residue))

        # Create ouput directories for each residue
        if not os.path.isdir(os.path.join(params.output.out_dir, residue)):
            os.makedirs(os.path.join(params.output.out_dir, residue))
            # Output filename for interpolated data
            # (Used to check existence of output)
        interpolate_csv = residue + interpolate_base_csv

        if not os.path.exists(
            os.path.join(params.output.out_dir, residue, interpolate_csv)
        ):

            # Iterate through the datasets
            for dataset_label, dataset_results in all_results.iteritems():

                current_dataset_results = dataset_results.loc[
                    (dataset_results.index == residue)
                ]

                if current_dataset_results.empty:
                    continue

                sorted_angles, sorted_map_values = normalise_and_sort_ringer_results(
                    current_dataset_results, params=params
                )

                interpolated_angles, interpolated_map_values = linear_interpolate_ringer_results(
                    sorted_angles,
                    sorted_map_values,
                    angle_sampling=params.settings.angle_sampling,
                )

                # Store these in a list
                residue_data_list.append((interpolated_angles, interpolated_map_values))

                # If it doesn't exist: Create dataframe to store results from
                # one residue, across multiple datasets
                if "single_residue_multiple_datasets" not in locals():
                    single_residue_multiple_datasets = pd.DataFrame(
                        columns=interpolated_angles
                    )

                # Populate dataframe with results from one residue,
                # across multiple datasets
                single_residue_multiple_datasets.loc[
                    "{}".format(dataset_label)
                ] = interpolated_map_values

            # Output CSV from one resiude, multiple datasets
            pd.DataFrame.to_csv(
                single_residue_multiple_datasets,
                os.path.join(params.output.out_dir, residue, interpolate_csv),
            )
            interpolated_datsets[residue] = single_residue_multiple_datasets

        else:
            logger.info("{}: Interpolated CSVs already generated,".format(residue))

    interpolated_datsets_df = pd.concat(interpolated_datsets, axis=0)
    interpolated_datsets_df.to_csv(params.output.interp_csv)

def normalise_and_sort_ringer_results(current_dataset_results, params):
    """Sorts ringer results by angle"""

    # Extract from current rhoesidue in dataset
    residue = current_dataset_results.index[0]
    start_ang = current_dataset_results.values[0, 2]
    ang_rel = (
        params.settings.angle_sampling * current_dataset_results.columns.values[3:] - 3
    )
    map_values = current_dataset_results.values[0, 3:]

    logger.debug("Start Angle: {}".format(start_ang))
    logger.debug("Angle Relative: {}".format(start_ang))

    logger.debug("Showing data for {}".format(residue))

    # Set angles
    ang = (start_ang + ang_rel) % 360
    logger.debug("Angles: {}".format(ang))

    ###################################################
    # Sort Angles
    ##################################################
    sorted_idx = sorted(range(len(ang)), key=lambda i: ang[i])
    logger.debug("Sorted Index: {}".format(sorted_idx))

    sorted_angles = [ang[i] for i in sorted_idx]
    logger.debug("Sorted angles: {}".format(sorted_angles))

    sorted_map_values = [map_values[i] for i in sorted_idx]
    logger.debug("Sorted map values: {}".format(sorted_map_values))

    return (sorted_angles, sorted_map_values)


def linear_interpolate_ringer_results(sorted_angles, sorted_map_values, angle_sampling):

    """ Interpolate ringer results to run across same angle range """

    # Extend the map values, and angles to include the first element at end
    # (over 360 deg),and the last elment at the start (below 0 deg)
    sorted_map_values.insert(0, sorted_map_values[-1])
    # now need to append 2nd value to end of list, as 1st values is the appended value
    sorted_map_values.append(sorted_map_values[1])
    sorted_angles.insert(0, sorted_angles[0] - angle_sampling)
    sorted_angles.append(sorted_angles[-1] + angle_sampling)

    # Generate a set of angles to interpolate to based on the angle sampling
    interpolated_angles = np.arange(1, 360, angle_sampling)

    # interpolate
    interpolated_map_values = np.interp(
        interpolated_angles, sorted_angles, sorted_map_values
    )

    # Offset to set peaks [60,180,300]
    offset_map_values_end = interpolated_map_values[150:]
    offset_map_values_start = interpolated_map_values[0:150]
    interpolated_map_values = np.concatenate(
        (offset_map_values_end, offset_map_values_start)
    )

    return (interpolated_angles, interpolated_map_values)
