import os
# Panda Dataframes
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

def interpolate_all_ringer_results(ref_set,
                                   all_results,
                                   angle_type,
                                   interpolate_base_csv,
                                   params):

    """Put all Ringer results on a common angle scale by interpolation
    
    Parameters
    -------------
    ref_set: pandas.DataFrame
        A dataframe conatining the reference dataset.
    all_results: dict
        Dictionary of dataframes containing ringer results, 
        labelled per dataset 
    params:
        python object derived from master phil file.
        contains the 
        
    Returns
    ----------
    None
    """
    datasets = all_results.keys()

    # Iterate through the residues
    for residue, data in ref_set.iterrows():
        residue_data_list = []

        logger.info("Interpolating ringer results for resiude {}".format(residue))

        # Create ouput directories for each residue
        if not os.path.isdir(os.path.join(params.output.out_dir, residue)):
            os.makedirs(os.path.join(params.output.out_dir, residue))
            # Output filename for interpolated data
            # (Used to check existence of output)
        interpolate_csv = residue + interpolate_base_csv
        print(interpolate_csv)
        if not os.path.exists(os.path.join(params.output.out_dir,
                                           residue,
                                           interpolate_csv)):

            # Iterate through the datasets
            for dataset_label, dataset_results in all_results.iteritems():

                current_dataset_results = dataset_results.loc[(
                    dataset_results.index == residue)]

                if current_dataset_results.empty:
                    continue

                sorted_angles, sorted_map_values = \
                    normalise_and_sort_ringer_results(current_dataset_results,
                                                      angle_type,
                                                      params=params)
                interpolated_angles, \
                interpolated_map_values = linear_interpolate_ringer_results(
                    sorted_angles,
                    sorted_map_values,
                    angle_sampling=params.settings.angle_sampling)

                # Store these in a list
                residue_data_list.append((interpolated_angles,
                                          interpolated_map_values))

                # If it doesn't exist: Create dataframe to store results from
                # one residue, across multiple datasets
                if'single_residue_multiple_datasets' not in locals():
                    single_residue_multiple_datasets = pd.DataFrame(
                        columns=interpolated_angles)

                # Populate dataframe with results from one residue,
                # across multiple datasets
                single_residue_multiple_datasets.loc['{}'.format(
                    dataset_label)] = interpolated_map_values

            # Output CSV from one resiude, multiple datasets
            pd.DataFrame.to_csv(single_residue_multiple_datasets,
                                    os.path.join(params.output.out_dir, residue,
                                                 interpolate_csv))
        else:
            logger.info(
                '{}: Interpolated CSVs already generated,'.format(residue))


def normalise_and_sort_ringer_results(current_dataset_results, angle_type, params):
    """Sorts ringer results by angle"""

    # Extract from current residue in dataset
    residue = current_dataset_results.index[0]

    for ang_dataset in current_dataset_results.iterrows():
        if ang_dataset[1][2] == angle_type:
            start_ang  = ang_dataset[1][3]
            map_values = ang_dataset[1].values[3:]
    ang_rel = params.settings.angle_sampling * (current_dataset_results.columns.values[3:]-3)


    logger.debug("Start Angle: {}".format(start_ang))
    logger.debug("Angle Relative: {}".format(start_ang))
    logger.debug('Showing data for {}'.format(residue))

    # Set angles
    ang = (start_ang+ang_rel)%360
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


def linear_interpolate_ringer_results(sorted_angles,
                                      sorted_map_values,
                                      angle_sampling):

    """ Interpolate ringer results to run across same angle range """

    # Extend the map values, and angles to include the first element at end 
    # (over 360 deg),and the last elment at the start (below 0 deg)
    sorted_map_values.insert(0,sorted_map_values[-1])
    # now need to append 2nd value to end of list, as 1st values is the appended value
    sorted_map_values.append(sorted_map_values[1])
    sorted_angles.insert(0,sorted_angles[0]-angle_sampling)
    sorted_angles.append(sorted_angles[-1]+angle_sampling)

    # Generate a set of angles to interpolate to based on the angle sampling
    interpolated_angles = np.arange(1, 360, angle_sampling)

    # interpolate
    interpolated_map_values = np.interp(interpolated_angles,
                                           sorted_angles,
                                           sorted_map_values)

    # Offset to set peaks [60,180,300]
    offset_map_values_end = interpolated_map_values[150:]
    offset_map_values_start = interpolated_map_values[0: 150]
    interpolated_map_values = np.concatenate((offset_map_values_end,
                                               offset_map_values_start))

    return (interpolated_angles, interpolated_map_values)

