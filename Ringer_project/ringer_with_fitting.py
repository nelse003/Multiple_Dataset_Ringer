

###############################################################################
# Packages
##############################################################################
import os
import sys
import glob
import pandas
import libtbx.phil
import numpy
import logging
#################################
import matplotlib

matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot

pyplot.style.use('ggplot')

from string import ascii_letters
##################################

###########################################################
# Function imports
##########################################################

# Ringer Processing with absolute electron density scaling (F000)
from process_ringer import process_with_ringer
# Sorting & Interpolating angles from ringer output
from interpolate_ringer import (linear_interpolate_ringer_results,
                                normalise_and_sort_ringer_results)
# Plotting functions
from plotting_ringer import (line_plot_ringer,
                             multiple_line_plot_ringer,
                             average_ringer_plots,
                             plot_correlation_vs_fitting,
                             plot_resloution_vs_dataset_score,
                             plot_RMSD_vs_dataset_score,
                             plot_RMSD_vs_resolution,
                             plot_RMSD_vs_residue_score,
                             peak_angle_histogram,
                             cluster_heatmap,
                             pairwise_heatmap)

# Hierarichal clustering functions
from hier_cluster_ringer import (hier_agg_cluster, find_pairwise_range)

# Fitting functions
from fitting_ringer import (fit_all_datasets, generate_RMSD,
                           calculate_euclidean_distance)

# Correlation functions
from correlation_ringer import correlation_single_residue

# TODO Add back in and test peak finding
# Peak finding
# from peak_count import find_peaks

###############################################################################
# Set up for passing arguments
############################################################################### 

blank_arg_prepend = {None: 'dir=', '.pdb': 'pdb=', '.mtz': 'mtz='}

master_phil = libtbx.phil.parse("""
input {
    dir = None
        .type = path
        .multiple = True
    pdb_style = "dimple.pdb"
        .type = str
        .multiple = False
    mtz_style = "dimple.mtz"
      .type = str
        .multiple = False
}
output {
    log = "ringer.log"
        .type = str
        .multiple = False
    out_dir = "output"
        .type = str
        .multiple = False
}
settings {
    # XXX mmtbx.ringer can only take this an integer, >1 XXX#
    angle_sampling = 2
        .type = int
        .multiple = False

    gen_heatmap = False
        .type = bool
        .multiple = False
}
""")

###############################################################################
# Logging
###############################################################################
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


################################################################################
#                                   FUNCTIONS                                  #
################################################################################

##############################################################################
# Main program 
###############################################################################

def run(params):

    """

    Notes
    ------------------
    Needs to be run with a input (params.input.dir) directory organised with
    pdb files in /dir/*pdb_style*
    
    Process ringer uses 2FOFCWT, PHI2FOFCWT columns. 
    Handle other columns by moving to input check.

    Example
    ------------------

    """

    # Dictionary to store all of the
    # ringer results for each of the 
    # datasets
    all_results = {}

    # Create an output directory if it doesn't already exist
    if not os.path.isdir(params.output.out_dir):
        os.makedirs(params.output.out_dir)

    # Resolution CSV filename & path
    resolution_csv_path = os.path.join(params.output.out_dir, 'dataset_resolution.csv')
    # DataFrame to store resolution
    dataset_resolution = pandas.DataFrame(index=params.input.dir,
                                          columns=['Resolution'])

    # Generate ringer results & resolution information
    dataset_counter = 0
    for dataset_dir in params.input.dir:
        # Label the dataset by the directory name
        dataset_label = os.path.basename(dataset_dir.rstrip('/'))
        pdb = glob.glob(os.path.join(dataset_dir, params.input.pdb_style))
        mtz = glob.glob(os.path.join(dataset_dir, params.input.mtz_style))

        if not pdb:
            continue

        if not mtz:
            continue

        pdb = pdb[0]
        mtz = mtz[0]

        if not os.path.exists(pdb):
            print('Skipping dir: No PDB Files found in {} matching {}'.format(
                dataset_dir, params.input.pdb_style))
            continue

        if not os.path.exists(mtz):
            print('Skipping dir: No MTZ Files found in {} matching {}'.format(
                dataset_dir, params.input.mtz_style))
            continue

        dataset_counter += 1

        # Process dataset with ringer and convert results to DataFrame
        ringer_csv, resolution = \
            process_with_ringer(pdb=pdb,
                                mtz=mtz,
                                angle_sampling=params.settings.angle_sampling,
                                output_dir=dataset_dir,
                                resolution_csv_path=resolution_csv_path)

        ringer_results = pandas.DataFrame.from_csv(ringer_csv, header=None)

        # Change order of residue name.
        # This is needed for the comparison:
        # current_dataset_results = dataset_results.loc[(dataset_results.index == residue)]
        # ~ line 238

        for i in range(0, len(ringer_results.index)):
            res_split = ringer_results.index.values[i].rsplit(' ')
            if len(res_split) > 2:
                ringer_results.index.values[i] = res_split[0] + ' ' + res_split[1]

        all_results[dataset_label] = ringer_results
        dataset_resolution.loc[dataset_label] = resolution

    print(dataset_counter)
    print("###########################################")
    # Resolution to CSV
    if not os.path.exists(resolution_csv_path):
        dataset_resolution.to_csv(resolution_csv_path)
        logger.info('Resolution CSV generated')
    else:
        logger.info('Resolution CSV already generated')

    # Pull out the "first" ringer results set as a reference
    ref_set = all_results.values()[0]

    # Map and angle types currently selected to analyse
    map_type = '2mFo-DFc'
    angle_type = 'chi1'

    # Name for storage of interpolated results (without residue)
    # TODO Need to swap out length with more appropriate length,
    # TODO as datasets can be skipped if no pdb/mtz
    interpolate_base_csv = '_{}_Datasets_{}_{}-ringer.csv'.format(len(params.input.dir), map_type, angle_type)

    # Choose a map_type/angle_type by reducing reference set
    ref_set = ref_set.loc[(ref_set[1] == map_type)]
    ref_set = ref_set.loc[(ref_set[2] == angle_type)]

    # Iterate through the residues
    for residue, data in ref_set.iterrows():
        residue_data_list = []

        # Create ouput directories for each residue
        if not os.path.isdir(os.path.join(params.output.out_dir, residue)):
            os.makedirs(os.path.join(params.output.out_dir, residue))
            # Output filename for interpolated data
            # (Used to check existence of output)
        interpolate_csv = residue + interpolate_base_csv

        # Output filename for correlation data
        correlation_csv = '{}_from {} datasets-correlation-ringer.csv'.format(
            residue, len(params.input.dir))

        if not os.path.exists(os.path.join(params.output.out_dir, residue, interpolate_csv)):

            # Iterate through the datasets
            for dataset_label, dataset_results in all_results.iteritems():

                current_dataset_results = dataset_results.loc[(dataset_results.index == residue)]
                sorted_angles, sorted_map_values = normalise_and_sort_ringer_results(current_dataset_results,
                                                                                     params=params)
                interpolated_angles, interpolated_map_values = linear_interpolate_ringer_results(sorted_angles,
                                                                        sorted_map_values,
                                                                        angle_sampling = params.settings.angle_sampling)

                # Store these in a list
                residue_data_list.append((interpolated_angles,
                                          interpolated_map_values))
                # If it doesn't exist: Create dataframe to store results from
                # one residue, across multiple datasets
                if not 'single_residue_multiple_datasets' in locals():
                    single_residue_multiple_datasets = pandas.DataFrame(
                        columns=
                        interpolated_angles)

                # Populate dataframe with results from one residue,
                # across multiple datasets
                single_residue_multiple_datasets.loc['{}'.format(dataset_label)] = interpolated_map_values


                # Print results for all of the datasets for this residue in the same graph
            multiple_line_plot_ringer(residue_data_list,
                                      title=residue,
                                      filename='all-{}-{}-dataset.png'.format(
                                          residue, len(residue_data_list)),
                                      out_dir=os.path.join(params.output.out_dir, residue))

            print(single_residue_multiple_datasets)

            # Output CSV from one resiude, multiple datasets
            pandas.DataFrame.to_csv(single_residue_multiple_datasets,
                                    os.path.join(params.output.out_dir, residue,
                                                 interpolate_csv))
        else:
            logger.info(
                '{}: Interpolated CSVs already generated,'.format(residue))

        # Generate correlation CSV
        if not os.path.exists(os.path.join(params.output.out_dir,
                                           residue, correlation_csv)):

            correlation_single_residue(input_csv=os.path.join(params.output.out_dir,
                                                              residue,
                                                              interpolate_csv),
                                       residue=residue,
                                       output_dir=os.path.join(params.output.out_dir, residue),
                                       out_filename=correlation_csv,
                                       params=params)
        else:
            logger.info('{}: Correlation CSV already generated,'.format(residue))

    ##########################################################################
    # Generate Average Ringer Plots
    ##########################################################################
    average_type = "Median"
    if not os.path.exists(os.path.join(params.output.out_dir, '{}_ringer_results'.format(average_type),
                                       '{}_ringer_{}_datasets.csv'.format(average_type,
                                                                          len(params.input.dir)))):
        average_ringer_plots(base_csv=interpolate_base_csv, ref_set=ref_set,
                             out_dir=params.output.out_dir, params=params,
                             average_type=average_type)
    else:
        logger.info('Average Ringer plots already generated')

    # Blue line plots
    # average_ringer_plots(base_csv=interpolate_base_csv, ref_set=ref_set,
    #                      out_dir=out_dir, params=params,
    #                      average_type=average_type,
    #                      bold_blue_map=None)


    ###########################################################################
    # Clustering for correlation
    ###########################################################################
    correlation_csv_end = '_from {} datasets-correlation-ringer.csv'.format(len(params.input.dir))
    min_corr, max_corr = find_pairwise_range(correlation_csv_end, ref_set, params.output.out_dir)

    hier_agg_cluster(correlation_csv_end, pairwise_type='correlation',
                     ref_set=ref_set, out_dir=params.output.out_dir, params=params)

    ##########################################################################
    # Generate heatmaps from clustering:Correlation
    ###########################################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format('correlation', '', '')

    if not os.path.exists(os.path.join(params.output.out_dir,
                                       "Adj_cluster_weight_heatmap_{}_{}_{}.png".format('correlation', '', ''))):
        cluster_heatmap(os.path.join(params.output.out_dir, clusters_weight_filename),
                        params.output.out_dir, pairwise_type='correlation', fit_type='', subset='')
    ###########################################################################
    # Generating histogram to show the location of the maximal points 
    # of peak in ringer plot. Shows three rotamer bins [60, 180, 300] 
    ###########################################################################

    if not os.path.exists(os.path.join(params.output.out_dir, 'Modal_peak_location.png')):

        max_peak_angle = []
        # Generate maximal values from interpolated map values
        for residue, data in ref_set.iterrows():
            interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue, len(params.input.dir), map_type,
                                                                       angle_type)

            interpolated_results = pandas.read_csv(os.path.join(params.output.out_dir, residue,
                                                                interpolate_csv), index_col=0)

            assert len(interpolated_results) == dataset_counter,\
                ('Input CSV data is length {} for {} datasets.'
                'Lengths should match'.format(len(interpolated_results), dataset_counter))

            angle_with_max_map = (interpolated_results.idxmax(axis=1).values).astype(numpy.float)
            max_peak_angle.append(angle_with_max_map)

            # Plot histogram
        peak_angle_histogram(max_peak_angle, params.output.out_dir)
    else:
        logger.info('Histogram of peak angles exists')

    # TODO Seperat into seperat files/ and set up phil file to run with and without these features

   
    ##########################################################################
    # Curve Fitting routine for all datasets
    ##########################################################################  
    fit_type = 'three_gaussian_offset'
    subset='Amplitudes'
    pairwise_type = 'euclidean distance'
    mean_bound = None
    fit_all_datasets(params.output.out_dir,ref_set,map_type,angle_type,
                    params,pairwise_type,fit_type, mean_bound = mean_bound)

    fit_base_filename = '_from_{}_datasets_{}.csv'.format(len(params.input.dir)
                                                          ,fit_type)
    ##########################################################################
    # Generate RMSD between fit and data. Plot histogram of all RMSD values,
    # Store RMSD values in single CSV for the fit type
    #########################################################################
    generate_RMSD(params.output.out_dir, ref_set, map_type, angle_type, fit_type,
                  fit_base_filename, params=params)

    ###########################################################################
    # Calculate Euclidean distance
    ###########################################################################
    euclidean_base_csv = calculate_euclidean_distance(params.output.out_dir, ref_set,
                         params, fit_type, subset=subset)
    ###########################################################################
    # Heirichal Clustering, Average linakge, for pairwise euclidean distances
    # for each residue
    ########################################################################### 
    hier_agg_cluster(euclidean_base_csv, pairwise_type, ref_set, params.output.out_dir,
                     params, fit_type, subset, incons_threshold = 3, depth = 10)
    hier_agg_cluster(euclidean_base_csv, pairwise_type, ref_set, params.output.out_dir,
                     params, fit_type, subset, incons_threshold = 4, depth = 10)
    hier_agg_cluster(euclidean_base_csv, pairwise_type, ref_set, params.output.out_dir,
                     params, fit_type, subset, incons_threshold = 3, depth = 5)
    hier_agg_cluster(euclidean_base_csv, pairwise_type, ref_set, params.output.out_dir,
                     params, fit_type, subset, incons_threshold = 4, depth = 5)

    ##########################################################################
    # Generate Heatmaps from clustering: Amplitudes
    ##########################################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                pairwise_type, fit_type, subset)
    if not os.path.exists(os.path.join(params.output.out_dir,
        "Adj_cluster_weight_heatmap_{}_{}_{}.png".format(pairwise_type,fit_type,subset))):
        cluster_heatmap(os.path.join(params.output.out_dir,clusters_weight_filename),
                    params.output.out_dir,pairwise_type,fit_type = fit_type, subset = subset)

    ##########################################################################
    # Clustering with Amplitudes & Means   
    ######################################################################### 
    subset='Amplitudes_Means'
    euclidean_base_csv = calculate_euclidean_distance(params.output.out_dir,ref_set,
                         params,fit_type,subset=subset)

    hier_agg_cluster(euclidean_base_csv,pairwise_type,ref_set,params.output.out_dir, params,
                    fit_type,subset)
     
    ##########################################################################
    # Generate Heatmaps from clustering: Amplitudes and Means
    ##########################################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                pairwise_type, fit_type, subset)
    if not os.path.exists(os.path.join(params.output.out_dir,
        "Adj_cluster_weight_heatmap_{}_{}_{}.png".format(pairwise_type,fit_type,subset))):
        cluster_heatmap(os.path.join(params.output.out_dir,clusters_weight_filename),
                        params.output.out_dir,pairwise_type,fit_type = fit_type,
                        subset = subset)
    #########################################################################
    # Explore cluster weights
    #########################################################################
    subset ='Amplitudes' 
    Adj_clusters_weight = pandas.read_csv(os.path.join(params.output.out_dir,
                                      'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                      pairwise_type, fit_type, subset))
                                      ,header=0, index_col =0)
    Adj_clusters_weight_corr = pandas.read_csv(os.path.join(params.output.out_dir,
                                      'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                      'correlation', '', ''))
                                      ,header=0, index_col =0)
    subset='Amplitudes_Means'
    Adj_clusters_weight_means = pandas.read_csv(os.path.join(params.output.out_dir,
                                      'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                      pairwise_type, fit_type, subset))
                                      ,header=0, index_col =0)
    ##########################################################################
    # Read RMSD data
    #########################################################################
    RMSD_filename = 'RMSD, with {}.csv'.format(fit_type)

    all_RMSD = pandas.read_csv(os.path.join(params.output.out_dir,
                    RMSD_filename),index_col = 0,
                    header = 0)
    
    ###########################################################################
    # Generate Scores
    ##########################################################################

    dataset_score = Adj_clusters_weight.sum()
    residue_score = (Adj_clusters_weight**2).sum(axis = 1)
    dataset_corr_score = Adj_clusters_weight_corr.sum()
    residue_corr_score = (Adj_clusters_weight_corr**2).sum(axis=1)
    dataset_means_score = Adj_clusters_weight_means.sum()
    residue_means_score = (Adj_clusters_weight_means**2 ).sum(axis =1)

    score_matrix= pandas.DataFrame(index=residue_score.index, columns=dataset_score.index)

    for dataset in dataset_score.index:
        for residue in residue_score.index:
            score_matrix.loc[residue][dataset]=(dataset_score.loc[dataset])*(residue_score.loc[residue])
  
    score_matrix.to_csv('score_matrix.csv')

    dataset_resolution = pandas.read_csv(resolution_csv_path, header =0, index_col =0)
        
    ###########################################################################
    # Reading in subset of bound ligands 
    ###########################################################################
    bound_ligands=pandas.read_csv('bound_ligands.csv')

    ##########################################################################
    # Plotting 
    ###########################################################################
    plot_correlation_vs_fitting(params.output.out_dir,dataset_corr_score,dataset_score)
    plot_resloution_vs_dataset_score(params.output.out_dir,dataset_score,dataset_means_score,
                                     dataset_resolution)
    plot_RMSD_vs_dataset_score(params.output.out_dir,all_RMSD,dataset_score,bound_ligands)
    plot_RMSD_vs_resolution(params.output.out_dir,all_RMSD,dataset_resolution)
    plot_RMSD_vs_residue_score(params.output.out_dir,all_RMSD,residue_score,
                               residue_means_score)
    ###########################################################################
    # Peak Finding Routine
    ###########################################################################
    find_peaks(params.output.out_dir, ref_set,params, map_type, angle_type)



# For allowing command manager
if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
