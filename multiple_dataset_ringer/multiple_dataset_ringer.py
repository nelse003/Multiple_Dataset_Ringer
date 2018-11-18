import glob
import libtbx.phil
import logging
import logging.config
import matplotlib
import numpy
import os
import pandas as pd
import sys
from string import ascii_letters

# Correlation functions
from correlation.correlation import correlation_single_residue
from fitting.fitting import calculate_euclidean_distance
# Fitting functions
from fitting.fitting import fit_all_datasets
from fitting.fitting import generate_RMSD
# Hierarichal clustering functions
from cluster.hier import hier_agg_cluster
from cluster.hier import find_pairwise_range
from plotting.plots import average_ringer_plots
from plotting.plots import cluster_heatmap
# Plotting functions
from plotting.plots import line_plot_ringer
from plotting.plots import multiple_line_plot_ringer
from plotting.plots import pairwise_heatmap
from plotting.plots import peak_angle_histogram
from plotting.plots import plot_RMSD_vs_dataset_score
from plotting.plots import plot_RMSD_vs_residue_score
from plotting.plots import plot_RMSD_vs_resolution
from plotting.plots import plot_correlation_vs_fitting
from plotting.plots import plot_resloution_vs_dataset_score
# Sorting & Interpolating angles from ringer output
from process.interpolate import interpolate_all_ringer_results
# Ringer Processing with absolute electron density scaling (F000)
from process.process import process_with_ringer
from process.process import process_all_with_ringer

matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot

pyplot.style.use('ggplot')

# TODO Add back in and test peak finding
#from peak_count import find_peaks

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
    column_labels = "FWT,PHWT"
        .type = str
        .multiple = False
        .help = mtz column labels of input files
}
output {
    log = "ringer.log"
        .type = str
        .multiple = False
    out_dir = "output"
        .type = str
        .multiple = False
    tmp_dir = "tmp"
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
        
    map_type = '2mFo-DFc'
        .type = str
        .help = Name of electron density map to be analysed with ringer
        
    angle_type = 'chi1'
        .type = str
        .help = Chi angle to be analysed across multiple datasets
        
    sample_only_ref_pdb = None
        .type = str
        .help = Run all ringer analyses against a single pdb structure
    qsub = False
        .type = bool
        .help = flag to run (ringer) with qsub
}
""")

###############################################################################
# Logging
###############################################################################
logger = logging.getLogger(__name__)

logging.config.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s - %(levelname)s - %(message)s'
        },
    },
    'handlers': {
        'default': {
            'level': 'INFO',
            'class': 'logging.StreamHandler'
        },
        'file': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': 'ringer.log'
        },
    },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'INFO',
            'propagate': True
        }
    }
})


def run(params):

    """ Run Ringer over multiple datasets.


    Notes
    ------------------
    Needs to be run with a input (params.input.dir) directory organised with
    pdb files in /dir/*pdb_style*
    
    Process ringer uses 2FOFCWT, PHI2FOFCWT columns. 
    Handle other columns by moving to input check.

    Example
    ------------------
    Command Line
    
    ccp4-python multiple_dataset_ringer.py \
    /hdlocal/home/enelson/DCP2B_ringer_test/*/ \
    input.pdb_style="refine.pdb" \
    input.mtz_style="refine.mtz" \
    output.out_dir="/hdlocal/home/enelson/DCP2B_ringer_test/output_test"

    """

    logger.info("Produce a dictionary of Dataframes, "
                "each containing ringer results for a dataset")
    all_results = process_all_with_ringer(params)

    # convert to df works, but is large as csv ~300mb for 1600 datasets.
    # may not be more use than seperate csvs.

    #all_results_df = pd.concat(all_results, keys=all_results.keys())

    datasets = all_results.keys()

    logger.info("Pull out the first ringer results set as a reference")
    ref_set = all_results.values()[0]

    logger.info("Choose a map_type: {}\n"
                "angle_type: {} \n"
                " by reducing reference set".format(params.settings.map_type,
                                                    params.settings.angle_type))

    ref_set = ref_set.loc[(ref_set[1] == params.settings.map_type)]
    ref_set = ref_set.loc[(ref_set[2] == params.settings.angle_type)]

    logger.info("Interpolating all ringer results to be on same angle range")
    interpolate_all_ringer_results(ref_set, all_results, params)

    interpolate_base_csv = \
        '_{}_Datasets_{}_{}-ringer.csv'.format(len(datasets),
                                               params.settings.map_type,
                                               params.settings.angle_type)

    logger.info("Producing ringer line plots for each residue: showing all datasets")
    for residue, data in ref_set.iterrows():

        interpolate_csv = residue + interpolate_base_csv

        single_residue_multiple_datasets = pd.read_csv(
            os.path.join(params.output.out_dir,
                         residue,
                         interpolate_csv),
            index_col=0)

        if not os.path.exists(os.path.join(params.output.out_dir,
                                           residue,
                                           'all-{}-{}-dataset.png'.format(
                                               residue, len(ref_set)))):

            multiple_line_plot_ringer(results_df=single_residue_multiple_datasets,
                                      title=residue,
                                      filename='all-{}-{}-dataset.png'.format(
                                          residue, len(ref_set)),
                                      out_dir=os.path.join(params.output.out_dir, residue))
        else:
            logger.info("{} :Ringer plot for already generated".format(residue))
    exit()
    # Generate correlations between datasets for each residue
    for residue, data in ref_set.iterrows():

        interpolate_csv = residue + interpolate_base_csv

        # Output filename for correlation data
        correlation_csv = '{}_from {} datasets-correlation-ringer.csv'.format(residue, len(datasets))

        # Generate correlation CSV
        if not os.path.exists(os.path.join(params.output.out_dir,
                                           residue, correlation_csv)):

            correlation_single_residue(input_csv=os.path.join(params.output.out_dir,
                                                              residue,
                                                              interpolate_csv),
                                       output_dir=os.path.join(params.output.out_dir, residue),
                                       out_filename=correlation_csv)
        else:
            logger.info('{}: Correlation CSV already generated,'.format(residue))


    ##########################################################################
    # Generate Average Ringer Plots
    ##########################################################################
    average_type = "Median"
    if not os.path.exists(os.path.join(params.output.out_dir,
                                       '{}_ringer_results'.format(average_type),
                                       '{}_ringer_{}_datasets.csv'.format(average_type,
                                                                          len(datasets)))):

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


    #################################
    # Clustering for correlation
    #################################
    correlation_csv_end = '_from {} datasets-correlation-ringer.csv'.format(len(datasets))
    min_corr, max_corr = find_pairwise_range(correlation_csv_end, ref_set, params.output.out_dir)

    hier_agg_cluster(correlation_csv_end,
                     datasets=datasets,
                     pairwise_type='correlation',
                     ref_set=ref_set,
                     out_dir=params.output.out_dir,
                     params=params)

    #################################################
    # Generate heatmaps from clustering:Correlation
    #################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format('correlation', '', '')

    if not os.path.exists(os.path.join(params.output.out_dir,
                                       "Adj_cluster_weight_heatmap_{}_{}_{}.png".format('correlation', '', ''))):
        cluster_heatmap(
            os.path.join(params.output.out_dir, clusters_weight_filename),
            params.output.out_dir,
            pairwise_type='correlation',
            fit_type='',
            subset='')

    ###################################################################
    # Generating histogram to show the location of the maximal points 
    # of peak in ringer plot. Shows three rotamer bins [60, 180, 300] 
    ###################################################################

    if not os.path.exists(os.path.join(params.output.out_dir,
                                       'Modal_peak_location.png')):

        max_peak_angle = []
        # Generate maximal values from interpolated map values
        for residue, data in ref_set.iterrows():
            interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue,
                                                                       len(datasets),
                                                                       params.settings.map_type,
                                                                       params.settings.angle_type)

            interpolated_results = pd.read_csv(
                os.path.join(params.output.out_dir, residue, interpolate_csv),
                index_col=0)

            angle_with_max_map = (interpolated_results.idxmax(axis=1).values).astype(numpy.float)
            max_peak_angle.append(angle_with_max_map)

        # Plot histogram
        peak_angle_histogram(max_peak_angle, params.output.out_dir)
    else:
        logger.info('Histogram of peak angles exists')

    # TODO Separate into seperate files/ and set up phil file to run with and without these features

    ##########################################################################
    # Curve Fitting routine for all datasets
    ##########################################################################  
    fit_type = 'three_gaussian_offset'
    subset='Amplitudes'
    pairwise_type = 'euclidean distance'
    mean_bound = None

    fit_all_datasets(params.output.out_dir,
                     ref_set,
                     params.settings.map_type,
                     params.settings.angle_type,
                     params,
                     pairwise_type,
                     fit_type,
                     mean_bound=mean_bound,
                     datasets=datasets)

    fit_base_filename = '_from_{}_datasets_{}.csv'.format(len(datasets),
                                                          fit_type)

    ###################################################################
    # Generate RMSD between fit and data.
    # Plot histogram of all RMSD values.
    # Store RMSD values in single CSV for the fit type
    ###################################################################
    generate_RMSD(params.output.out_dir,
                  ref_set,
                  params.settings.map_type,
                  params.settings.angle_type,
                  fit_type,
                  fit_base_filename,
                  datasets=all_results.keys())

    ###################################################################
    # Calculate Euclidean distance
    ###################################################################
    euclidean_base_csv = calculate_euclidean_distance(params.output.out_dir,
                                                      ref_set,
                                                      params,
                                                      fit_type,
                                                      datasets=datasets,
                                                      subset=subset)

    ###################################################################
    # Heirichal Clustering, Average linakge, for pairwise euclidean
    # distances for each residue
    ###################################################################

    hier_agg_cluster(euclidean_base_csv,
                     pairwise_type,
                     ref_set,
                     params.output.out_dir,
                     params,
                     datasets=datasets,
                     fit_type=fit_type,
                     subset=subset,
                     incons_threshold=3,
                     depth=10)

    # Note was also previously tested with
    # incons_threshold=4, depth=10
    # incons_threshold=3, depth=5
    # incons_threshold=5, depth=5


    ##########################################################################
    # Generate Heatmaps from clustering: Amplitudes
    ##########################################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                pairwise_type, fit_type, subset)

    if not os.path.exists(os.path.join(
            params.output.out_dir,
            "Adj_cluster_weight_heatmap_{}_{}_{}.png".format(
                pairwise_type,
                fit_type,
                subset))):

        cluster_heatmap(os.path.join(params.output.out_dir,
                                     clusters_weight_filename),
                        params.output.out_dir,
                        pairwise_type,
                        fit_type = fit_type,
                        subset = subset)

    ##########################################################################
    # Clustering with Amplitudes & Means   
    ######################################################################### 
    subset='Amplitudes_Means'
    euclidean_base_csv = calculate_euclidean_distance(params.output.out_dir,
                                                      ref_set,
                                                      params,
                                                      fit_type,
                                                      datasets=datasets,
                                                      subset=subset)

    hier_agg_cluster(euclidean_base_csv,
                     pairwise_type,
                     ref_set,
                     params.output.out_dir,
                     params,
                     datasets=datasets,
                     fit_type=fit_type,
                     subset=subset)
     
    ##########################################################################
    # Generate Heatmaps from clustering: Amplitudes and Means
    ##########################################################################
    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format(
                                pairwise_type, fit_type, subset)

    if not os.path.exists(os.path.join(
            params.output.out_dir,
            "Adj_cluster_weight_heatmap_{}_{}_{}.png".format(pairwise_type,
                                                             fit_type,subset))):

        cluster_heatmap(os.path.join(params.output.out_dir,
                                     clusters_weight_filename),
                        params.output.out_dir,
                        pairwise_type,
                        fit_type=fit_type,
                        subset=subset)

    #########################################################################
    # Explore cluster weights
    #########################################################################
    subset ='Amplitudes' 
    Adj_clusters_weight = pd.read_csv(
        os.path.join(params.output.out_dir,
                     'Adj_clusters_weight_{}_{}_{}.csv'.format(
                         pairwise_type, fit_type, subset)),
        header=0,
        index_col=0)

    Adj_clusters_weight_corr = pd.read_csv(
        os.path.join(params.output.out_dir,
                     'Adj_clusters_weight_{}_{}_{}.csv'.format('correlation', '', '')),
        header=0,
        index_col=0)

    subset='Amplitudes_Means'
    Adj_clusters_weight_means = pd.read_csv(
        os.path.join(params.output.out_dir,
                     'Adj_clusters_weight_{}_{}_{}.csv'.format(
                         pairwise_type, fit_type, subset)),
        header=0,
        index_col=0)

    ##########################################################################
    # Read RMSD data
    #########################################################################
    RMSD_filename = 'RMSD, with {}.csv'.format(fit_type)

    all_RMSD = pd.read_csv(os.path.join(params.output.out_dir,
                                            RMSD_filename),
                               index_col=0,
                               header=0)
    
    ###########################################################################
    # Generate Scores
    ##########################################################################
    dataset_score = Adj_clusters_weight.sum()
    residue_score = (Adj_clusters_weight**2).sum(axis=1)
    dataset_corr_score = Adj_clusters_weight_corr.sum()
    residue_corr_score = (Adj_clusters_weight_corr**2).sum(axis=1)
    dataset_means_score = Adj_clusters_weight_means.sum()
    residue_means_score = (Adj_clusters_weight_means**2).sum(axis=1)

    score_matrix = pd.DataFrame(index=residue_score.index,
                                    columns=dataset_score.index)

    for dataset in dataset_score.index:
        for residue in residue_score.index:
            score_matrix.loc[residue][dataset] = \
                (dataset_score.loc[dataset]) * (residue_score.loc[residue])
  
    score_matrix.to_csv('score_matrix.csv')

    dataset_resolution = pd.read_csv(resolution_csv_path,
                                         header=0,
                                         index_col=0)
        
    ######################################
    # Reading in subset of bound ligands 
    ######################################
    # bound_ligands=pd.read_csv('bound_ligands.csv')

    ###############
    # Plotting 
    ###############
    plot_correlation_vs_fitting(params.output.out_dir, dataset_corr_score,
                                dataset_score)

    # TODO rewrite resolution based code
    # plot_resloution_vs_dataset_score(params.output.out_dir, dataset_score,
    #                                  dataset_means_score, dataset_resolution)

    # TODO rewrite boun ligands lsit generation
    # plot_RMSD_vs_dataset_score(params.output.out_dir, all_RMSD,
    #                            dataset_score, bound_ligands)

    # TODO rewrite resolution based code
    # plot_RMSD_vs_resolution(params.output.out_dir, all_RMSD,
    #                         dataset_resolution)

    plot_RMSD_vs_residue_score(params.output.out_dir, all_RMSD, residue_score,
                               residue_means_score)

    # TODO Peak fitting code: Assess, and install peak extra if needed
    #######################
    # Peak Finding Routine
    #######################
    # find_peaks(params.output.out_dir, ref_set,params,
    #            params.settings.map_type, params.settings.angle_type)


# For allowing command manager
if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)
