import os, sys, copy, glob

import pandas,numpy,time

# Hierarchy
from scipy.cluster.hierarchy import linkage, fcluster

import logging

# Imported functions
from plotting_ringer import plot_dendrogram, number_clusters_histogram

###############################################################
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


###################################################################
# Hierarichial agglomerative clustering used for selecting out
# datasets that are consitently different (Inconsitency condition)
###################################################################

def hier_agg_cluster(base_input_csv, pairwise_type, ref_set, out_dir, params,
                     datasets, fit_type = '', subset ='',
                     incons_threshold = 3, depth = 10):

    """ Hierarichial agglomerative clustering"""

    # Generate range over which metric ranges
    min_range, max_range = find_pairwise_range(base_input_csv, ref_set, out_dir)

    num_cluster_all=[]
    clusters_weight = pandas.DataFrame(index = ref_set.index.values,
                                       columns = datasets)

    clusters_weight_filename = 'Adj_clusters_weight_{}_{}_{}.csv'.format(
        pairwise_type, fit_type ,subset)

    for residue, data in ref_set.iterrows():
        input_csv = '{}'.format(residue) + base_input_csv

        # Generate Linkage matrix and dendrogram 
        dendrogram_filename = '{}_{}_{}_{}_dendrogram.png'.format(residue,
                                                                  pairwise_type,
                                                                  fit_type,
                                                                  subset)

        cluster_number_hist = 'number_cluster_{}_{}_{}_{}_{}.png'.format(
            pairwise_type, fit_type, subset, incons_threshold, depth)

        linkage_matrix, dataset_labels = generate_linkage_matrix(
            os.path.join(out_dir, residue, input_csv))

        if not os.path.exists(os.path.join(out_dir,
                                           residue,
                                           dendrogram_filename)):
            start=time.time()
            plot_dendrogram(linkage_matrix,
                            out_dir,
                            residue,
                            pairwise_type,
                            dendrogram_filename,
                            dataset_labels=dataset_labels)
            end= time.time()
            duration = end-start
            logger.info('{}: Dendrogram ({}) generated in {} seconds.'.format(
                        residue, pairwise_type, duration))
        else:
            logger.info('{}: Dendrogram already generated for {} '
                        'with {}'.format(residue, pairwise_type, fit_type))

        #######################################################################
        # Clustering using inconsitency matrix.
        # Generates a weight for each residue dataset pair, 
        # given by 1/(cluster size)
        #####################################################################
        if not os.path.exists(os.path.join(out_dir,clusters_weight_filename)):

            cluster_groups = fcluster(linkage_matrix, incons_threshold, depth = depth)
            clusters = pandas.DataFrame(data = cluster_groups)
            clusters.columns = ['Cluster_number']
            number_of_clusters = clusters.values.max()
            # group into clusters
            grouped = clusters.groupby('Cluster_number')

            # look at cluster with less than 20 members
            #small_clusters = grouped.filter(lambda x: len(x) < 20) 
            num_cluster_all.append(number_of_clusters)
            for i in range(1,number_of_clusters+1):
                group_size = len(grouped.groups[i])
                weight = 1/group_size
                for j in range(0,group_size):
                    clusters_weight.loc[residue:,grouped.groups[i][j]] = weight/number_of_clusters  
        else:
            logger.info('Cluster weights ({}) already generated'.format(pairwise_type))
        # Generate Heatmap

        heatmap_filename = '{}_{}_{}_{}_heatmap.png'.format(residue,pairwise_type,
                                                            fit_type,subset)
        if params.settings.gen_heatmap == True:
            if not os.path.exists(os.path.join(out_dir,residue,heatmap_filename)):
                start = time.time()
                pairwise_heatmap(os.path.join(out_dir,residue,input_csv),
                                residue=residue,out_dir = out_dir,
                                pairwise_type = pairwise_type,
                                plot_filename = heatmap_filename, 
                                fit_type = fit_type,
                                max_scale = max_range, min_scale = min_range,
                                params = params)
                end = time.time()
                duration = end-start
                logger.info('{}: Heatmap ({} with {}) for {} datasets generated in {} seconds.'.format(
                            residue,pairwise_type,fit_type,len(params.input.dir)
                            ,duration))
            else:
                logger.info('{}: Heatmap already generated for '
                            '{} with {} and {}'.format(residue,
                                                       pairwise_type,
                                                       fit_type, subset))
        else:
            logger.info('Skip Heatmap')

    # Send cluster weights to file
    if not os.path.exists(os.path.join(out_dir,clusters_weight_filename)):
        clusters_weight.to_csv(os.path.join(out_dir,clusters_weight_filename))

    # Histogram
    if not os.path.exists(os.path.join(out_dir,cluster_number_hist)):
        if not num_cluster_all:
            logger.info("Number of cluster not generated")
        else:
            number_clusters_histogram(num_cluster_all,cluster_number_hist,
                                  out_dir)

def generate_linkage_matrix(pairwise_csv):
    """ Generate Linkage matrix """

    # Load csv into pandas DataFrame
    data = pandas.read_csv(pairwise_csv, index_col=0)
    dataset_labels = data.columns.values
    # Generate linkage matrix 
    linkage_matrix = linkage(data.values, "single")

    return linkage_matrix, dataset_labels

def find_pairwise_range(pairwise_csv_end, ref_set, out_dir):
    """ Find minimum value across all pairwise correlations"""

    pairwise_min=[]
    pairwise_max=[]

    for residue, data in ref_set.iterrows():

        pairwise_csv='{}'.format(residue) + pairwise_csv_end
        # Load csv into pandas DataFrame
        data = pandas.read_csv(os.path.join(out_dir, residue, pairwise_csv), index_col=0)
        pairwise_min.append(min(data.min()))
        pairwise_max.append(max(data.max()))

    return min(pairwise_min), max(pairwise_max)

