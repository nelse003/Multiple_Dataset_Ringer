#!/usr/bin/env pandda.python

###############################################################################
# Packages
##############################################################################
# System tasks
import os, sys, copy, glob
# Numpy
import numpy
# Pandas Dataframes
import pandas
# Dendrogram
from scipy.cluster.hierarchy import dendrogram
# Spearman rank coefficent for plot label
from scipy.stats import spearmanr

####################################
# Plotting prettily using matplotlib
####################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
matplotlib.rcParams['savefig.facecolor'] = 'white'
matplotlib.rcParams['savefig.dpi'] = 80
matplotlib.rcParams['axes.facecolor'] = 'white'
matplotlib.rcParams['xtick.labelsize'] = 18
matplotlib.rcParams['ytick.labelsize'] = 18
import matplotlib.lines as mlines
##################################

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


###############################################################################
# Functions
##############################################################################

# Set some parameters for plots: Only left and top axes
def myfig(**kwargs):
    fig = pyplot.figure(**kwargs)
    ax = fig.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', 
            color='0.9')
    ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', 
            color='0.9')
    ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', 
            color='0.9')
    ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', 
            color='0.9')
    return fig, ax

def line_plot_ringer(sorted_angles,sorted_map_values,title,filename,out_dir):
    """Plot single ringer graph """
    fig = myfig()
    pyplot.plot(sorted_angles, sorted_map_values)
    pyplot.title(title)
    pyplot.xlabel('Angle')
    pyplot.tight_layout()
    pyplot.savefig(os.path.join(out_dir,filename))
    pyplot.close(fig)

def multiple_line_plot_ringer(all_data_list,title, filename, out_dir):
    """Plot multiple ringer plots  """
    fig = myfig()
    pyplot.title(title)

    for i in range(0,len(all_data_list)):
        sorted_angles=all_data_list[i][0]
        sorted_map_values=all_data_list[i][1]
        pyplot.plot(sorted_angles,sorted_map_values)

    pyplot.xlabel('Angle')
    pyplot.tight_layout()
    pyplot.savefig(os.path.join(out_dir,filename))
    pyplot.close()

def multiple_line_plot_bold(bold_angles, bold_map_values,all_data_list, title,
                            filename, out_dir, average_type, 
                            bold_blue_map_values = None, bold_blue_map = None):
    """ Plot multiple ringer plots, with a bold plot ontop"""
    fig = myfig(figsize=(10,8))
    pyplot.title(title, fontsize = 24)
    pyplot.xlabel('Angle',fontsize = 20)
    pyplot.ylabel('Map Value', fontsize = 20)
    pyplot.xlim(0,360)

    # Legend
    black_line = mlines.Line2D([],[], color ='k', ls ='--',
                                label = '{} Ringer Plot'.format(average_type))
    grey_line = mlines.Line2D([],[], color ='k', alpha = 0.05, label = 'Ringer plot for each datset')
    if bold_blue_map_values is not None:
        blue_line =mlines.Line2D([],[], color ='b',ls = '-.', 
                                 label = '{} Ringer plot'.format(bold_blue_map))
        legend = pyplot.legend(handles =[black_line,grey_line,blue_line], 
                               loc = 'best', fontsize =20)
    else:
        legend = pyplot.legend(handles =[black_line,grey_line], 
                               loc= 'best', fontsize = 20)
    legend.get_frame().set_facecolor('w')
     # Background Plots     
    for i in range(0,len(all_data_list)):
        sorted_angles=all_data_list[i][0]
        sorted_map_values=all_data_list[i][1]
        pyplot.plot(sorted_angles,sorted_map_values, color= 'k', alpha = 0.05)
    # Foreground Plot   
    pyplot.plot(bold_angles, bold_map_values, color = 'k', linewidth = 2, ls ='--')
    if bold_blue_map_values is not None:
        pyplot.plot(sorted_angles, bold_blue_map_values, color = 'b', 
                    linewidth = 2, ls = '-.')

    #pyplot.gca().patch.set_facecolor('0.95')
    pyplot.tight_layout()
    pyplot.savefig(os.path.join(out_dir,filename), dpi = 300, format = 'png')
    pyplot.close()


def average_ringer_plots(base_csv,ref_set,out_dir,params,average_type= 'Mean',
                         bold_blue_map = None):
    """Calculate average ringer plot for each residue across all datasets"""
    all_average_results = pandas.DataFrame() #index = ref_set.index.values

    for residue, data in ref_set.iterrows():
        interpolate_csv = residue + base_csv
        results = pandas.read_csv(os.path.join(out_dir,residue,interpolate_csv),
                                  header = 0, index_col = 0)
        if average_type == 'Mean':
            average_results = results.mean(axis = 0)
        elif average_type == 'Median':
            average_results = results.median(axis =0)
        else:
            logger.warning("Average type not Mean or Median")

        average_results.name = residue
        average_results.index = average_results.index.values.astype(int)

        output_dir = os.path.join(out_dir,'{}_ringer_results'.format(average_type))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        angles = average_results.index.values.astype(int)
        plot_filename = '{}_ringer_{}_with_{}_datasets.png'.format(average_type,
                                                                   residue,
                                                            len(params.input.dir))
        all_average_results = all_average_results.append(average_results)

        if bold_blue_map is not None:
            bold_blue_map_values = results.loc[bold_blue_map]
            plot_filename = 'blue_{}_ringer_{}_with_{}_datasets.png'.format(average_type,
                                                                   residue,
                                                            len(params.input.dir))
        else:
            bold_blue_map_values = None

        # Multiple lines plotted in background

        # Generate list of map values
        all_data_list = []
        for dataset in results.index:
            map_values = results.loc[dataset].values
            all_data_list.append((angles,map_values))

        multiple_line_plot_bold(bold_angles = angles,
                                bold_map_values = average_results.values,
                                all_data_list = all_data_list,
                                title = '{} ringer {} with {} datasets'.format(
                                        average_type,residue,len(params.input.dir)),
                                filename = plot_filename,
                                out_dir = output_dir,
                                average_type = average_type, 
                                bold_blue_map_values = bold_blue_map_values,
                                bold_blue_map = bold_blue_map)

    all_average_results.to_csv(os.path.join(output_dir,
                            '{}_ringer_{}_datasets.csv'.format(average_type,len(params.input.dir))))


def augmented_dendrogram(*args, **kwargs):
    """Improved dendrogram plotting for Hierarchal clustering"""
    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        for i, d in zip(ddata['icoord'], ddata['dcoord']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            pyplot.plot(x, y, 'ro')
            pyplot.annotate("%.3g" % y, (x, y), xytext=(0, -8),
                         textcoords='offset points',
                         va='top', ha='center')
    return ddata


def plot_dendrogram(linkage_matrix,out_dir,residue,pairwise_type,plot_filename,
                    dataset_labels):
    # Plotting dendorgram
    fig = myfig(figsize=(10,10))
    pyplot.clf()
    ddata = augmented_dendrogram(linkage_matrix,color_threshold=1,p=15,
                                truncate_mode='lastp',labels=dataset_labels)
    pyplot.gcf().subplots_adjust(bottom=0.25)
    pyplot.xticks(rotation=90)
    pyplot.title("{} {} Dendrogram".format(residue,pairwise_type))
    pyplot.ylabel("Cophenetic distance")
    pyplot.savefig(os.path.join(out_dir,residue,plot_filename))
    pyplot.close()

def plot_correlation_vs_fitting(out_dir,dataset_corr_score,dataset_score):

    fig = myfig() 
    fig, ax = pyplot.subplots(1)
    ax.scatter(dataset_corr_score,dataset_score)
    pyplot.xlabel('Dataset score: correlation')
    pyplot.ylabel('Dataset score: Amplitude only fit')
    ax.set_xlim(left =0)
    ax.set_ylim(bottom =0)
    spearman_dataset, p_val = spearmanr(dataset_corr_score,dataset_score)
    pyplot.title('Spearman corr = {}'.format(spearman_dataset))
    pyplot.savefig(os.path.join(out_dir,
                   'dataset_score_correlation_fitting_vs_correlation.png'))
    pyplot.close()

def plot_resloution_vs_dataset_score(out_dir,dataset_score,dataset_means_score,
                                     dataset_resolution):
    fig, (ax1,ax2) = pyplot.subplots(2, sharex=True, sharey = True)
    ax1.scatter(dataset_resolution,dataset_score,label='Amplitude Fit')
    ax2.scatter(dataset_resolution,dataset_means_score,
                   label = 'Amplitude & Means fit')
    ax1.set_ylim(bottom = 0)
    ax1.set_title('Amplitude fit')
    ax2.set_title('Amplitude & Means fit')
    pyplot.xlabel('Resolution')
    fig.text(0.04,0.5,'Dataset score',va='center', rotation='vertical')
    pyplot.savefig(os.path.join(out_dir,'Resolution_dataset_score.png'))
    pyplot.close()

def plot_RMSD_vs_dataset_score(out_dir,all_RMSD,dataset_score,bound_ligands):
    
    fig,ax1 = myfig()
    ax1.scatter(all_RMSD.sum(axis = 1), dataset_score,
                label = 'No ligand bound datasets',marker = '+', 
                color  = 'b')
    ax1.scatter(all_RMSD.sum(axis = 1).loc[bound_ligands.columns.values],
                dataset_score.loc[bound_ligands.columns.values], color = 'k',
                label = 'Ligand bound dataset')
    ax1.set_title('Amplitude Fit')
    ax1.set_ylim(bottom =0)
    ax1.set_xlabel('Sum of RMSD across residues for each dataset')
    ax1.set_ylabel('Dataset Score')
    ax1.legend(loc = 'best')
    pyplot.savefig(os.path.join(out_dir,'RMSD_vs_Dataset_score.png'), dpi = 300)
    pyplot.close()

def plot_RMSD_vs_resolution(out_dir,all_RMSD,dataset_resolution):
    
    pyplot.scatter(dataset_resolution,all_RMSD.sum(axis = 1))
    pyplot.xlabel('Resolution')
    pyplot.ylabel('Sum of RMSD over all residues for each dataset')
    pyplot.savefig(os.path.join(out_dir,'Resolution_RMSD.png'))
    pyplot.close()

def plot_RMSD_vs_residue_score(out_dir,all_RMSD,residue_score,residue_means_score):
    fig, (ax1,ax2) = pyplot.subplots(2, sharex=True,sharey=True,
                                         figsize=(10,20))
    ax1.scatter(all_RMSD.sum(axis = 0),residue_score,
                   label = 'Amplitude Fit')
    ax1.set_title('Amplitude Fit')
    ax1.set_ylim(bottom =0)
    ax1.set_xlim(left =0)
    ax2.scatter(all_RMSD.sum(axis = 0),residue_means_score,
                   label = 'Amplitude & Means Fit')
    ax2.set_title('Amplitude & Means fit')
    ax2.set_ylim(bottom = 0)
    ax2.set_xlim(left = 0)
    pyplot.xlabel('Sum of RMSD over datasets for each residue')
    ax1.set_ylabel('Residue score')
    ax2.set_ylabel('Residue score')
    pyplot.savefig(os.path.join(out_dir,'RMSD_vs_residue_score.png'))
    pyplot.close()

def RMSD_histogram(all_RMSD,out_dir,RMSD_filename,fit_type):
    
    RMSD_hist_filename = "RMSD_hist_{}.png".format(fit_type)
    # If all_RMSD is null, i.e. the earlier part of function has not run this 
    # instance, due to file already existing
    if pandas.isnull(all_RMSD).any().any():
        all_RMSD = pandas.read_csv(os.path.join(out_dir,RMSD_filename),
                                   header =0, index_col =0)
    if not os.path.exists(os.path.join(out_dir,RMSD_hist_filename)):

        fig, ax =  myfig()

        bins = numpy.arange(0,1.05,0.05)
        counts, bins, patches = pyplot.hist(numpy.clip(numpy.concatenate(
                                            all_RMSD.values),bins[0],bins[-1]),
                                            bins=bins)
        # set xtick labels, to align with values clipped at 1.0
        xticks =[str(b) for b in bins[1:]]
        xticks[-1] = '>1.0'

        num_ticks=len(xticks)
        pyplot.xlim([0,1.05])
        pyplot.xticks(0.05*numpy.arange(num_ticks)+0.025)
        ax.set_xticklabels(xticks)

        pyplot.xlabel('RMSD {}  vs ringer data'.format(fit_type))
        pyplot.ylabel('Frequency')
        pyplot.savefig(os.path.join(out_dir,RMSD_hist_filename))
        pyplot.close()
    else:
        logger.info('RMSD already generated for fit type: {}'.format(fit_type))

def peak_angle_histogram(max_peak_angle,out_dir):
    
        fig, ax = myfig()
        counts, bins, patches = pyplot.hist(numpy.concatenate(
                                            max_peak_angle,axis=0),
                                            bins=60, normed=True)
        ax.set_xlim([0,360])
        ax.set_label('Angle')
        ax.set_ylabel('Relative Frequency')
        pyplot.title('Angles with maximal map value')
        pyplot.savefig(os.path.join(out_dir,"Modal_peak_location"))
        pyplot.close(fig)

def number_clusters_histogram(num_cluster_all,cluster_number_hist,out_dir):
    
    fig,ax=pyplot.subplots()

    n, bins, patches= pyplot.hist(num_cluster_all)

    pyplot.xlabel('Number of clusters')
    pyplot.ylabel('Number of residues')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.grid(False)
    
    pyplot.savefig(os.path.join(out_dir,cluster_number_hist))
    pyplot.close()

def plot_fit(interpolated_angles,interpolated_row,y_fit,angle_type,map_type,
             residue,out_dir,dataset,fit_type):

    fig, ax = myfig()
        
    ax.scatter(interpolated_angles,interpolated_row)
    ax.plot(interpolated_angles, y_fit, label= '3 gaussian fit')
    ax.set_xlim([0,360])
    pyplot.xlabel('Angle:{}'.format(angle_type))
    pyplot.ylabel('Map Value:{}'.format(map_type))
    pyplot.title('{} from {}'.format(residue,dataset))
    ax.legend(loc = 'upper right')

    #output image filename        
    filename='{}_{}_multi_gaussian_{}.png'.format(residue,dataset,
                                                          fit_type)
    # output directory
    fig.savefig(os.path.join(out_dir,residue,filename))
    pyplot.close()

def pairwise_heatmap(pairwise_csv, residue, out_dir,pairwise_type,
                     params, plot_filename,fit_type='', subset='',
                     min_scale = -1, max_scale = 1):
    """Plot heatmap from csv file with pairwise comparisons between datasets"""
 
    # Load csv into pandas DataFrame
    data = pandas.read_csv(pairwise_csv, index_col=0)
    assert (len(data) == len(params.input.dir)),(
           'Input CSV data is length {} for {} datasets.'
           'Lengths should match'.format(len(data)+1,len(params.input.dir)))

    dataset_labels = data.columns.values

    # Make figure & plot: scale minimum & maximum to minimum across all datasets
    
    fig, ax = pyplot.subplots()
    heatmap= ax.pcolor(data.values,cmap='RdBu',vmin = min_scale,vmax = max_scale)

    # Set title font size
    title_font_size = 0.1*int(len(params.input.dir))
    if title_font_size < 15:
        title_font_size =15

    # Scale title fontsize to numbe of datasets
    pyplot.title('{}'.format(residue),fontsize=title_font_size)

    # Format
    fig = pyplot.gcf()
    # Scale image size to dataset (with a minimum size)
    image_size = 0.1*int(len(params.input.dir))
    if image_size < 12:
        image_size = 12

    fig.set_size_inches(image_size,image_size)

    #Turn off the frame
    ax.set_frame_on(False)

    # Set tick positions
    ax.set(xticks=numpy.arange(len(dataset_labels))+0.5,xticklabels=dataset_labels,yticks = numpy.arange(len(dataset_labels))+0.5,yticklabels=dataset_labels)      #ax.set_xticks()

    #Set a table like display
    ax.invert_yaxis()

    # Set labels
    column_labels=dataset_labels
    row_labels=dataset_labels

    # Rotate x ticks
    pyplot.xticks(rotation=90)

    #set axis font size
    axis_font_size = 0.1*len(params.input.dir)
    if axis_font_size < 15:
        axis_font_size = 15

    # Colourbar
    cbar = fig.colorbar(heatmap,shrink=0.5)
    cbar.ax.tick_params(labelsize= axis_font_size)
    cbar.set_label('Pairwise_{}'.format(pairwise_type),rotation=270,fontsize = axis_font_size, labelpad = 40)
    #Output file
    filename = "{}_{}_{}_{}_heatmap.png".format(residue,pairwise_type,fit_type,subset)
    pyplot.savefig(os.path.join(out_dir,residue,filename))
    pyplot.close(fig)


def cluster_heatmap(input_csv, out_dir, pairwise_type,fit_type, subset):
    """Plot heatmap from csv file with pairwise comparisons between datasets"""

    # Load csv into pandas DataFrame
    data = pandas.read_csv(input_csv, index_col=0)

    dataset_labels = data.columns.values
    residue_labels = data.index.values

    # Make figure & plot: scale minimum & maximum to minimum across all datasets
    fig, ax = pyplot.subplots()
    heatmap= ax.pcolormesh(data.values,cmap='binary',vmin = 0,vmax = 1)

    # Set title font size
    title_font_size = 0.1*int(len(data.columns.values))
    if title_font_size < 15:
        title_font_size =15

    # Scale title fontsize to numbe of datasets
    pyplot.title('Cluster weights',fontsize=title_font_size)
    pyplot.xlabel('Datasets',fontsize=title_font_size)
    pyplot.ylabel('Residues', fontsize=title_font_size)

    # Format
    fig = pyplot.gcf()
    # Scale image width to dataset (with a minimum size)
    image_width = 0.1*int(len(data.columns.values))
    if image_width < 12:
        image_width = 12

    # Scale image height to dataset (with a minimum size)
    image_height = 0.1*int(len(data.index.values))
    if image_height < 12:
        image_height = 12
    fig.set_size_inches(image_width,image_height)

    #Turn off the frame
    ax.set_frame_on(False)

    # Set tick positions
    ax.set(xticks=numpy.arange(len(dataset_labels))+0.5,xticklabels=dataset_labels,yticks = numpy.arange(len(residue_labels))+0.5,yticklabels=residue_labels)      #ax.set_xticks()

    #Set a table like display
    ax.invert_yaxis()

    # Set labels
    column_labels=dataset_labels
    row_labels=residue_labels

    # Rotate x ticks
    pyplot.xticks(rotation=90)

    #set axis font size
    axis_font_size = 0.1*len(data.columns.values)
    if axis_font_size < 15:
        axis_font_size = 15

    # Colourbar
    cbar = fig.colorbar(heatmap,shrink=0.5)
    cbar.ax.tick_params(labelsize= axis_font_size)
    cbar.set_label('Cluster weight',rotation=270,fontsize = axis_font_size)
    # Output file
    filename = "Adj_cluster_weight_heatmap_{}_{}_{}.png".format(pairwise_type,fit_type,subset)

    pyplot.savefig(os.path.join(out_dir,filename))
    pyplot.close(fig)
