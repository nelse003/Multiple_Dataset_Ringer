#!/usr/bin/env pandda.python

###############################################################################
# Packages
##############################################################################
import os, sys, copy, glob

import pandas, numpy
################################################################################
#                                   FUNCTIONS                                  #
################################################################################
def correlation_single_residue(input_csv,residue,output_dir,params,
                               out_filename):
    """Generate correlation matrix as csv"""
    # read in csv 
    data = pandas.read_csv(input_csv,index_col=0)
    assert (len(data) == len(params.input.dir)),(
           'Input CSV data is length {} for {} datasets.'
           'Lengths should match'.format(len(data)+1,len(params.input.dir)))

    dataset_labels=data.index.values
    all_map_values=data.values

    # Generate correlation coefficents
    correlation_matrix = numpy.corrcoef(all_map_values)
    # Store in labelled data frame
    correlation_data = pandas.DataFrame(correlation_matrix, 
                                        index = dataset_labels, 
                                        columns = dataset_labels)
    # Correlation data as CSV 
    correlation_data.to_csv(os.path.join(output_dir,out_filename))
    

