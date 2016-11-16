#!/usr/bin/env pandda.python

import os,sys,copy,glob
import pandas

from detect_peaks import detect_peaks


###############################################################################
# Test detect peaks on residues
###############################################################################
def find_peaks(out_dir,ref_set,params,map_type,angle_type,peak_filename,
               mpd = 1,threshold = 0.0, mph = None):

    all_residues = ref_set.index.values 

    all_residues_peaks = pandas.DataFrame() 

    for residue, data in ref_set.iterrows():
        interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue,
                          len(params.input.dir),map_type,angle_type)
        interpolated_results = pandas.read_csv(os.path.join(out_dir,residue,
                                               interpolate_csv),header=0,
                                               index_col=0)
        interpolated_angles = interpolated_results.columns.values
        all_peaks = pandas.DataFrame(index = params.input.dir, 
                                     columns = ['Residue','Number of peaks',
                                                'Peak Locations'])
        for dir in params.input.dir:
            dataset_label = os.path.basename(dir)
            interpolated_map_value = interpolated_results.loc[dataset_label]
            peaks = interpolated_angles[detect_peaks(
                    interpolated_map_value, mpd, mph, threshold)]
            num_peaks = len(peaks)
            all_peaks.loc[dataset_label] = residue,num_peaks,peaks

        all_residues_peaks = all_residues_peaks.append(all_peaks)

    all_residues_peaks.to_csv(os.path.join(out_dir, peak_filename))
 

