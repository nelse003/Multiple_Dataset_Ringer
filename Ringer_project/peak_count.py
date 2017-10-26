#!/usr/bin/env pandda.python

import os,sys,copy,glob
import pandas
import peakutils
import numpy as np

import sys
from numpy import NaN, Inf, arange, isscalar, asarray, array

from detect_peaks import detect_peaks
from scipy.signal import find_peaks_cwt
from peakdetect import peakdetect

###############################################################################
# Logging
###############################################################################
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

def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

    
###############################################################################
# Test detect peaks on residues
###############################################################################
def find_peaks(out_dir,ref_set,params,map_type,angle_type,peak_filename,
               threshold = 0.0):

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
                                                'Peak Locations','Peak_value'])
        for dir in params.input.dir:
            dataset_label = os.path.basename(dir)
            interpolated_map_value = interpolated_results.loc[dataset_label]
            peaks, dips = peakdetect(interpolated_map_value,lookahead=10, delta =0.1)
            peak_index = []
            peak_value = []
           
            for peak in peaks:
                if peak[1] < threshold:
                    peaks.remove(peak)
            for peak in peaks:
                peak_index.append(peak[0])
                peak_value.append(peak[1])
            #get angle instead of positional value
            peak_index=interpolated_angles[peak_index]

            num_peaks = len(peak_index)
            all_peaks.loc[dataset_label] = residue,num_peaks,peak_index,peak_value

        all_residues_peaks = all_residues_peaks.append(all_peaks)

    all_residues_peaks.to_csv(os.path.join(out_dir, peak_filename))
###############################################################################
# Use Peakutils python package to find peaks
############################################################################### 

def find_peakutils(out_dir, ref_set, params, map_type, angle_type, 
                   peak_filename, min_dist = 0.0, threshold = 0.2):
    
    all_residues = ref_set.index.values 
    all_residues_peaks = pandas.DataFrame() 

    exception_counter = 0
    try_counter = 0

    for residue, data in ref_set.iterrows():
        interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue,
                          len(params.input.dir),map_type,angle_type)
        interpolated_results = pandas.read_csv(os.path.join(out_dir,residue,
                                               interpolate_csv),header=0,
                                               index_col=0)
        interpolated_angles = interpolated_results.columns.values
        interpolated_angles = interpolated_angles.astype(np.float)
        all_peaks = pandas.DataFrame(index = params.input.dir, 
                                     columns = ['Residue','Number of peaks',
                                                'Peak Locations','indexes',
                                                'CWT'])
        for dir in params.input.dir:
            dataset_label = os.path.basename(dir)
            interpolated_map_value = interpolated_results.loc[dataset_label]
            peak_loc = find_peaks_cwt(interpolated_map_value,np.arange(5,60))
            cwt_peaks = interpolated_map_value[peak_loc]
            from IPython import embed; embed()
            indexes = peakutils.indexes(interpolated_map_value,
                                        thres = threshold,
                                        min_dist = min_dist)
            try:
                peaks = peakutils.interpolate(interpolated_angles,
                                              interpolated_map_value, 
                                              ind = indexes)
            except:
                try:
                    expected_peaks =[60,180,300]
                    peaks = peakutils.interpolate(interpolated_angles,
                                                  interpolated_map_value,
                                                  ind =expected_peaks)
                except:
                    peaks = indexes
                    exception_counter+=1
             
            num_peaks = len(peaks)
            all_peaks.loc[dataset_label] = residue, num_peaks, np.array(map(str,peaks),dtype =object), indexes, cwt_peaks

        all_residues_peaks = all_residues_peaks.append(all_peaks)

    all_residues_peaks.to_csv(os.path.join(out_dir, peak_filename))        
    #fail_percentage = 100*(exception_counter/try_counter)
    #logger.info("Peak fitting with interpolation failed in {}%".format(fail_percentage))

