#!/usr/bin/env pandda.python

#######################################################################
# Packages
#######################################################################
import os
import sys
import copy
import glob
import pandas
import numpy
import time
import logging

from scipy.spatial.distance import euclidean
from scipy.optimize import curve_fit
from sklearn.metrics import mean_squared_error

# Plotting functions import
from plotting_ringer import RMSD_histogram
from plotting_ringer import plot_fit

#######################################################################
# Logging
#######################################################################

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

#######################################################################
# Fitting Functions
#######################################################################


def single_gaussian(x, a, x0, sigma):
    """Gaussain to fit data to (3 parameters)"""
    return a*numpy.exp(-(x-x0)**2/(2*sigma**2))


def three_gaussian_offset(x, a1, x01, sigma_1, a2, x02, sigma_2, a3, x03,
                          sigma_3, offset):

    """Three gaussian curve with offset (10 parameters)"""

    return (single_gaussian(x, a1, x01, sigma_1) +
            single_gaussian(x, a2, x02, sigma_2) +
            single_gaussian(x, a3, x03, sigma_3) +
            offset )


def three_gaussian_fix(x, a1, sigma_1, a2, sigma_2, a3, sigma_3, offset):

    """Three gaussian curve, fixed means [60,180,300] with offset. 
    
    Has 7 Parameters
    """
    return three_gausian_offset(x, a1, 60, sigma_1, a2, 180, sigma_2,
                                a3, 300, sigma_3, offset)


def three_normal_fix(x, sigma_1, sigma_2, sigma_3, offset):

    """Three normal curves,fixed means [60,180,300] with offset. 4 parameters"""

    amplitude = 1.0/math.sqrt((2*sigma_1**2)*math.pi)
    return three_gaussian_fix(x, amplitude, sigma_1, amplitude, sigma_2,
                              amplitude, sigma_3, offset)


def fit_all_datasets(out_dir, ref_set, map_type, angle_type, params,pairwise_type,
                    fit_type = 'three_gaussian_offset', subset= None,
                    mean_bound = None, datasets=None):

    """ Non Linear Least squares fitting routine, single & multiple gaussian
    
    Notes
    ------------
    
    Output: CSV file with fit parameters (one file per residue) 
    Images of the fitting (One per dataset per residue)  
    
    Returns
    ------------
    
    """

    for residue, data in ref_set.iterrows():

        parameters_csv_filename = '{}_from_{}_datasets_{}.csv'.format(residue,
                                len(datasets), fit_type)

        if not os.path.exists(os.path.join(out_dir,residue,parameters_csv_filename)):
            start =time.time()

            # Reading in interpolated angle data
            interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue,
                                                                       len(datasets),
                                                                       map_type,
                                                                       angle_type)
            interpolated_results=pandas.read_csv(os.path.join(out_dir,
                                                              residue,
                                                              interpolate_csv),
                                                 index_col=0)

            # select fit type
            if fit_type in ('three_gaussian_offset',
                            'positive_amplitude_three_gaussian_offset'):

                # Dataframe to store parameters for multiple datasets
                fit_parameters=pandas.DataFrame(columns=['a1','mean_1','sigma_1',
                            'a2','mean_2','sigma_2','a3','mean_3','sigma_3','offset'])

                initialise_three=[3, 60, 10, 3, 180, 10, 3, 300, 10, -1]

            elif fit_type == 'three_gaussian_fix':
                # Dataframe to store parameters for multiple datasets
                fit_parameters=pandas.DataFrame(columns=['a1','sigma_1','a2',
                               'sigma_2','a3','sigma_3','offset'])
                intialise_three=[3,10,3,10,10,-1]

            else:
                logger.info('fit_type is not a specified fit type')

            # Initial parameter to be used in first instance only      
            initialised = False

            # Bounds for means for three gaussian offset
            if mean_bound is not None:
                bounds = ([-numpy.inf,60-mean_bound,-numpy.inf,
                           -numpy.inf,180-mean_bound,-numpy.inf,
                           -numpy.inf,300-mean_bound,-numpy.inf,
                           -numpy.inf],
                          [numpy.inf,60+mean_bound,numpy.inf,
                           numpy.inf,180+mean_bound,numpy.inf,
                           numpy.inf,300+mean_bound,numpy.inf,
                           numpy.inf])
            else:
                bounds = (-numpy.inf,numpy.inf)

            for dataset in interpolated_results.index:

                interpolated_row = interpolated_results.loc[dataset].values
                interpolated_angles = interpolated_results.columns.values.astype(int)

                # Non linear least sqaures Curve fit on data:
                # popt returnd best fit values for parameters of the model
                # pcov contains covariance matrix for the fit
                if fit_type == 'three_gaussian_offset':

                    popt, pcov = curve_fit(f=three_gaussian_offset,
                                           xdata=interpolated_angles,
                                           ydata=interpolated_row,
                                           p0=initialise_three,
                                           ftol=1e-4,
                                           maxfev=20000)

                    # TODO work out development environment which allows
                    # TODO running scipy > 0.17

                    # Non compatible with scipy 0.16.1,
                    # due to keyword bounds. This is version distributed
                    # with ccp4-python and cannot easily be updated.
                    # Previously nick's version of panddas python
                    # had updated scipy

                    # popt, pcov = curve_fit(f=three_gaussian_offset,
                    #                        xdata=interpolated_angles,
                    #                        ydata=interpolated_row,
                    #                        p0=initialise_three,
                    #                        ftol=1e-4,
                    #                        bounds = bounds,
                    #                        method='trf')

                    y_fit = three_gaussian_offset(interpolated_angles, popt[0],
                                                  popt[1], popt[2], popt[3],
                                                  popt[4], popt[5], popt[6],
                                                  popt[7], popt[8], popt[9])

                if fit_type == 'three_gaussian_fix':
                    popt, pcov = curve_fit(three_gaussian_fix, interpolated_angles,
                                           interpolated_row, p0 = initialise_three,
                                           ftol=1e-4)
                    y_fit = three_gaussian_offset(interpolated_angles, popt[0],popt[1],
                                         popt[2], popt[3], popt[4], popt[5], popt[6])
                if fit_type == 'positive_amplitude_three_gaussian_offset':
                    # Set Amplitude lower bounds to zero
                    bounds[0][0]=0;bounds[0][3]=0;bounds[0][6]=0
                    # Run fit
                    popt, pcov = curve_fit(three_gaussian_offset, interpolated_angles,
                                           interpolated_row, p0 = initialise_three,
                                           ftol=1e-4, bounds = bounds)
                    # Fitted Data
                    y_fit = three_gaussian_offset(interpolated_angles, popt[0],popt[1],
                                         popt[2], popt[3], popt[4], popt[5],
                                         popt[6], popt[7], popt[8], popt[9])

                # Set Inital parameters to the first fit for all other fits    
                if not initialised:
                    intialise_three=popt
                    initialised = True

                # Plot_fit               
                plot_fit(interpolated_angles,interpolated_row,y_fit,angle_type,map_type,
                         residue,out_dir,dataset,fit_type)
                # CSV with parameters
                if fit_type in ('three_gaussian_offset',
                                'positive_amplitude_three_gaussian_offset'):
                    fit_df=pandas.DataFrame(popt.reshape(1,10),columns=['a1',
                            'mean_1','sigma_1','a2','mean_2','sigma_2','a3',
                            'mean_3','sigma_3','offset'], index =[dataset])

                if fit_type == 'three_gaussian_fix':
                    fit_df=pandas.DataFrame(popt.reshape(1,7),columns=['a1',
                            'sigma_1','a2','sigma_2','a3','sigma_3','offset']
                            ,index =[dataset])

                # Adding fit parameters to list for all datasets
                fit_parameters = fit_parameters.append(fit_df)
                fit_parameters.to_csv(os.path.join(out_dir, residue,
                                                   parameters_csv_filename))


            end=time.time()
            duration = end-start
            logger.info('{} fits  in {} seconds'.format(residue, duration))
        else:
            logger.info('{} Fitting already undertaken'.format(residue))


def generate_RMSD(out_dir, ref_set, map_type, angle_type, fit_type,
                  fit_base_filename, datasets):

    """ Generate RMSD between data and fits
     
    Generate RMSD between fit and data. 
    Plot histogram of all RMSD values.
    Store RMSD values in single CSV for the fit type.
     
     """

    # Dataframe to store all RMSD values
    all_RMSD = pandas.DataFrame(index=datasets, columns=ref_set.index.values)
    RMSD_filename = 'RMSD, with {}.csv'.format(fit_type)

    if not os.path.exists(os.path.join(out_dir, RMSD_filename)):
        for residue, data in ref_set.iterrows():

            fit_filename = residue + fit_base_filename

            # Pandas Dataframe for RMSD        
            RMSD_results = pandas.DataFrame(columns=['RMSD'])

            # Retrieve list of datasets
            interpolate_csv = '{}_{}_Datasets_{}_{}-ringer.csv'.format(residue,
                              len(datasets),map_type,angle_type)

            interpolated_results = pandas.read_csv(os.path.normpath(
                os.path.join(out_dir, residue, interpolate_csv)), index_col=0)

            interpolated_angles = interpolated_results.columns.values.astype(int)
            # Read in fit parameters
            fit_parameters = pandas.read_csv(os.path.normpath(os.path.join(out_dir,residue,
                                            fit_filename)),header = 0,index_col = 0)
            for dataset in interpolated_results.index:
                # fit parameters
                current_fit = fit_parameters.loc[dataset]

                # map value data that was fitted
                y = interpolated_results.loc[dataset].values

                # Turn fit parameters into fitted values for eacg fit type
                if fit_type in ('three_gaussian_offset',
                                'positive_amplitude_three_gaussian_offset'):
                    y_fit = three_gaussian_offset(interpolated_angles,current_fit[0],
                                                  current_fit[1], current_fit[2],
                                                  current_fit[3], current_fit[4],
                                                  current_fit[5], current_fit[6],
                                                  current_fit[7], current_fit[8],
                                                  current_fit[9])

                elif fit_type == "three_gaussian_fix":
                    y_fit = three_gaussian_fix(interpolated_angles,current_fit[0],
                                                  current_fit[1], current_fit[2],
                                                  current_fit[3], current_fit[4])
                else:
                    logger.warning("Fit type not recongised")

                # Generate RMSD between fit and data
                RMSD=mean_squared_error(y,y_fit)**0.5
                # RMSD check
                c=y_fit-y
                d=c.dot(c)/len(c)
                output=d**0.5
                numpy.testing.assert_almost_equal(output,RMSD,decimal=7,
                              err_msg="RMSD not calculated correctly")

                RMSD_result=pandas.DataFrame(data = RMSD, index = [dataset],
                                                columns =['RMSD'])
                RMSD_results = RMSD_results.append(RMSD_result)


            all_RMSD.loc[:,residue]= RMSD_results.values

        # Store RMSD for all residues
        all_RMSD.to_csv(os.path.join(out_dir, RMSD_filename))

    #RMSD Histogram
    RMSD_histogram(all_RMSD,out_dir, RMSD_filename, fit_type)


def calculate_euclidean_distance(out_dir, ref_set, params,
                                 fit_type, datasets, subset=None):

    '''Calculates and store pairwise euclidean distances'''

    euclidean_base_csv = 'euclidean_distances_from_{}_datasets_{}'.format(len(params.input.dir), fit_type)
    # Selecting a subset of parameters to perform clustering on
    if subset is not None:
        # set pairwise type so output filename changes 
        euclidean_base_csv = euclidean_base_csv.rsplit('.')[0] + subset + '.csv'

    for residue, data in ref_set.iterrows():
        euclidean_csv='{}'.format(residue) + euclidean_base_csv

        if not os.path.exists(os.path.join(out_dir,residue,euclidean_csv)):
            start=time.time()
            # Read data from csv
            parameters_csv_filename = '{}_from_{}_datasets_{}.csv'.format(
                                      residue,len(params.input.dir),fit_type)
            fit_parameters=pandas.read_csv(os.path.join(out_dir,residue,parameters_csv_filename),   
                                  index_col=0, header=0)

            assert (len(fit_parameters) == len(datasets)),(
                    'Input CSV data is length {} for {} datasets.'
                    'Lengths should match'.format(len(fit_parameters), len(params.input.dir)))

            euclidean_distance=pandas.DataFrame(index=fit_parameters.index, columns=fit_parameters.index)

            if subset == 'Amplitudes':
                fit_parameters = fit_parameters[['a1','a2','a3']]
            elif subset == 'Amplitudes_Means':
                fit_parameters = fit_parameters[['a1','a2','a3','mean_1','mean_2','mean_3']]
            else:
                logger.warning('Incorrect identifier for subset of parameters on which to perform clustering')

            #Calcualte euclidean distances 
            for dataset_col in fit_parameters.index:
                for dataset_row in fit_parameters.index:
                    euclidean_distance.loc[dataset_col][dataset_row] = euclidean(
                        fit_parameters.loc[dataset_col].values,fit_parameters.loc[dataset_row].values)

            # Output CSV
            euclidean_distance.to_csv(os.path.join(out_dir,residue,euclidean_csv))
            end=time.time()
            duration=end-start

            logger.info('{}: Euclidean distance '
                        'for {} datasets in '
                        '{} seconds'.format(residue,
                                            len(datasets),
                                            duration))
        else:
            logger.info('{}: Euclidean distance already calculated for these {} datasets'.format(residue,len(params.input.dir)))

    return euclidean_base_csv

