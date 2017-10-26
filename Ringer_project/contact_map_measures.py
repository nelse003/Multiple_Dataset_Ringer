#/usr/bin/env_pandda.python

import pconpy
import os

from pdb_to_dist import (expand_dist_mat, print_conformer_matrix) 

def measure_test():
 
    datasets = ['270','283','463']

    measures ={'CA': 8.0, 
               'CB': 8.0,
               'cmass': 8.0,
               'sccmass': 8.0 , 
               'minvdw': 8.0 }

    out_dir = '/hdlocal/home/enelson/test-data/JMJD2DA-Zenobia-Screen/Alternate_conformer_maps/measure_test'

    for dataset in datasets:
        pdb_filepath = '/hdlocal/home/enelson/test-data/JMJD2DA-Zenobia-Screen/JMJD2DA-x' + dataset
        pdb_filename = 'JMJD2DA-x' + dataset + '.dimple.pdb'
        pdb_file = os.path.join(pdb_filepath,pdb_filename)
        img_filename_base = 'Contact_map_JMJD2DA-x' + dataset

        for measure in measures.keys():
            img_filename = measure + '_' + img_filename_base

            contact_matrix = expand_dist_mat(pdb_file, measure, 
                                             dist_thres = measures.get(measure))
            print_conformer_matrix(out_dir, img_filename, 
                                   contact_map = contact_matrix, 
                                   title = measure + ' Thres: ' + 
                                           str(measures.get(measure)) + 
                                           'Dataset: ' + dataset) 
