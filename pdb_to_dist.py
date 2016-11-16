import pconpy
import pandas as pd
import pylab
import numpy as np
import os

import matplotlib.patches as mpatches
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['savefig.facecolor'] = 'white'
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['savefig.edgecolor']= 'none'

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

def return_dist_matrix(pdb_fn, measure = 'CA', dist_thres = 8.0, chain_ids='A'):
	
    # Returns a list of biopython residues
    residues = pconpy.get_residues(pdb_fn, chain_ids)
	
    # Calculates a distance matrix (as a masked array)
    # See pconpy for options 
    dist_mat = pconpy.calc_dist_matrix(residues, measure, dist_thres)
	
    # returns masked array 
    # Boolean True/False if threshold given
    # Distances if no threshold given
    return dist_mat
	
def get_residue_list(pdb_fn,chain_ids = 'A'):

    """ Get list of residues form pdb file in format """

    residues = pconpy.get_residues(pdb_fn, chain_ids ='A')

    residue_list = []

    for residue in residues:
        res_name = residue.get_resname()
        res_id = residue.get_id()
        res_id = res_id[1]
        
        residue_name = str(res_id) + ' ' + str(res_name) + ' ' + chain_ids
        residue_list.append(residue_name)

    return residue_list


def expand_dist_mat(pdb_fn, measure = 'CA', dist_thres = 8.0, chain_ids='A'):
    
    dist_mat = return_dist_matrix(pdb_fn, measure, dist_thres)

    residue_list = get_residue_list(pdb_fn, chain_ids)

    expanded_matrix = pd.DataFrame(dist_mat, index = residue_list, 
                                   columns = residue_list)

    return expanded_matrix

def import_conformer_matrix(csv_file_path):
    
    conformer_matrix = pd.read_csv(csv_file_path, header = 0, index_col = 0)
    
    return conformer_matrix
    
def spatial_conformer_matrix(csv_file_path,pdb_fn, measure = 'CA',
                             dist_thres = 8.0, chain_ids='A'):

    residue_list = get_residue_list(pdb_fn,chain_ids)
    contact_matrix = expand_dist_mat(pdb_fn, measure = 'CA', dist_thres = 8.0, chain_ids='A')
    # Import conformer matrix, and fill uncalculated residues with NaN
    conformer_matrix = import_conformer_matrix(csv_file_path)
    conformer_matrix = conformer_matrix.reindex(residue_list)
    conformer_matrix = conformer_matrix.reindex(columns=residue_list)

    conformer_spatial_mat =  contact_matrix.multiply(conformer_matrix, fill_value = 0.0)

    return (conformer_spatial_mat, contact_matrix)

def print_conformer_matrix(out_dir, img_filename, conformer_spatial_mat = None,
                           contact_map= None,title=None):
    ax,fig = pylab.gca(), pylab.gcf()
    
    # Colour maps
    cmap1 = mpl.cm.Greys
    cmap2 = mpl.cm.bwr
    cmap1._init()
    cmap2._init()
    alphas = np.linspace(0,1, cmap2.N+3)
    cmap1._lut[:,-1] = alphas
    cmap2._lut[:,-1] = alphas

    
    # Axis labels
    ax.set_xlabel('Residues')
    ax.set_ylabel('Residues')

    if contact_map is None and conformer_spatial_map is None:
        logger.error('No matrix provided')
        return

    if contact_map is not None:
        map_obj1 = pylab.pcolormesh(contact_map.values,
                                   shading = 'Flat', edgecolors = "None", 
                                   cmap = cmap1)
        black_patch = mpatches.Patch(color = 'Black', label = 
                                'Contacting Residue')
        lgd = pylab.legend(handles = [black_patch],
                       bbox_to_anchor=(0.5,-0.2),loc = 'lower center')
        # Axis Limits
        ax.set_xlim(0,len(contact_map))
        ax.set_ylim(0,len(contact_map))
    if conformer_spatial_mat is not None:
        map_obj2 = pylab.pcolormesh(conformer_spatial_mat.values, 
                               shading= "Flat", edgecolors ="None", 
                               cmap = cmap2)
        
        red_patch = mpatches.Patch(color='Red', 
                               label='>1 alternate conformer & contacting')
        lgd = pylab.legend(handles = [red_patch],
                       bbox_to_anchor=(0.5,-0.2),loc = 'lower center')
          
        # Axis Limits
        ax.set_xlim(0,len(conformer_spatial_mat))
        ax.set_ylim(0,len(conformer_spatial_mat))
    
    if title is not None:
        ax.set_title(title, fontweight = "bold")

    if contact_map is not None and conformer_spatial_mat is not None:
        lgd = pylab.legend(handles=[red_patch,black_patch],
                       bbox_to_anchor=(0.5,-0.3),loc = 'lower center')
    
    art =[]
    art.append(lgd)
    pylab.savefig(os.path.join(out_dir,img_filename),dpi=600, 
                  additional_artists=art, bbox_inches = 'tight')

    logger.info("Generated image:" + img_filename)

    pylab.close()

def print_all_alt_conformer(csv_file_path, pdb_fn, out_dir, img_filename, 
                           title = None, measure = 'CA', dist_thres = 8.0, 
                           chain_ids='A'):

    residue_list = get_residue_list(pdb_fn,chain_ids)
    conformer_matrix = import_conformer_matrix(csv_file_path)
    conformer_matrix = conformer_matrix.reindex(residue_list)
    conformer_matrix = conformer_matrix.reindex(columns=residue_list)
   
    conformer_matrix = conformer_matrix.fillna(value = 0.0)
     
    ax,fig = pylab.gca(), pylab.gcf()
    
    # Colour maps
    cmap1 = mpl.cm.Greys
    cmap1._init()
    alphas = np.linspace(0,1, cmap1.N+3)
    cmap1._lut[:,-1] = alphas

    # Axis Limits
    ax.set_xlim(0,len(conformer_matrix))
    ax.set_ylim(0,len(conformer_matrix))
    
    map_obj1 = pylab.pcolormesh(conformer_matrix.values,shading = 'Flat', 
                                edgecolors = "None",cmap = cmap1)

    black_patch = mpatches.Patch(color = 'Black', label = 
                                '<1 alternate conformer')
    art =[]
    lgd = pylab.legend(handles =[black_patch],
                       bbox_to_anchor=(0.5,-0.3),loc = 'lower center')
    art.append(lgd)
    pylab.savefig(os.path.join(out_dir,img_filename),dpi=300, 
                  additional_artists=art, bbox_inches = 'tight')

    logger.info("Generated image:" + img_filename)
    pylab.close()
