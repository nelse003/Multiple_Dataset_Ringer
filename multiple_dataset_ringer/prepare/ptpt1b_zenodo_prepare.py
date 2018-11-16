import os
import glob
import time
from itertools import izip

def run_all_phenix_map_to_structure_factors(input_dir,
                                            pandda_processed_dir,
                                            resolution_mtz_stlye="_mrflagsref_idxs.mtz",
                                            ccp4_map_style="-aligned-map.ccp4",
                                            output_mtz_style="_aligned_from_panddas_ccp4.mtz"):

    for dataset in os.listdir(input_dir):

        resolution_mtz = os.path.join(input_dir,
                                      dataset,
                                      dataset + resolution_mtz_stlye)
        ccp4_file = os.path.join(pandda_processed_dir,
                                 dataset,
                                 dataset + ccp4_map_style)
        output_mtz = os.path.join(input_dir,
                                  dataset,
                                  dataset + output_mtz_style)

        run_phenix_map_to_structure_factors(resolution_mtz=resolution_mtz,
                                            ccp4_file=ccp4_file,
                                            output_mtz=output_mtz)

def run_phenix_map_to_structure_factors(resolution_mtz,
                                        ccp4_file,
                                        output_mtz):

    if os.path.exists(output_mtz):
        return output_mtz

    d_min = get_high_resolution_limit(resolution_mtz)
    os.system("phenix.map_to_structure_factors {} "
              "d_min={} output_file_name={}".format(ccp4_file,
                                                    d_min,
                                                    output_mtz))
    if os.path.exists(output_mtz):
        return output_mtz
    else:
        raise ValueError("phenix map to structure "
                         "factors failed for {}".format(ccp4_file))

def prepare_all_missing_reflections(input_dir,
                                    mtz_style="*.mtz",
                                    pdb_style="*.pdb"):

    """Prepare all datasets with missing reflections
    
    Notes
    -----------
    
    This is trivially parallelisable and needs to be parallelised 
    
    """

    for dataset in os.listdir(input_dir):
        os.chdir(os.path.join(input_dir,dataset))

        mtz = os.path.join(input_dir,
                           dataset,
                           glob.glob(mtz_style)[0])

        pdb = os.path.join(input_dir,
                           dataset,
                           glob.glob(pdb_style)[0])

        prepare_missing_reflections(mtz, pdb)



def prepare_missing_reflections(mtz, pdb):
    """ Prepare mtz file with missing reflections up to 999A
    
    Parameters
    -------------
    mtz: str
        path to mtz file
    pdb: str
        path to pdb file
        
    Returns
    -------------
    output_mtz: str
        path to output mtz file
    """

    run_cad(mtz)
    output_mtz = run_phenix_maps(mtz, pdb)

    return output_mtz

def run_phenix_maps(mtz, pdb):
    """ Run phenix.maps.
    
    Recalculate the 2Fo-DFc coefficients with phenix.maps
    
    which will create a file called OUTPUT_map_coeffs.mtz 
    containing columns 2FOFCWT_fill and PH2FOFCWT_fill, 
    which are the 2FOFC amplitudes and phases with 
    missing values filled with their estimates
    
    Parameters
    -------------
    mtz: str
        path to mtz file
    pdb: str
        path to pdb file

    Returns
    ------------
    output_mtz: str
        path to output mtz file
    
    """
    output_mtz = pdb.rstrip(".pdb") + "_map_coeffs.mtz"

    if os.path.exists(output_mtz):
        return output_mtz

    os.system("phenix.maps {} {}".format(mtz,pdb))

    if os.path.exists(output_mtz):
        return output_mtz
    else:
        raise ValueError("phenix.maps failed to generate output mtz"
                         " {} from {} and {}".format(mtz, pdb, output_mtz))

def run_cad(mtz):
    """ Run CAD to add reflections to an mtz file
    
    Parameters
    -------------
    mtz: str
        path to mtz file
    
    Returns
    ------------
    mtz_with_reflections: str
        path to mtz file with added reflections
    
    """
    mtz_with_reflections = mtz.rstrip(".mtz") + "with_missing_reflections.mtz"

    if os.path.exists(mtz_with_reflections):
        return mtz_with_reflections

    high_res = get_high_resolution_limit(mtz)
    cad_script = "cad hklin1 {} hklout {} <<eof\n"\
                 " monitor BRIEF\n" \
                 " labin file 1 - \n" \
                 "  ALL\n"\
                 " resolution file 1 999.0 {}\n"\
                 "eof".format(mtz, mtz_with_reflections, high_res)
    os.system(cad_script)

    if os.path.exists(mtz_with_reflections):
        return mtz_with_reflections
    else:
        raise ValueError("cad script failed for {}".format(mtz))

def get_high_resolution_limit(mtz):
    """ Get high resolution limit for a mtz file
    
    Uses mtzdmp to parse mtz file for high resolution limit. 
    Requires CCP4 to be sourced.
    
    Parameters
    -----------
    mtz: str
        path to mtz file to get high resolution limit from

    Returns
    -----------
    high_res: float
        high resolution limit from mtz file
    
    """
    cwd = os.getcwd()
    os.chdir(os.path.dirname(mtz))

    mtzdmp_txt_file = "{}_dmp.txt".format(mtz.rstrip(".*"))
    os.system("mtzdmp {} > {}".format(mtz, mtzdmp_txt_file))
    os.chdir(cwd)

    with open(mtzdmp_txt_file) as mtzdmp_file:
        contents = mtzdmp_file.readlines()
        res_range_line_number = 0
        res_line = None
        for i, line in enumerate(contents):
            if " *  Resolution Range :" in line:
                res_range_line_number = i
                break
        res_line = contents[i+2]

    if res_line.split()[3] > res_line.split()[5]:
        return res_line.split()[5]
    else:
        raise ValueError("High resolution limit {} "
                         "is not smaller then "
                         "low res limit {}".format(res_line.split()[5],
                                                   res_line.split()[3]))


def prepare_datasets_folder(input_models_dir,
                            mtz_dir,
                            mtz_style="_mrflagsref_idxs.mtz",
                            pdb_style="-pandda-input.pdb"):

    """ Generate pandda input like dataset forlders using symlinks
    
    From a folder containing all pdb files, and one contain all mtz
    files create a directory strucutre with one folder per dataset
    with an mtz and pdb file
    
    Parameters
    -----------
    input_models_dir: str
        path to the input model directory
    mtz_dir: str
        path to mtz directory
    mtz_style: str
        filename common to all mtz
    pdb_style: str
        filename common to all pdbs

    Returns
    ----------
    None
    """

    pdb_datasets = set(d[0:11] for d in os.listdir(input_models_dir))
    mtz_datasets = set(d[0:11] for d in os.listdir(mtz_dir))

    mtz_style = mtz_style.rstrip("*").lstrip("*")
    pdb_style = pdb_style.rstrip("*").lstrip("*")

    for dataset in pdb_datasets.intersection(mtz_datasets):

        if not os.path.exists(os.path.join(input_dir,dataset)):
            os.mkdir(os.path.join(input_dir,dataset))

        os.symlink(os.path.join(mtz_dir, dataset + mtz_style),
                   os.path.join(input_dir, dataset,
                                dataset + mtz_style))

        os.symlink(os.path.join(input_models_dir, dataset + pdb_style),
                   os.path.join(input_dir, dataset, dataset + pdb_style))

def symlink_pdb_mtz_only(input_dir,
                         output_dir,
                         mtz_style="_mrflagsref_idxs.mtz",
                         pdb_style="-pandda-input.pdb",
                         link_pdb=True,
                         link_mtz=True):
    """Make symlinks to a new dataset folder for pdb and mtz files"""

    datasets = set(d[0:11] for d in os.listdir(input_dir))

    mtz_style = mtz_style.rstrip("*").lstrip("*")
    pdb_style = pdb_style.rstrip("*").lstrip("*")

    for dataset in datasets:

        if not os.path.exists(os.path.join(output_dir, dataset)):
            os.makedirs(os.path.join(output_dir, dataset))

        if link_mtz:
            os.symlink(os.path.join(input_dir, dataset, dataset + mtz_style),
                       os.path.join(output_dir, dataset,
                                    dataset + mtz_style))
        if link_pdb:
            os.symlink(os.path.join(input_dir, dataset, dataset + pdb_style),
                       os.path.join(output_dir, dataset, dataset + pdb_style))



def merge_ptpt1b(ground_state_model_dir, bound_state_model_dir, merged_dir):
    """ Run giant.merge_datasets on ptp1b models"""

    os.chdir(merged_dir)

    for ground_state_pdb, bound_state_pdb \
        in izip(os.listdir(ground_state_model_dir),
                os.listdir(bound_state_model_dir)):

        if ground_state_pdb.split('_')[0] == bound_state_pdb.split('_')[0]:

            ground_state_pdb_path = os.path.join(ground_state_model_dir,
                                                 ground_state_pdb)
            bound_state_pdb_path = os.path.join(bound_state_model_dir,
                                                bound_state_pdb)

            log = os.path.join(merged_dir,
                               "{}-merge-confs.log".format(
                                   ground_state_pdb.split('_')[0]))

            output_pdb = os.path.join(merged_dir,
                                      "{}_merged.pdb".format(
                                          ground_state_pdb.split('_')[0]))

            os.system('giant.merge_conformations input.major={} input.minor={} '
                      'output.log={} '
                      'output.pdb={}'.format(ground_state_pdb_path,
                                             bound_state_pdb_path,
                                             log,
                                             output_pdb))



# ground_state_model_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/ground_state_models"
# bound_state_model_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/bound_state_models"
# merged_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/merged_models"

#input_dir ="/media/nelse003/Data/ringer_test_set/PTP1B/datasets"
input_dir = "/hdlocal/home/enelson/PTP1B/datasets_aligned"

if not os.path.exists(input_dir):
    os.mkdir(input_dir)

# if not os.path.exists(merged_dir):
#     os.mkdir(merged_dir)

# input_models_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/pandda_input_models"
# mtz_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/mtzs"

input_models_dir = "/hdlocal/home/enelson/PTP1B/pandda_input_models"
mtz_dir = "/hdlocal/home/enelson/PTP1B/mtzs"

mtz = "/hdlocal/home/enelson/PTP1B/datasets/PTP1B-y0001/PTP1B-y0001_mrflagsref_idxs.mtz"
pdb = "/hdlocal/home/enelson/PTP1B/datasets/PTP1B-y0001/PTP1B-y0001-pandda-input.pdb"
ccp4_file = "/hdlocal/home/enelson/PTP1B/pandda_04_11_18_test2/processed_datasets/PTP1B-y0001/PTP1B-y0001-aligned-map.ccp4"
output_mtz = os.path.join(os.path.dirname(mtz),"PTP1B-y0001_aligned_from_panddas_ccp4.mtz")

# run_all_phenix_map_to_structure_factors(input_dir="/hdlocal/home/enelson/PTP1B/datasets",
#                                         pandda_processed_dir="/hdlocal/home/enelson/PTP1B/pandda_04_11_18_test2/processed_datasets",
#                                         resolution_mtz_stlye="_mrflagsref_idxs.mtz",
#                                         ccp4_map_style="-aligned-map.ccp4",
#                                         output_mtz_style="_aligned_from_panddas_ccp4.mtz")

# run_phenix_map_to_structure_factors(resolution_mtz=mtz,
#                                     ccp4_file=ccp4_file,
#                                     output_mtz=output_mtz)

symlink_pdb_mtz_only(input_dir="/hdlocal/home/enelson/PTP1B/datasets",
                     output_dir="/hdlocal/home/enelson/PTP1B/datasets_aligned",
                     mtz_style="-pandda-input_map_coeffs.mtz",
                     pdb_style="-pandda-input.pdb",
                     link_pdb=False)

# prepare_all_missing_reflections(input_dir,
#                                 mtz_style="*_mrflagsref_idxs.mtz",
#                                 pdb_style="*-pandda-input.pdb")







