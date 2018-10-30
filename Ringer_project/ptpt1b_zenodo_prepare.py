import os
from itertools import izip

# ground_state_model_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/ground_state_models"
# bound_state_model_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/bound_state_models"
# merged_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/merged_models"

#input_dir ="/media/nelse003/Data/ringer_test_set/PTP1B/datasets"
input_dir = "/hdlocal/home/enelson/PTP1B/datasets_single_pdb"

if not os.path.exists(input_dir):
    os.mkdir(input_dir)

# if not os.path.exists(merged_dir):
#     os.mkdir(merged_dir)

# input_models_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/pandda_input_models"
# mtz_dir = "/media/nelse003/Data/ringer_test_set/PTP1B/mtzs"

input_models_dir = "/hdlocal/home/enelson/PTP1B/pandda_input_models"
mtz_dir = "/hdlocal/home/enelson/PTP1B/mtzs"

pdb_datasets = set(d[0:11] for d in os.listdir(input_models_dir))
mtz_datasets = set(d[0:11] for d in os.listdir(mtz_dir))

for dataset in pdb_datasets.intersection(mtz_datasets):

    if not os.path.exists(os.path.join(input_dir,dataset)):
        os.mkdir(os.path.join(input_dir,dataset))

    os.symlink(os.path.join(mtz_dir, dataset + "_mrflagsref_idxs.mtz"),
               os.path.join(input_dir, dataset,
                            dataset + "_mrflagsref_idxs.mtz"))

    os.symlink(os.path.join(input_models_dir, dataset + "-pandda-input.pdb"),
               os.path.join(input_dir, dataset, dataset + "-pandda-input.pdb"))



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







