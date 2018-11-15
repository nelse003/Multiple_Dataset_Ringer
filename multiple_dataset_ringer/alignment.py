import os
import copy
from giant.structure.select import protein
from giant.structure.align import align_structures_flexible
from giant.dataset import CrystallographicModel

def align_ref_to_every_model(input_folder,
                             output_folder,
                             ref_model_path=None):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Select a reference model if not provided
    if ref_model_path is None:
        ref_model_folder = os.path.join(input_folder,
                                        os.listdir(input_folder)[0])
        models = [f for f in os.listdir(ref_model_folder) if f.endswith('.pdb')]
        ref_model_path = os.path.join(ref_model_folder, models[0])

    # Read in reference model
    ref_model = CrystallographicModel.from_file(ref_model_path)

    dataset_folders = [os.path.join(input_folder, f)
                       for f in os.listdir(input_folder)
                       if os.path.isdir(os.path.join(input_folder,f))]
    failures=[]
    # Loop over all models
    for dataset_folder in dataset_folders:

        print(os.path.basename(dataset_folder),
              os.path.basename(os.path.dirname(ref_model_path)))

        output_dataset_folder = os.path.join(output_folder,
                                             os.path.basename(dataset_folder))

        if not os.path.exists(output_dataset_folder):
            os.makedirs(output_dataset_folder)

        output_file_name = os.path.join(output_dataset_folder,
                                        "aligned_ref_{}_to_{}.pdb".format(
                                            os.path.basename(
                                                os.path.dirname(ref_model_path)),
                                            os.path.basename(dataset_folder)))

        if os.path.exists(output_file_name):
            continue

        # Read in comparision model
        models = [f for f in os.listdir(dataset_folder) if f.endswith('.pdb')]
        model_path = os.path.join(dataset_folder, models[0])
        model = CrystallographicModel.from_file(model_path)

        # Deepcopy reference model
        local_ref_model = copy.deepcopy(ref_model)

        # Align reference model to current model
        try:
            local_ref_model.align_to(model.hierarchy)
        except AssertionError:
            failures.append(os.path.basename(dataset_folder))
            continue
        # Write pdb file of reference data moved
        local_ref_model.hierarchy.write_pdb_file(
            file_name=output_file_name,
            crystal_symmetry=local_ref_model.crystal_symmetry)

    print(failures)

ref_model_path = "/hdlocal/home/enelson/PTP1B/datasets/PTP1B-y0001/PTP1B-y0001-pandda-input.pdb"

align_ref_to_every_model(input_folder="/hdlocal/home/enelson/PTP1B/datasets",
    output_folder="/hdlocal/home/enelson/PTP1B/datasets_aligned",
    ref_model_path=ref_model_path)

#failures=['PTP1B-y0112', 'PTP1B-y1904', 'PTP1B-y1608', 'PTP1B-y1646', 'PTP1B-y1320', 'PTP1B-y1781', 'PTP1B-y1548', 'PTP1B-y1486', 'PTP1B-y0426', 'PTP1B-y1043', 'PTP1B-y0572', 'PTP1B-y1594', 'PTP1B-y1787', 'PTP1B-y1159', 'PTP1B-y0884', 'PTP1B-y0654', 'PTP1B-y0117', 'PTP1B-y0211', 'PTP1B-y0725', 'PTP1B-y1525', 'PTP1B-y0965', 'PTP1B-y0678', 'PTP1B-y0772', 'PTP1B-y0072', 'PTP1B-y1589', 'PTP1B-y0205', 'PTP1B-y1819', 'PTP1B-y1288', 'PTP1B-y1271', 'PTP1B-y0911', 'PTP1B-y0363', 'PTP1B-y1922', 'PTP1B-y1373', 'PTP1B-y0845', 'PTP1B-y1938', 'PTP1B-y1465', 'PTP1B-y1312', 'PTP1B-y0847', 'PTP1B-y0829', 'PTP1B-y0085', 'PTP1B-y0111', 'PTP1B-y0718', 'PTP1B-y0522', 'PTP1B-y0105', 'PTP1B-y1759', 'PTP1B-y0676', 'PTP1B-y1554', 'PTP1B-y1885', 'PTP1B-y0656', 'PTP1B-y0538', 'PTP1B-y1711', 'PTP1B-y1438', 'PTP1B-y0076', 'PTP1B-y1417', 'PTP1B-y1153', 'PTP1B-y0972', 'PTP1B-y0175', 'PTP1B-y0070', 'PTP1B-y1710', 'PTP1B-y0652', 'PTP1B-y0517', 'PTP1B-y0296', 'PTP1B-y0193', 'PTP1B-y0571', 'PTP1B-y0116', 'PTP1B-y1264', 'PTP1B-y0660', 'PTP1B-y1009', 'PTP1B-y0528', 'PTP1B-y0766', 'PTP1B-y0846', 'PTP1B-y0579', 'PTP1B-y0446', 'PTP1B-y0963', 'PTP1B-y0077', 'PTP1B-y0049', 'PTP1B-y1294', 'PTP1B-y1703', 'PTP1B-y0877', 'PTP1B-y1512', 'PTP1B-y0591', 'PTP1B-y1618', 'PTP1B-y0812', 'PTP1B-y1103', 'PTP1B-y0118', 'PTP1B-y0056', 'PTP1B-y1842', 'PTP1B-y0249', 'PTP1B-y1304', 'PTP1B-y1298', 'PTP1B-y0573', 'PTP1B-y0816', 'PTP1B-y0822', 'PTP1B-y0648', 'PTP1B-y1314', 'PTP1B-y1652', 'PTP1B-y1125', 'PTP1B-y1136', 'PTP1B-y0180', 'PTP1B-y1492', 'PTP1B-y0580', 'PTP1B-y1629', 'PTP1B-y0518', 'PTP1B-y0650', 'PTP1B-y1335', 'PTP1B-y0559', 'PTP1B-y0914', 'PTP1B-y0071', 'PTP1B-y0835', 'PTP1B-y1865', 'PTP1B-y0375', 'PTP1B-y1257', 'PTP1B-y1943', 'PTP1B-y1011', 'PTP1B-y1957', 'PTP1B-y1612', 'PTP1B-y1763']



