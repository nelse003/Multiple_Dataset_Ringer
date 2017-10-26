
import os, sys
import pandas
import iotbx.pdb
from scitbx.array_family import flex

# Load input csv
input_data = pandas.DataFrame.from_csv('./score_matrix.csv')

# Load the input structure into the hierarchy object (google "cctbx hierarchy")
input_structure = iotbx.pdb.hierarchy.input('./input_structure.pdb').hierarchy
# Reset the B-factors in the structure to 0.0
input_structure.atoms().set_b(0.0*input_structure.atoms().extract_b())

# Create output folder
output_folder = './output_structures_score'
if not os.path.exists(output_folder): os.mkdir(output_folder)

# Iterate through the datasets
for dataset, residue_info in input_data.transpose().iterrows():
    print dataset

    # Get a new copy of the structure
    new_struc = input_structure.deep_copy()

    for rg in new_struc.residue_groups():
        # Create the right label for the residue
        label = '{} {} {}'.format(rg.resseq, rg.unique_resnames()[0], rg.parent().id).strip()
        if label not in residue_info.index:
            print 'No clustering to process for "{}"'.format(label)
            continue
        else:
            print 'Labelling: ', label, ' -', residue_info[label]
            rg.atoms().set_b(flex.double(rg.atoms_size(), residue_info[label]))

    new_struc.write_pdb_file(os.path.join(output_folder, '{}-labelled.pdb'.format(dataset)))


from IPython import embed; embed(); raise SystemExit()
