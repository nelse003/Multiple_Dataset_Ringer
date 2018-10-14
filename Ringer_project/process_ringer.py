
##########################################################
#Packages

import os
import sys
import copy
import glob

# For command manager         
from bamboo.common.command import CommandManager

# Package for getting summary from crystal (mtz): resolution
from iotbx.reflection_file_reader import any_reflection_file
from cctbx import miller
from giant.jiffies.quick_insert_f000 import run as insert_f000
from giant.jiffies.quick_insert_f000 import master_phil as f000_phil


def process_with_ringer(pdb, mtz, angle_sampling, output_dir=None,
                        output_base=None, column_labels="FWT,PHWT"):

    """Analyse a pdb-mtz pair with mmtbx.ringer

    Parameters
    -----------


    Returns
    -----------


    """

    assert os.path.exists(pdb), 'PDB File does not exist'
    assert os.path.exists(mtz), 'MTZ File does not exist'

    if not output_dir:  output_dir = os.path.dirname(pdb)
    if not output_base: output_base = os.path.splitext(os.path.basename(pdb))[0]

    # Check/create output directory
    if not os.path.exists(output_dir): os.mkdir(output_dir)

    output_csv =os.path.join(output_dir,output_base + '.csv')

    ###################################################
    #Run Ringer with F000 adjusted ,map
    ###################################################

    # When providing absolute electron density values
    abs_mtz = mtz.replace('.mtz', '.with-F000.mtz')

    # Only run if results don't already exist
    if not os.path.exists(os.path.join(output_csv)):
        if not os.path.exists(abs_mtz):

            f000_params = f000_phil.extract()
            f000_params.input.pdb = pdb
            f000_params.input.mtz = mtz
            f000_params.options.column.label = column_labels

            print(f000_params)

            insert_f000(f000_params)


        # Initialise and populate command object
        ringer = CommandManager(program='mmtbx.ringer')
        ringer.add_command_line_arguments(pdb, abs_mtz)
        ringer.add_command_line_arguments("scaling=volume")
        ringer.add_command_line_arguments('sampling_angle={}'.format(angle_sampling))
        ringer.add_command_line_arguments('output_base={}'.format(os.path.join(output_dir, output_base)))
        ringer.add_command_line_arguments("map_label=2FOFCWT-ABS,PHI2FOFCWT-ABS")
        ringer.add_command_line_arguments('difference_map_label=None')

        # Print and run
        ringer.print_settings()
        ringer.run()
        # Write log
        ringer.write_output(os.path.join(output_dir, output_base+'.log'))

    # Check the output csv file exists
    assert os.path.exists(output_csv), 'Ringer output CSV does not exist: {}'.format(output_csv)

    return output_csv

