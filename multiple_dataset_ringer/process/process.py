import os
import sys
import copy
import glob
import pandas as pd
import logging
from itertools import izip

# Package for getting summary from crystal (mtz): resolution
from iotbx.reflection_file_reader import any_reflection_file
from iotbx.pdb import hierarchy
from cctbx import miller

# For command manager
from bamboo.common.command import CommandManager

from giant.jiffies.quick_insert_f000 import run as insert_f000
from giant.jiffies.quick_insert_f000 import master_phil as f000_phil
from giant.structure.select import protein

logger = logging.getLogger(__name__)


def rename_chain_if_differs(pdb, ref_pdb):
    """ Rename protein pdb chain to that of reference pdb if different.

    Parameters
    ----------
    pdb: str
        path of pdb that is being compared to reference pdb
    ref_pdb: str
        path of reference pdb

    Returns
    -------
    pdb: str
        path to pdb, edited name iof pdb had the same chain

    Notes
    --------
    Runs phenix.pdbtools because initial attempt
    to use hierarchy editing failed.

    Example command is:

    phenix.pdbtools dimple.pdb rename_chain_id.old_id=D \
    rename_chain_id.new_id=A output.file_name="dimple_edited.pdb"

    https://www.phenix-online.org/documentation/reference/pdbtools.html#manipulations-on-a-model-in-a-pdb-file-including

    """
    orig_pdb = pdb
    edited_pdb = os.path.join(os.path.dirname(pdb), "dimple_edited.pdb")
    if os.path.exists(edited_pdb):
        pdb = edited_pdb

    pdb_in = hierarchy.input(file_name=pdb)
    ref_pdb_in = hierarchy.input(file_name=ref_pdb)

    for chain, ref_chain in izip(
        protein(pdb_in.hierarchy, copy=False).only_model().chains(),
        protein(ref_pdb_in.hierarchy, copy=False).only_model().chains(),
    ):

        if chain.id != ref_chain.id:

            os.system(
                "phenix.pdbtools {} rename_chain_id.old_id={}"
                " rename_chain_id.new_id={}"
                " output.file_name={}".format(
                    orig_pdb, chain.id, ref_chain.id, edited_pdb
                )
            )

            pdb = edited_pdb

    return pdb


def rename_residue_labels(ringer_results):
    """ Change order of residue labels in dataframe
    
    This is needed for the comparison:
    
    current_dataset_results = dataset_results.loc[(dataset_results.index == residue)
    
    Parameters
    ----------
    ringer_results: pandas.DataFrame
        DataFrame containing a single dataset ringer results

    Returns
    -------
    ringer_results: pandas.DataFrame
        DataFrame containing a single dataset ringer results

    """
    for i in range(0, len(ringer_results.index)):
        res_split = ringer_results.index.values[i].rsplit(" ")
        if len(res_split) > 2:
            ringer_results.index.values[i] = res_split[0] + " " + res_split[1]

    return ringer_results


def process_all_with_ringer(params):

    # Dictionary to store all of the
    # ringer results for each of the
    # datasets
    all_results = {}

    # Create an output directory if it doesn't already exist
    if not os.path.isdir(params.output.out_dir):
        os.makedirs(params.output.out_dir)

    # Generate ringer results & resolution information
    ref_pdb = None

    for qsub_number, dataset_dir in enumerate(params.input.dir):

        print(qsub_number,dataset_dir)

        print(os.path.realpath(dataset_dir))
        # TODO PDB is not a single file!!

        files = glob.glob(os.path.join(os.path.realpath(dataset_dir), params.input.pdb_style))

        print(type(files))
        print(files)

        files.extend(glob.glob(os.path.join(os.path.realpath(dataset_dir), params.input.mtz_style)))

        print(type(files))
        print(files)

        if not pdb:
            continue

        if not mtz:
            continue

        pdb = pdb[0]
        mtz = mtz[0]

        if ref_pdb is None:
            ref_pdb = pdb

        if params.settings.sample_only_ref_pdb is not None:
            pdb = ref_pdb

        pdb = rename_chain_if_differs(pdb, ref_pdb)

        if not os.path.exists(pdb):
            print(
                "Skipping dir: No PDB Files found in {} matching {}: {}".format(
                    dataset_dir, params.input.pdb_style, pdb
                )
            )
            continue

        if not os.path.exists(mtz):
            print(
                "Skipping dir: No MTZ Files found in {} matching {}: {}".format(
                    dataset_dir, params.input.mtz_style, mtz
                )
            )
            continue

        # Process dataset with ringer and convert results to DataFrame

        ringer_csv = process_with_ringer(
            pdb=pdb,
            mtz=mtz,
            angle_sampling=params.settings.angle_sampling,
            output_dir=dataset_dir,
            column_labels=params.input.column_labels,
            qsub=params.settings.qsub,
            qsub_number=qsub_number,
            tmp_dir=params.output.tmp_dir,
        )

    if params.settings.qsub:

        if not os.path.exists(params.output.tmp_dir):
            os.mkdir(params.output.tmp_dir)

        with open(
            os.path.join(params.output.tmp_dir, "ringer_master.sh"), "w"
        ) as ringer_master:
            ringer_master.write("# PBS -joe -N ringer_master")
            ringer_master.write("./ringer_$SGE_TASK_ID.sh")

        os.system(
            "qsub -t 1:{!s} -tc {!s} {}".format(
                str(len(params.input.dir) + 2),
                100,
                os.path.join(params.output.tmp_dir, "ringer_master.sh"),
            )
        )

    for dataset_dir in params.input.dir:

        # print(os.path.realpath(os.path.join(dataset_dir, params.input.pdb_style)))

        pdb = glob.glob(
            os.path.realpath(os.path.join(dataset_dir, params.input.pdb_style))
        )

        # TODO Allow to continue when pdb does not exist in folder

        pdb = pdb[0]

        print(pdb)

        dataset_label = os.path.basename(dataset_dir.rstrip("/"))
        output_base = os.path.splitext(os.path.basename(pdb))[0]
        ringer_csv = os.path.join(dataset_dir, output_base + ".csv")
        ringer_results = pd.DataFrame.from_csv(ringer_csv, header=None)
        ringer_results = rename_residue_labels(ringer_results)
        all_results[dataset_label] = ringer_results

    return all_results


def process_with_ringer(
    pdb,
    mtz,
    angle_sampling,
    output_dir=None,
    output_base=None,
    column_labels="FWT,PHWT",
    qsub=False,
    qsub_number=None,
    tmp_dir=None,
):

    """Analyse a pdb-mtz pair with mmtbx.ringer

    Parameters
    -----------
    pdb: str
        path to pdb file
    mtz: str
        path to mtz file
    angle_sampling: int
        angle sampling rate used in mmtbx.ringer
    output_dir: str
        path to output directory
    output_base: str
        string to use for naming the output file
    column_labels: str
        column labels used in mmtbx.ringer, these columns will be 
        extracted from mtz file a supplied path
    qsub: bool
        flag to indicate use of qsub
    qsub_number: int
        number for the qsub run
    tmp_dir: str
        path of temporary directory
        

    Returns
    -----------
    output_csv: str
        path to the output csv generated by mmtbx.ringer

    """

    assert os.path.exists(pdb), "PDB File does not exist"
    assert os.path.exists(mtz), "MTZ File does not exist"

    if not output_dir:
        output_dir = os.path.dirname(pdb)
    if not output_base:
        output_base = os.path.splitext(os.path.basename(pdb))[0]

    # Check/create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_csv = os.path.join(output_dir, output_base + ".csv")

    ###################################################
    # Run Ringer with F000 adjusted ,map
    ###################################################

    # When providing absolute electron density values
    abs_mtz = mtz.replace(".mtz", ".with-F000.mtz")

    # Only run if results don't already exist
    if not os.path.exists(os.path.join(output_csv)):
        if not os.path.exists(abs_mtz):

            f000_params = f000_phil.extract()
            f000_params.input.pdb = pdb
            f000_params.input.mtz = mtz
            f000_params.options.column.label = column_labels
            insert_f000(f000_params)

        # Initialise and populate command object
        ringer = CommandManager(program="mmtbx.ringer")
        ringer.add_command_line_arguments(pdb, abs_mtz)
        ringer.add_command_line_arguments("scaling=volume")
        ringer.add_command_line_arguments("sampling_angle={}".format(angle_sampling))
        ringer.add_command_line_arguments(
            "output_base={}".format(os.path.join(output_dir, output_base))
        )
        ringer.add_command_line_arguments("map_label=2FOFCWT-ABS,PHI2FOFCWT-ABS")
        ringer.add_command_line_arguments("difference_map_label=None")

        # Print and run
        ringer.print_settings()
        if qsub:
            with open(
                os.path.join(tmp_dir, "ringer_{}.sh".format(qsub_number + 1)), "w"
            ) as qsub_file:
                qsub_file.write("#!/bin/bash\n")
                qsub_file.write("module load phenix\n")
                qsub_file.write(ringer.as_command() + "\n")

        else:
            print("AAAAA")
            ringer.run()
            ringer.write_output(os.path.join(output_dir, output_base + ".log"))
            assert os.path.exists(
                output_csv
            ), "Ringer output CSV does not exist: {}".format(output_csv)

    return output_csv
