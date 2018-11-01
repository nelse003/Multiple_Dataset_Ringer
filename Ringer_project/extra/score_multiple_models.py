from giant.jiffies import score_model

import os, sys, copy, re, glob

#################################
import matplotlib
matplotlib.use('Agg')
matplotlib.interactive(0)
from matplotlib import pyplot
pyplot.style.use('ggplot')
#################################

import libtbx.phil

#################################
bar = '=======================++>'

PROGRAM = 'giant.score_multiple_model'
DESCRIPTION = """
    A tool to quickly score a ligand model against crystallographic electron density

    1) Simple usage (for ligand called LIG or UNL):
        > giant.score_model refined.pdb refined.mtz

    2) Score a residue in the file (replace XXX with ligand 3-letter ID)
        > giant.score_model ... res_names=XX1,XX2,XX3

    3) Define a "fitted" model to compare the refined model against
        > giant.score_model ... pdb1=fitted.pdb mtz1=fitted.mtz
"""
blank_arg_prepend = {'.pdb':'pdb1=', '.mtz':'mtz1='}

residue_plot_phil =  """
plot {
    remove_blank_entries = False
        .type = bool
    print_axis_values = True
        .type = bool
    parameters {
        rscc {
            title = 'Model\nQuality\n(RSCC)'
                .type = str
            axis_min = 0.60
                .type = float
            axis_max = 0.85
                .type = float
            axis_invert = True
                .type = bool
        }
        rszd {
            title = 'Model\nAccuracy\n(RSZD)'
                .type = str
            axis_min = 1.50
                .type = float
            axis_max = 4.00
               .type = float
            axis_invert = False
                .type = bool
        }
        rszo {
            title = 'Density\nPrecision\n(RSZO/OCC)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 2.00
                .type = float
            axis_invert = True
                .type = bool
        }
        b_factor_ratio {
            title = 'B-Factor\nStability\n(B-factor Ratio)'
                .type = str
            axis_min = 1.00
                .type = float
            axis_max = 3.00
                .type = float
            axis_invert = False
                .type = bool
        }
        rmsd {
            title = 'Coordinate\nStability\n(RMSD)'
                .type = str
            axis_min = 0.00
                .type = float
            axis_max = 1.50
                .type = float
            axis_invert = False
                .type = bool
        }
    }
}
"""

master_phil = libtbx.phil.parse("""
input {
    pdb1 = None
        .type = str
        .multiple = False
    mtz1 = None
        .type = str
        .multiple = False

    pdb2 = None
        .type = str
        .multiple = False
    mtz2 = None
        .type = str
        .multiple = False
    label = ''
        .type = str
        .multiple = False
}
selection {
    res_names = LIG,UNL,DRG
        .type = str
        .help = "Comma-separated list of residue names to score -- if None then scores all residues"
}
output {
    out_dir = ./
        .type = path
}
"""+residue_plot_phil, process_includes=True)

#######################################

def prepare_output_directory(params):
    if not os.path.exists(params.output.out_dir): os.mkdir(params.output.out_dir)
    images_dir = os.path.join(params.output.out_dir, 'residue_plots')
    if not os.path.exists(images_dir): os.mkdir(images_dir)
    return params.output.out_dir, images_dir

def run(params):

    assert params.input.pdb1, 'No pdb1 provided'
    assert params.input.mtz1, 'No mtz1 provided'
    # REMOVE THIS WHEN IMPLEMENTED
    assert params.selection.res_names
    # REMOVE THIS WHEN IMPLEMENTED
    if params.selection.res_names:      params.selection.__inject__("res_names_list", params.selection.res_names.split(','))
    else:                               params.selection.__inject__("res_names_list", None)

    output_dir, images_dir = prepare_output_directory(params)
    scores_file = os.path.join(output_dir, 'residue_scores.csv')

    print bar
    print 'Scoring model...'
    data_table = score_model(   params = params,
                                pdb1   = params.input.pdb1,
                                mtz1   = params.input.mtz1,
                                pdb2   = params.input.pdb2,
                                mtz2   = params.input.mtz2,
                                label_prefix = params.input.label
                            )
    print '...Done'
    print bar

    data_table.to_csv(scores_file)
    print 'Output written to {}'.format(scores_file)
    print bar

    ###################################################################
    # Image parameters
    ###################################################################
    columns = format_parameters_for_plot(params=params.plot.parameters)

    ###################################################################
    # Output Images
    ###################################################################
    all_images = []
    print 'Generating Output Images...'
    for label, row in data_table.iterrows():
        print 'Making: {}...'.format(label)
        image_path = os.path.join(images_dir,'{}.png'.format(label))
        make_residue_radar_plot(path = image_path,
                                data = row.to_frame().T,
                                columns = columns,
                                remove_blank_entries = params.plot.remove_blank_entries,
                                print_axis_values    = params.plot.print_axis_values  )
        all_images.append(image_path)

    from IPython import embed; embed()

    print '...Done.'
    print bar

#######################################

if __name__=='__main__':
    from giant.jiffies import run_default
    run_default(
        run                 = run,
        master_phil         = master_phil,
        args                = sys.argv[1:],
        blank_arg_prepend   = blank_arg_prepend,
        program             = PROGRAM,
        description         = DESCRIPTION)

