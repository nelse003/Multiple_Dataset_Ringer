import libtbx.phil

blank_arg_prepend = {None: "dir=", ".pdb": "pdb=", ".mtz": "mtz="}
master_phil = libtbx.phil.parse(
    """
input {
    dir = None
        .type = path
        .multiple = True
    pdb_style = "dimple.pdb"
        .type = str
        .multiple = False
    mtz_style = "dimple.mtz"
        .type = str
        .multiple = False
    column_labels = "FWT,PHWT"
        .type = str
        .multiple = False
        .help = mtz column labels of input files
}
output {
    log = "ringer.log"
        .type = str
        .multiple = False
    out_dir = "output"
        .type = str
        .multiple = False
    tmp_dir = "tmp"
        .type = str
        .multiple = False
}
settings {
    # XXX mmtbx.ringer can only take this an integer, >1 XXX#
    angle_sampling = 2
        .type = int
        .multiple = False

    gen_heatmap = False
        .type = bool
        .multiple = False
        
    map_type = '2mFo-DFc'
        .type = str
        .help = Name of electron density map to be analysed with ringer
        
    angle_type = 'chi1'
        .type = str
        .help = Chi angle to be analysed across multiple datasets
        
    sample_only_ref_pdb = None
        .type = str
        .help = Run all ringer analyses against a single pdb structure
    qsub = False
        .type = bool
        .help = flag to run (ringer) with qsub
}
"""
)
