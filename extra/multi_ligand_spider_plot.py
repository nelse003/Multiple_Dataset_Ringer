import giant.score_model as score_model
import itertools

#########################################

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
"""

#######################################
# Functions for combined model scoring
######################################

def get_site_list(inspect_csv, pdb1):
   #XXX Not used    
    """ Get site for a the PDb file's dataset given the pandaa_inspect-event.csv"""

    pandda_inspect_event = pandas.read_csv(inspect_csv)
    dataset_id = os.path.split(os.path.dirname(pdb1))[1]
    dataset_events = pandda_analyse_event.loc[pandda_analyse_event['dtag'] == dataset_id]
    site_list = dataset_events["site_idx"].tolist()

    return site_list  

def prepare_heirarchy(pdb1):

    # Extract Strucutre
    pdb_hierarchy = strip_pdb_to_input(pdb1, remove_ter=True, remove_end=True).hierarchy
    sanitise_hierarchy(pdb_hierarchy)

    return pdb_hierarchy

def check_multiple_ligands(pdb_hierarchy,params):

    """ Check whether pdb file contains multiple lignads"""

    # Residues to look for
    res_names = params.selection.res_names_list
    # Pull out residues to analyse
    if res_names: rg_for_analysis = [rg for rg in pdb_hierarchy.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
    else:         rg_for_analysis = pdb_hierarchy.residue_groups()

    if len(rg_for_analysis) >1 and res_names == ["LIG","UNL","DRG"]:
        return True
    else:
        return False  

def seperate_multiple_ligands_into_sites(pdb_hierarchy):
    
    """ For a given pdb hierarchy split ligands into sites """

    res_names = params.selection.res_names_list

    if res_names: rg_for_analysis = [rg for rg in pdb_hierarchy.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]    

    from IPython import embed();
        
    #for each pair of 

    itertools.combinations

    return site_dict

def combine_ligands(pdb_hierarchy,site_dict,params):
    
    """ Combine residues based on site information so that Edstats/ score_model can process as one block"""

    # TODO :Include site information so that only ligands in one site are combined.
    # This is required to deal with cases where one pdb would contain ligands 
    # from multiple different sites
    # 
    # Issue: Identifying ligand to site record:
    #
    # Solution Ideas:
    #
    # Spatial location of ligands within (say) 5A of one another

    # Extract residues to look for
    res_names = params.selection.res_names_list

    # Pull out residues to analyse
    if res_names: 
        rg_for_analysis = [rg for rg in pdb_hierarchy.residue_groups() if [n for n in rg.unique_resnames() if n in res_names]]
    else:
        # If ligands are not being analysed 
        return pdb_hierarchy

    # TODO loop over sites in site dict
    #set all residue id of ligands to the same id
    for rg in rg_for_analysis:
        rg.resseq = rg_for_analysis[0].resseq      

    return pdb_hierarchy
###########################################

def run():

    pdb_hierarchy = prepare_hierarchy(params.input.pdb1)

    ###########################################################################
    # TESTING combined ligands
    ###########################################################################
    if check_multiple_ligands(pdb_hierarchy,params):
        print "Multiple ligands per file exist"
        site_dict = seperate_multiple_ligands_into_sites(pdb_hierarchy)

        combined_pdb_hierarchy = combine_ligands(pdb_hierarchy,site_dict,params)
        combined_pdb = os.path.dirname(pdb1)+'/'+os.path.splitext(os.path.basename(pdb1))[0]+'_combined_ligs'+'.pdb'
        combined_pdb_hierarchy.write_pdb_file(combined_pdb)

        print "Generated merged pdb file {}".format(combined_pdb)

        print 'Scoring combined model...'
        data_table = score_model(   params = params,
                                pdb1   = combined_pdb,
                                mtz1   = params.input.mtz1,
                                pdb2   = params.input.pdb2,
                                mtz2   = params.input.mtz2,
                                label_prefix = params.input.label + 'combined'
                            )
        print '...Done'
        print bar

        ###################################################################
        # Output Images
        ###################################################################
        all_images = []
        print 'Generating combined Output Images...'
        for label, row in data_table.iterrows():
            print 'Making: {}...'.format(label)
            image_path = os.path.join(images_dir,'{}.png'.format(label))
            make_residue_radar_plot(path = image_path,
                                    data = row.to_frame().T,
                                    columns = columns,
                                    remove_blank_entries = params.plot.remove_blank_entries,
                                    print_axis_values    = params.plot.print_axis_values  )
            all_images.append(image_path)
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
                                               

