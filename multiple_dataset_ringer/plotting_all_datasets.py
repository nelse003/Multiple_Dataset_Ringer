import sys
from multiple_dataset_ringer import master_phil
from multiple_dataset_ringer import blank_arg_prepend
from process.process import process_all_with_ringer
import pickle

def run(params):
    """
    Quick way to get all datasets from existing run to test plotting runs
    
    Example
    ------------
    
    ccp4-python plotting_all_datasets.py 
    /hdlocal/home/enelson/PTP1B/datasets_single_pdb/* 
    input.pdb_style="*pandda-input.pdb" 
    input.mtz_style="*_mrflagsref_idxs.mtz" 
    output.out_dir="/hdlocal/home/enelson/PTP1B/output_single_pdb"
    settings.sample_only_ref_pdb="yes"
    
    """

    all_results = process_all_with_ringer(params)
    pickle.dump(all_results,"/hdlocal/home/enelson/PTP1B/datasets_single_pdb/all_results.pickle")

if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)