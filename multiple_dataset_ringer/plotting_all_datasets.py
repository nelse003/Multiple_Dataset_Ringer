import os
import sys
from multiple_dataset_ringer import master_phil
from multiple_dataset_ringer import blank_arg_prepend
from process.process import process_all_with_ringer
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import re
from collections import OrderedDict
from collections import defaultdict
from plotting.plots import multiple_line_plot_ringer
from iotbx.pdb import hierarchy
import numpy as np

def read_occ_b_factor(pdb_path, chain):
    """Extract occupancy and B factor of chain of 
    interest from one PDB file into a dataframe"""

    # Read in single PDB file
    pdb_in = hierarchy.input(file_name=pdb_path)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()
    lig_sel = sel_cache.selection("chain {}".format(chain))
    lig_hierarchy = pdb_in.hierarchy.select(lig_sel)

    lig_occ_b = []
    # Get occupancy & B factor of ligand
    for model in lig_hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                rg_occ = []
                rg_b = []

                for ag in rg.atom_groups():
                    for atom in ag.atoms():
                         rg_occ.append(atom.occ)
                         rg_b.append(atom.b)

                lig_occ_b.append([chain.id,
                                  rg.resseq,
                                  ag.resname,
                                  np.mean(rg_occ),
                                  np.mean(atom.b)])

    occ_b_df = pd.DataFrame(lig_occ_b,
                            columns=["Chain", "Resid", "Residue",
                                     "Occupancy", "B_factor"])

    return occ_b_df


def run(params):
    """
    Quick way to get all datasets from existing run to test plotting runs & metrics
    
    Example
    ------------
    
    ccp4-python plotting_all_datasets.py 
    /hdlocal/home/enelson/PTP1B/datasets_aligned/* 
    input.pdb_style="*pandda-input.pdb" 
    input.mtz_style="*_mrflagsref_idxs.mtz" 
    output.out_dir="/hdlocal/home/enelson/PTP1B/output_aligned"
    
    """
    all_results_pickle_path = "/hdlocal/home/enelson/PTP1B/datasets_aligned/all_results.pickle"

    if not os.path.exists(all_results_pickle_path):
        all_results = process_all_with_ringer(params)
        with open(all_results_pickle_path,'w') as all_results_pickle:
            pickle.dump(all_results, all_results_pickle)

    with open(all_results_pickle_path, 'r') as all_results_pickle:
        all_results = pickle.load(all_results_pickle)
    all_results_df = pd.concat(all_results)
    ref_set = all_results.values()[0]
    ref_set_chi1 = ref_set.loc[(ref_set[2] == params.settings.angle_type)]

    datasets = [os.path.basename(dataset)
                for dataset in params.input.dir
                if os.path.isdir(os.path.join(params.input.dir, dataset))]

    #base_dir = os.path.dirname(params.input.dir[0])
    # needs to pulll from non aligned to get variation in b factor,
    # although this may not be fully refined

    base_dir = "/hdlocal/home/enelson/PTP1B/datasets"
    b_occ_all = []
    for dataset in datasets:
        print(dataset)
        pdb_path = os.path.join(base_dir,
                                dataset,
                                "{}-{}".format(dataset, params.input.pdb_style.lstrip('*')))
        b_occ_all.append(read_occ_b_factor(pdb_path, 'A'))

    b_occ_all_df = pd.concat(b_occ_all, keys=datasets)

    b_factor_mean = b_occ_all_df.groupby(['Resid']).agg(np.mean)['B_factor'].values
    b_factor_std = b_occ_all_df.groupby(['Resid']).agg(np.std)['B_factor'].values
    res = b_occ_all_df.groupby(['Resid']).agg(np.mean).index.values.astype(int)

    plt.errorbar(x=res, y=b_factor_mean, yerr=b_factor_std)
    plt.xlabel('Residue')
    plt.ylabel('B Factor')
    plt.savefig('bfactor_errorbar.png')
    plt.close()

    interpolate_base_csv = \
        '_{}_Datasets_{}_{}-ringer.csv'.format(len(datasets),
                                               params.settings.map_type,
                                               params.settings.angle_type)
    interpolate_results = {}

    std_mad = OrderedDict()
    mean_mad = OrderedDict()
    norm_mean_mad = OrderedDict()
    norm_results = OrderedDict()
    norm_median_mad = OrderedDict()
    norm_std_mad = OrderedDict()
    norm_median_mad_std_mad = OrderedDict()
    norm_median_mean_mad = OrderedDict()

    for residue, data in ref_set_chi1.iterrows():
        print(residue)
        interpolate_csv = residue + interpolate_base_csv

        single_residue_multiple_datasets = pd.read_csv(
            os.path.join(params.output.out_dir,
                         residue,
                         interpolate_csv),index_col=0)

        interpolate_results[residue] = single_residue_multiple_datasets

        norm_results[residue] = (single_residue_multiple_datasets -
                single_residue_multiple_datasets.min().min()) / \
               (single_residue_multiple_datasets.max().max() -
                single_residue_multiple_datasets.min().min())

        # multiple_line_plot_ringer(results_df=norm_results[residue],
        #                           title=residue,
        #                           filename='all-{}-{}-norm_feature_scaled_dataset.png'.format(
        #                               residue, len(ref_set)),
        #                           out_dir=os.path.join(params.output.out_dir, residue))

        norm_median_mad[residue] = (single_residue_multiple_datasets -
                single_residue_multiple_datasets.median().mean()) / \
                single_residue_multiple_datasets.mad().mean()

        # multiple_line_plot_ringer(results_df=norm_results[residue],
        #                           title=residue,
        #                           filename='all-{}-{}-norm_median_mad_dataset.png'.format(
        #                               residue, len(ref_set)),
        #                           out_dir=os.path.join(params.output.out_dir, residue))

        std_mad[residue] = single_residue_multiple_datasets.mad().std()
        mean_mad[residue] = single_residue_multiple_datasets.mad().mean()
        norm_std_mad[residue] = norm_results[residue].mad().std()
        norm_median_mad_std_mad[residue] = norm_median_mad[residue].mad().std()
        norm_mean_mad[residue] = norm_results[residue].mad().mean()
        norm_median_mean_mad[residue] = norm_median_mad[residue].mad().mean()


    residues = list(std_mad.keys())
    residue_nums = [int(re.findall(r'\d+', residue)[0]) for residue in residues]

    b_dict = b_occ_all_df.groupby(['Resid']).agg(np.mean)['B_factor'].to_dict()
    b_ordered_dict = OrderedDict(sorted(b_dict.iteritems(), key=lambda x: x[0]))

    dd = defaultdict(list)

    for d in (b_ordered_dict, std_mad):
        for key, value in d.iteritems():
            key_short = int(re.findall(r'\d+', key)[0])
            dd[key_short].append(value)

    # Create a temporary copy of dictionary
    copy_dd = dict(dd)

    # Iterate over the temporary dictionary and delete corresponding key from original dictionary
    for (key, value) in copy_dd.items():
        if len(value) != 2:
            del dd[key]

    ################################
    # Normalised version

    dd = defaultdict(list)

    for d in (b_ordered_dict, norm_std_mad):
        for key, value in d.iteritems():
            key_short = int(re.findall(r'\d+', key)[0])
            dd[key_short].append(value)

    # Create a temporary copy of dictionary
    copy_dd = dict(dd)

    # Iterate over the temporary dictionary and delete corresponding key from original dictionary
    for (key, value) in copy_dd.items():
        if len(value) != 2:
            del dd[key]
    norm_std_mad_b = pd.DataFrame.from_dict(dd, orient='index')
    ########################################

    std_mad_b = pd.DataFrame.from_dict(dd, orient='index')

    # plt.text(0.1, 0.5,
    #         'Standard deviation across angles of '
    #         'median absolute deviation across datasets',
    #         horizontalalignment='left',
    #         verticalalignment='center',
    #         rotation='horizontal',
    #         transform=ax.transAxes)

    plt.scatter(x=std_mad_b[0], y=std_mad_b[1])
    plt.xlabel('B Factor')
    # plt.ylabel('Standard deviation across angles of'
    #            'median absolute deviation across datasets')
    plt.savefig('b_scatter_std_mad.png')
    plt.close()

    std_mad_div_b = std_mad_b[1]/std_mad_b[0]
    plt.scatter(x=dd.keys(), y=std_mad_div_b)
    plt.ylim(0,0.001)
    plt.ylabel('Standard deviation across angles of '
               'median absolute deviation across datasets / B factor')
    plt.xlabel('Residues (with sampled Chi1 angle)')
    plt.savefig('norm_std_mad_div_b.png')
    plt.close()

    norm_std_mad_b = pd.DataFrame.from_dict(dd, orient='index')
    plt.scatter(x=norm_std_mad_b[0], y=norm_std_mad_b[1])
    plt.xlabel('B Factor')
    plt.ylabel('Standard deviation across angles of'
               'median absolute deviation across'
               'feature normalised datasets')
    plt.savefig('b_scatter_norm_std_mad.png')
    plt.close()

    norm_std_mad_div_b = norm_std_mad_b[1]/norm_std_mad_b[0]
    plt.scatter(x=dd.keys(), y=norm_std_mad_div_b)
    plt.ylim(0,0.001)
    plt.ylabel('Standard deviation across angles of '
               'median absolute deviation across '
               'feature normalised datasets / B factor',fontsize=12)
    plt.xlabel('Residues (with sampled Chi1 angle)')
    plt.savefig('norm_std_mad_div_b.png')
    plt.close()

    plt.scatter(list(std_mad.values()), std_mad_b[0])
    plt.title("Correlation fo B factor with variability")
    plt.xlabel('Standard deviation across angles of '
               'median absolute deviation across datasets.',
               fontsize=12)
    plt.ylabel('Mean B factor across datasets')
    plt.savefig('b_factor_std_corr.png')
    plt.close()

    plt.scatter(residue_nums, list(std_mad.values()))
    plt.title('Variability metrics for ringer plots in PTP1B')
    plt.xlabel('Residue')
    plt.ylabel('Standard deviation across angles of '
               'median absolute deviation across datasets.',fontsize=12)
    plt.savefig('std_mad.png', dpi=300)
    plt.close()

    plt.scatter(residue_nums, list(mean_mad.values()))
    plt.title('Variability metrics for ringer plots in PTP1B')
    plt.xlabel('Residue')
    plt.ylabel('Mean across angles of mean absolute deviation '
               'across datasets.', fontsize=12)
    plt.savefig('mean_mad.png', dpi=300)
    plt.close()

    plt.scatter(residue_nums, list(norm_std_mad.values()))
    plt.title('Variability metrics for ringer plots in PTP1B')
    plt.xlabel('Residue')
    plt.ylabel('Standard deviation across angles of median '
               'absolute deviation across '
               'feature normalised datasets', fontsize=12)
    plt.savefig('norm_std_mad.png', dpi=300)
    plt.close()

    plt.scatter(residue_nums, list(norm_median_mad_std_mad.values()))
    plt.xlabel('Residue')
    plt.ylabel('Standard deviation across angles of '
               'median absolute deviation across datasets'
               ' normalised using medain scaling', fontsize=12)
    plt.savefig('norm_median_mad_std_mad.png')
    plt.close()

    plt.scatter(residue_nums, list(norm_mean_mad.values()))
    plt.xlabel('Residue')
    plt.ylabel('Mean across angles of '
               'median absolute deviation across datasets'
               'using feature normalised datasets', fontsize=12)
    plt.savefig('norm_mean_mad.png')
    plt.close()

    plt.scatter(residue_nums, list(norm_median_mean_mad.values()))
    plt.xlabel('Residue')
    plt.ylabel('Mean across angles of '
               'median absolute deviation across datasets'
               ' normalised using medain scaling', fontsize=12)
    plt.savefig('norm_median_mad_mean.png')
    plt.close()


    from IPython import embed; embed()

if __name__ == '__main__':
    from giant.jiffies import run_default

    run_default(run=run, master_phil=master_phil, args=sys.argv[1:], blank_arg_prepend=blank_arg_prepend)