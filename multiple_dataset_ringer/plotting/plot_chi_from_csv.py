import os
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

def plot_ringer_dataset(csv, datasets, plot_path, plot_all=True, plot_median=True, colours=None):
    """

    Parameters
    ----------
    csv: str
        path to csv file from ringer vd chi angle
    dataset: str, list
        name of datset to highlight, or list thereor
    plot_all: bool
        whether to plot all datasets in background
    plot_path:

    Returns
    -------
    None
    """
    ringer_df = pd.read_csv(csv)
    ringer_df = ringer_df.set_index(keys='Unnamed: 0')

    angles = ringer_df.columns.values
    median = ringer_df.median(axis=0)

    # Convert str to list if needed
    datasets = [datasets] if isinstance(datasets, str) else datasets

    if colours is None:
        raise ValueError("Need Colours")

    ax = plt.subplot(111)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    for index, ringer_result in ringer_df.iterrows():

        if plot_all:
            plt.plot(angles,
                     ringer_result.values,
                     color='k',
                     alpha=0.05,
                     label='_nolegend_',
                     )
    for index, ringer_result in ringer_df.iterrows():

        for i, dataset in enumerate(datasets):
            if dataset == index:

                plt.plot(angles,
                         ringer_result.values,
                         label=dataset,
                         color=colours[i],
                         linewidth=2,
                         )
    if plot_median:
        plt.plot(angles,
                 median,
                 color='y',
                 label="median",
                 linewidth=2,
                 )

    plt.xlim(0,360)
    plt.xlabel("Angle")
    plt.ylabel("Electron Density")
    plt.legend()
    plt.savefig(plot_path,dpi=300)
    plt.close()

if __name__ == "__main__":
    folder = "/dls/science/groups/i04-1/elliot-dev/ringer_tcru"
    csv_ext = "_616_Datasets_2mFo-DFc_chi1-ringer.csv"
    base_dir_plots = "/dls/science/groups/i04-1/elliot-dev/ringer_tcru_plots"

    datasets = ["TCRUFPPS-x0051","TCRUFPPS-x0106"]

    if not os.path.exists(base_dir_plots):
        os.makedirs(base_dir_plots)

    # residue : (folder, csv_path)
    residue_folders = {residue_folder:
                           (os.path.join(folder,residue_folder),
                            os.path.join(folder,residue_folder,
                                         "{residue_folder}{csv_ext}".format(residue_folder=residue_folder,
                                                                            csv_ext=csv_ext)))
                         for residue_folder
                         in os.listdir(folder)
                         if os.path.isdir(os.path.join(folder,residue_folder))
                       }

    #for dataset in datasets:
    for residue, (folder, csv_path) in residue_folders.iteritems():

        if not os.path.exists(os.path.join(base_dir_plots,"TCRUFPPS-x0051_106")):
            os.makedirs(os.path.join(base_dir_plots,"TCRUFPPS-x0051_106"))

        plot_ringer_dataset(csv=csv_path,
                            datasets=["TCRUFPPS-x0051",
                                      "TCRUFPPS-x0106"],
                            plot_all=True,
                            colours=['m','c'],
                            plot_path=os.path.join(base_dir_plots,
                                                   "TCRUFPPS-x0051_106",
                                                   "{residue}.png".format(residue=residue)))


