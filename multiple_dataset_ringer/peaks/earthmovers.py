import os
import pandas as pd
from scipy.stats import wasserstein_distance
from matplotlib import pyplot as plt

def wasserstein_distance_all_residues(folder, chi, out_csv, dataset=None, median=True):

    """
    Get wassertien distance compared to median for each residue

    Parameters
    ----------
    folder: str
        folder with ringer results
    chi: str
        chi1, chi2, chi3 string mathcing which chi angle
    out_csv: str
        output csv path

    Returns
    -------

    """

    df_list = []
    for residue_folder in os.listdir(folder):

        if not os.path.isdir(os.path.join(folder, residue_folder)):
            continue

        for file_name in os.listdir(os.path.join(folder,
                                                  residue_folder)):


            if ".csv" in file_name and chi in file_name:
                csv = os.path.join(folder,
                                   residue_folder,
                                   file_name)

                ringer_df = pd.read_csv(csv)
                ringer_df = ringer_df.set_index(keys='Unnamed: 0')

                if median:
                    ringer_compare = ringer_df.median(axis=0)
                elif dataset is not None:
                    print(ringer_df.iloc[dataset])
                    ringer_compare = ringer_df.iloc[dataset]

                earthmover_dict = {}
                for index, ringer_result in ringer_df.iterrows():
                    earthmover_dict[index] = wasserstein_distance(ringer_compare,
                                                                  ringer_result.values)

                earthmover_df = pd.DataFrame.from_dict(earthmover_dict,
                                                       orient='index',
                                                       columns=[residue_folder])

                df_list.append(earthmover_df)

    out_df = pd.concat(df_list, axis=1)
    out_df.to_csv(out_csv)

def plot_wasserstein_distance(csv):

    wasserstein_df = pd.read_csv(csv, index_col=0)

    for index, wasserstein in wasserstein_df.iterrows():
        if "TCRUFPPS-x0051" == index:
            plt.plot(wasserstein.values, color='r')
        if "TCRUFPPS-x0106" == index:
            plt.plot(wasserstein.values, color='b')
        else:
            plt.plot(wasserstein.values, color='k',alpha=0.01)
        # elif plot_all == True:
        #     plt.plot(angles, ringer_result.values, color='k', alpha=0.05)
    plt.show()


if __name__ == "__main__":

    """
    
    Notes
    ---------
    Can't use ccp4-python as scipy not up to date.
    """
    # wasserstein_csv = "/dls/science/groups/i04-1/elliot-dev/ringer_tcru/eartmovers_median.csv"
    #
    # if not os.path.exists(wasserstein_csv):
    #     wasserstein_distance_all_residues(folder="/dls/science/groups/i04-1/elliot-dev/ringer_tcru",
    #                                       chi="chi1",
    #                                       out_csv=wasserstein_csv)

    wasserstein_csv = "/dls/science/groups/i04-1/elliot-dev/ringer_tcru/eartmovers_tcrufpps_x0106.csv"
    if not os.path.exists(wasserstein_csv):
        wasserstein_distance_all_residues(folder="/dls/science/groups/i04-1/elliot-dev/ringer_tcru",
                                          chi="chi1",
                                          dataset="TCRUFPPS-x0106",
                                          out_csv=wasserstein_csv)
    #plot_wasserstein_distance(csv=wasserstein_csv)

    wasserstein_df = pd.read_csv(wasserstein_csv, index_col=0)
    # print(wasserstein_df.loc['TCRUFPPS-x0051'].nlargest(10))
    # print(wasserstein_df.loc['TCRUFPPS-x0106'].nlargest(10))

    print(wasserstein_df['PHE A50'].nsmallest(30))

    ratio_df = wasserstein_df.div(wasserstein_df.mean())
    # print(ratio_df.loc['TCRUFPPS-x0106'].nlargest(10))
    # print(ratio_df.loc['TCRUFPPS-x0051'].nlargest(10))