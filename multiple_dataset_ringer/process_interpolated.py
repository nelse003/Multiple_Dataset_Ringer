import pandas as pd
from scipy.signal import find_peaks
from matplotlib import pyplot as plt

interpolated_df = pd.read_csv("/home/nelse003/PycharmProjects/Multiple_Dataset_Ringer/output/interpolated_datasets.csv", index_col=[0,1])

print(interpolated_df.index[0])

x = interpolated_df.loc[('ARG A105', 'PTP1B-y0004')])
find_peaks(x, height=0.3)