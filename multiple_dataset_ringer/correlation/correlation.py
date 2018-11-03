import os
import pandas as pd
import numpy as np

def correlation_single_residue(input_csv,
                               output_dir,
                               out_filename):

    """Generate correlation matrix as csv"""

    # read in csv 
    data = pd.read_csv(input_csv, index_col=0)

    dataset_labels=data.index.values
    all_map_values=data.values

    # Generate correlation coefficents
    correlation_matrix = np.corrcoef(all_map_values)
    # Store in labelled data frame
    correlation_data = pd.DataFrame(correlation_matrix,
                                        index = dataset_labels, 
                                        columns = dataset_labels)
    # Correlation data as CSV 
    correlation_data.to_csv(os.path.join(output_dir,out_filename))
    

