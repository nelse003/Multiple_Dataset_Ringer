import pandas
import os

###############################################################################
# Logging
###############################################################################
import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh = logging.FileHandler('ringer_script.log')
fh.setLevel(logging.DEBUG)
fh.setFormatter(formatter)
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)
###############################################################################

def peak_to_array(out_dir,conformer_dir,peak_filename,conformer_base_filename): 
        
    if not os.path.exists(os.path.join(out_dir,conformer_dir)):
        os.mkdir(os.path.join(out_dir,conformer_dir))

    all_peaks = pandas.read_csv(os.path.join(out_dir,peak_filename),
                                    header = 0, index_col =0)

    # Split into separate dataframe (in dictionary with dataset name as index)
    PeaksDataFrameDict = {elem : pandas.DataFrame for elem in all_peaks.index.unique()}
    # Dictiorary of arrays to store results from 
    PeaksArrayDataFrameDict ={elem: pandas.DataFrame for elem in all_peaks.index.unique()}

    dataset_count=1

    for key in PeaksDataFrameDict.keys():
        if not os.path.exists(os.path.join(out_dir,conformer_dir,key+
                                           conformer_base_filename)):

            PeaksDataFrameDict[key] = all_peaks[:][all_peaks.index == key]

            currentDF = PeaksDataFrameDict[key]
            currentDF = currentDF.set_index(['Residue'])
            
            ArrayDF = pandas.DataFrame(
                      index = PeaksDataFrameDict[key]['Residue'].values, 
                      columns = PeaksDataFrameDict[key]['Residue'].values)

            for residue_i in currentDF.index.values:
               for residue_j in currentDF.index.values:

                   # For beta branched residues (ILE,VAL,THR), alternate 
                   # conformations only exist if more than two peaks shown

                   if ("ILE" or "VAL" or "THR") in residue_i:
                       NumAltConf_i = currentDF.at[residue_i,'Number of peaks']-2
                   else:
                       NumAltConf_i = currentDF.at[residue_i,'Number of peaks']-1

                   if ("ILE" or "VAL" or "THR") in residue_j:
                       NumAltConf_i = currentDF.at[residue_j,'Number of peaks']-2
                   else:
                       NumAltConf_j = currentDF.at[residue_j,'Number of peaks']-1

                   # Convert to 1 for alternate conformers present, 
                   # 0 for no alternate conformers                   

                   if NumAltConf_i*NumAltConf_j > 0:
                       ArrayDF.at[residue_i,residue_j] = 1
                   else:
                       ArrayDF.at[residue_i,residue_j] = 0     

            PeaksArrayDataFrameDict[key] = ArrayDF
            # Output array to CSV
            ArrayDF.to_csv(os.path.join(out_dir,conformer_dir,
                           key+conformer_base_filename))

            logger.info('Converting Peak CSV to Array for Dataset ' + 
                        str(dataset_count)
                        + ' of '+ str(len(PeaksArrayDataFrameDict)))
        else:
            logger.info('Conformer Array already generated for dataset ' +
                        str(dataset_count))
    
        dataset_count+=1
