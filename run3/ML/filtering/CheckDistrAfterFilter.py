import pandas as pd
import sys
from hipe4ml import plot_utils
import matplotlib.pyplot as plt

sys.path.append('/home/luca/alice/DmesonAnalysis/utils')
from DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error

dfPaths = [#'../../../analyses/hf_run3/flowD/input/LHC23/Data_Dpluspp13TeV_run3_LHC22b1b_pT_1_50.parquet.gzip',
           '/home/luca/HF_data/TrainingBDTapass7/filtered/Prompt_Dstar_pprun3_LHC24d3_pT_0_100.parquet.gzip',
           '/home/luca/HF_data/TrainingBDTapass7/filtered/FD_Dstar_pprun3_LHC24d3_pT_0_100.parquet.gzip',
           '/home/luca/HF_data/TrainingBDTapass7/filtered/Data_Dstar_pprun3_LHC22pass7_pT_0_100.parquet.gzip']
           
dfDirs = ['', '', '']
dfTrees = ['O2hfcanddstlite', 'O2hfcanddstlite', 'O2hfcanddstlite']
dfLabels = ['Bkg', 'Prompt', 'FD']
colsToPlot = ['fPt', 'fDeltaM']
presels = 'fPt > 0 and fPt < 24'
listOfKeys = []
Dfs = []
for i, (df, dirName, treeName, label) in enumerate(zip(dfPaths, dfDirs, dfTrees, dfLabels)):
    print(f'Loading {df}')
    data = LoadDfFromRootOrParquet(df, dirName, treeName)
    data.query(presels, inplace=True)
    if i == 0:
        listOfKeys = data.columns
    Dfs.append(data[colsToPlot])

plot_utils.plot_distr(Dfs, colsToPlot, 1000, dfLabels, figsize=(20, 20), density=True, histtype='stepfilled', grid=False, log=True, alpha=0.3)
plt.savefig('checkDistr.pdf')
print('Done')
    