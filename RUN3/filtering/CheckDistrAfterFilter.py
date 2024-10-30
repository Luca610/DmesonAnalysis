import pandas as pd
import sys
from hipe4ml import plot_utils
import matplotlib.pyplot as plt

sys.path.append('../../utils')
from DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error

dfPaths = [#'../../../analyses/hf_run3/flowD/input/LHC23/Data_Dpluspp13TeV_run3_LHC22b1b_pT_1_50.parquet.gzip',
           '/data/shared/DsReso_run3/data/Data_Dpluspp13TeV_run3_LHC22pass6_pT_1_50.parquet.gzip',
           '/data/shared/DsReso_run3/mc/Prompt_Dpluspp13TeV_run3_LHC22b1b_pT_1_50.parquet.gzip',
           '/data/shared/DsReso_run3/mc/FD_Dpluspp13TeV_run3_LHC22b1a_pT_1_50.parquet.gzip']
           
dfDirs = ['DF_2261906078563584', 'DF_2849422059001', 'DF_2847588306001']
dfTrees = ['O2hfcanddplite', 'O2hfcanddplite', 'O2hfcanddplite']
dfLabels = ['Bkg', 'Prompt', 'FD']
colsToPlot = ['fChi2PCA', 'fDecayLength', 'fDecayLengthXY', 'fDecayLengthNormalised',
       'fDecayLengthXYNormalised', 'fPtProng0', 'fPtProng1', 'fPtProng2',
       'fImpactParameter0', 'fImpactParameter1', 'fImpactParameter2',
       'fNSigTpcPi0', 'fNSigTpcKa0',
       'fNSigTpcPi1', 'fNSigTpcKa1',
       'fNSigTpcPi2', 'fNSigTpcKa2',
       'fM', 'fPt', 'fCpa', 'fCpaXY',
       'fMaxNormalisedDeltaIP']
presels = 'fPt > 1 and fPt < 50'

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
    