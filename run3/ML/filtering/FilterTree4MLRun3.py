'''
python script to filter tree from task output and save output trees in parquet files for ML studies
run: python FilterTree4ML.py cfgFileName.yml
'''

import sys
import argparse
import numpy as np
import yaml
sys.path.append('/home/luca/alice/DmesonAnalysis/utils')
from DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error
import pandas as pd
import uproot
import os

# fFlagMcMatchRec
flagSignal = {'D0': 0, 'Dplus': 1, 'Dstar':1, 'Ds': 4} #bitwise
# fOriginMcRec
bitPrompt = 1
bitFD = 2

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='cfgFileName.yml',
                    help='input config yaml file name')

args = parser.parse_args()
print('Opening input file')
with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

channel = cfg['channel']
if channel not in ['Ds', 'D0', 'Dplus', 'Dstar']:
    print('Error: only Ds, D0, Dplus, Dstar channels are implemented! Exit')
    sys.exit()

inFileNames = cfg['infile']['filename']
inDirName = cfg['infile']['dirname']
inTreeName = cfg['infile']['treename']
isMC = cfg['infile']['isMC']
outTreeName = cfg['outfile']['treename']
outSuffix = cfg['outfile']['suffix']
outDirName = cfg['outfile']['dirpath']
preSelections = cfg['skimming']['preselections']
colsToKeep = cfg['skimming']['colstokeep']
PtMins = cfg['skimming']['pt']['min']
PtMaxs = cfg['skimming']['pt']['max']

if not os.path.exists(outDirName):
    print(f"Creating output directory: {outDirName}")
    os.makedirs(outDirName)
dataframe = pd.DataFrame()
if not isinstance(inFileNames, list):
        inFileNames = [inFileNames]
if not inDirName:
    for file in inFileNames:
        f = uproot.open(file)
        dirnames = f.keys(recursive = False)
        print(f"keys found in file: {dirnames}")
        for dirname in dirnames:
            if 'parent' in dirname:
                print('got parentFiles')
                continue
            else:
                dataframe = pd.concat([dataframe, LoadDfFromRootOrParquet(file, dirname, inTreeName)], ignore_index=True)
                print(f"processed key: {dirname}")
else:  
    dataframe = LoadDfFromRootOrParquet(inFileNames, inDirName, inTreeName)
dataFramePtIntSel = pd.DataFrame()
dataFramePtIntSelSigPrompt = pd.DataFrame()
dataFramePtIntSelSigFD = pd.DataFrame()
dataFramePtIntSelBkg = pd.DataFrame()

print('Applying selections')
for PtMin, PtMax, presel in zip(PtMins, PtMaxs, preSelections):
    if presel:
        dataFramePtCutSel = dataframe.query(f'fPt > {PtMin} and fPt < {PtMax} and {presel}')
    else:
        dataFramePtCutSel = dataframe.query(f'fPt > {PtMin} and fPt < {PtMax}')
        
    if dataFramePtCutSel.empty:
        print(f'No data to save for pT range {PtMin:.0f} - {PtMax:.0f}')
        continue
    
    if channel == 'Dstar':
        dataFramePtCutSel['fDeltaM'] = dataFramePtCutSel['fM'] - dataFramePtCutSel['fMD0']  
    print(f'Number of events for pT range {PtMin:.0f} - {PtMax:.0f}: {dataFramePtCutSel.shape[0]}')
    print(dataFramePtCutSel[colsToKeep].head(5))

    bitsForSel, labelsContr = ({} for _ in range(2))
    if isMC:
        dataFramePtCutSelSigPrompt = dataFramePtCutSel.query(f'abs(fFlagMcMatchRec) == {flagSignal[channel]} and fOriginMcRec == {bitPrompt}', inplace=False)
        dataFramePtCutSelSigFD = dataFramePtCutSel.query(f'abs(fFlagMcMatchRec) == {flagSignal[channel]} and fOriginMcRec == {bitFD}', inplace=False)
        dataFramePtCutSelBkg = dataFramePtCutSel.query(f'abs(fFlagMcMatchRec) != {flagSignal[channel]}', inplace=False)
        
        if not dataFramePtCutSelSigPrompt.empty:
            print(f'Number of signal prompt events for pT range {PtMin:.0f} - {PtMax:.0f}: {dataFramePtCutSelSigPrompt.shape[0]}')
            dataFramePtIntSelSigPrompt = pd.concat([dataFramePtIntSelSigPrompt,dataFramePtCutSelSigPrompt], ignore_index=True) 
        else:
            print(f'No signal prompt data to save for pT range {PtMin:.0f} - {PtMax:.0f}')
        
        if not dataFramePtCutSelSigFD.empty:
            print(f'Number of signal FD events for pT range {PtMin:.0f} - {PtMax:.0f}: {dataFramePtCutSelSigFD.shape[0]}')
            dataFramePtIntSelSigFD = pd.concat([dataFramePtIntSelSigFD,dataFramePtCutSelSigFD], ignore_index=True) 

        else:
            print(f'No signal FD data to save for pT range {PtMin:.0f} - {PtMax:.0f}')
        
        if not dataFramePtCutSelBkg.empty:
            print(f'Number of background events for pT range {PtMin:.0f} - {PtMax:.0f}: {dataFramePtCutSelBkg.shape[0]}')
            dataFramePtIntSelBkg = pd.concat([dataFramePtIntSelBkg,dataFramePtCutSelBkg], ignore_index=True) 

        else:
            print(f'No background data to save for pT range {PtMin:.0f} - {PtMax:.0f}')

    else:
        print(f'Number of events for pT range {PtMin:.0f} - {PtMax:.0f}: {dataFramePtCutSel.shape[0]}')
        dataFramePtIntSel = pd.concat([dataFramePtIntSel,dataFramePtCutSel], ignore_index=True)
                                           
print('Saving data to parquet')
if isMC:
    if not dataFramePtIntSelSigPrompt.empty:
        dataFramePtIntSelSigPrompt[colsToKeep].to_parquet(
                f'{outDirName}/Prompt{outSuffix}_pT_{PtMins[0]:.0f}_{PtMaxs[-1]:.0f}.parquet.gzip',
                compression='gzip')
    if not dataFramePtIntSelSigFD.empty:
        dataFramePtIntSelSigFD[colsToKeep].to_parquet(
                f'{outDirName}/FD{outSuffix}_pT_{PtMins[0]:.0f}_{PtMaxs[-1]:.0f}.parquet.gzip',
                compression='gzip')
    if not dataFramePtIntSelBkg.empty:
        dataFramePtIntSelBkg[colsToKeep].to_parquet(
                f'{outDirName}/Bkg{outSuffix}_pT_{PtMins[0]:.0f}_{PtMaxs[-1]:.0f}.parquet.gzip',
                compression='gzip')
else:
    if not dataFramePtIntSel.empty:
        dataFramePtIntSel[colsToKeep].to_parquet(
                f'{outDirName}/Data{outSuffix}_pT_{PtMins[0]:.0f}_{PtMaxs[-1]:.0f}.parquet.gzip',
                compression='gzip')