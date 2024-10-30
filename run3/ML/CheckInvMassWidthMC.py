import os
import sys
sys.path.append('/home/luca/alice/DmesonAnalysis/utils')
import argparse
import pickle
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm

from DfUtils import LoadDfFromRootOrParquet
from ROOT import TFile
from hipe4ml.tree_handler import TreeHandler

def sigma_par(pT):
    return (0.00796 + 0.00176*pT) 

def main(): #pylint: disable=too-many-statements
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    PromptHandler = TreeHandler(inputCfg['input']['prompt'], inputCfg['input']['treename'])

    PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]
    PromptHandler.slice_data_frame('fPt', PtBins, True)
    
    fig, ax = plt.subplots(nrows = 4, ncols = 3, figsize=(16,9))
    muList = [] 
    sigmaList = [] 
    meanPtList = []
    pdgMass = 1.87
    for iBin, PtBin in enumerate(PtBins):
        row = iBin//3
        col = iBin%3
        meanPtList.append((PtBin[0] + PtBin[1])/2 )
        df = PromptHandler.get_slice(iBin)
        mass = df['fM']
        count, bins, ignored = ax[row, col].hist(mass, bins=50, density=True, alpha=0.6, color='g')
        fitMin = pdgMass - 4*sigma_par(meanPtList[iBin])
        fitMax = pdgMass + 4*sigma_par(meanPtList[iBin])
        filtered_data = mass[(mass >= fitMin) & (mass <= fitMax)]
        mu, std = norm.fit(filtered_data)
        muList.append(mu) 
        sigmaList.append(std)
        xmin, xmax = ax[row, col].get_xlim()  # Ottiene i limiti dell'asse x
        x = np.linspace(fitMin, fitMax, 100)
        p = norm.pdf(x, mu, std)
        ax[row, col].plot(x, p, 'k', linewidth=2)
    plt.scatter(meanPtList, sigmaList)
    fig.savefig("WifthDplusPeakTrkTuner1p5.png")
    
    print(muList, sigmaList)
main()