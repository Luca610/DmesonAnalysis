import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import yaml
from FitterUtils import fitInvMass

def LoadSparse(fileName, dirName, sparseName):

    print('Loading THnSparses from file', fileName)
    infileData = ROOT.TFile(fileName)
    if not infileData:
        print(f'File {fileName} not found!')
        return None
    indirData = infileData.Get(dirName)
    if not indirData:
        print(f'Directory {dirName} not found!')
        return None
    sparse = indirData.Get(sparseName)
    if not sparse:
        print(f'ERROR: sparse {sparseName} not found!')
        return None
    return sparse

def ProjectSparse(config):
    '''
    Project ThnSparse in pT and bkgScore bins and save TH1s in another file
    '''
    # configuration
    with open(config, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    inputFile = inputCfg['inputfile']
    dirName = inputCfg['dirname']
    sparseName = inputCfg['sparsename']
    if not os.path.exists(os.path.dirname(inputCfg['projfile'])):
        os.makedirs(os.path.dirname(inputCfg['projfile']))
    outFile = ROOT.TFile(inputCfg['projfile'], 'recreate')
    cutVars = inputCfg['cutvars']
    assert len(cutVars['pt']['min']) == len(cutVars['pt']['max']) == len(cutVars['bkgBDT']['min']) == len(cutVars['bkgBDT']['max'])
    # load sparse and axis
    sparse = LoadSparse(inputFile, dirName, sparseName)
    ptAxis = sparse.GetAxis(cutVars['pt']['axisnum'])
    bkgBDTAxis = sparse.GetAxis(cutVars['bkgBDT']['axisnum'])
    # for loop for projection
    for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['pt']['min'], cutVars['pt']['max'])):
        print(f'Projecting pT range {ptMin:.1f}-{ptMax:.1f}')
        directory = outFile.mkdir(f"pT {ptMin:.1f}-{ptMax:.1f}")
        directory.cd()
        ptAxis.SetRangeUser(ptMin, ptMax)
        bkgBDTlims = np.linspace(cutVars['bkgBDT']['min'][iPt], cutVars['bkgBDT']['max'][iPt], inputCfg['nsteps'] + 1)
        for bkgBDTmax in bkgBDTlims[1:]:
            bkgBDTAxis.SetRangeUser(0, bkgBDTmax)
            hMass = sparse.Projection(0)
            hMass.SetName(f'hMass_pT{ptMin:.1f}-{ptMax:.1f}_bkgBDT{bkgBDTmax:.3f}')
            hMass.Write()
        outFile.cd()

def FitForOptimization(config):
    with open(config, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    inputFileName = inputCfg['projfile']
    cutVars = inputCfg['cutvars']
    outDir = inputCfg['outputfolder']
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    ROOT.gROOT.SetBatch(True)
    inputFile = ROOT.TFile(inputFileName, "read")
    for iPt, key in enumerate(inputFile.GetListOfKeys()):  #Loop on folders (1 for pT Bin)
        if (key.ReadObj().InheritsFrom("TDirectory")):
            folder = key.ReadObj()
        else: 
            continue
        ptAvg = (cutVars['pt']['min'][iPt] + cutVars['pt']['max'][iPt])/ 2
        if ptAvg > 10:
            ptAvg = 10
        c = ROOT.TCanvas("c", f"Cut optimization for {cutVars['pt']['min'][iPt]}<pT<{cutVars['pt']['max'][iPt]}",1400, 800)
        c.Divide(5,4)
        bkgBDTlims = np.linspace(cutVars['bkgBDT']['min'][iPt], cutVars['bkgBDT']['max'][iPt], inputCfg['nsteps'] + 1)
        meanList = []
        sMeanList = []
        sigmaList = []
        sSigmaList = []
        significanceList = []
        nSigList = []
        sNSigList = []
        nBkgList = []
        sNBkgList = []
        for iBDT, subkey in enumerate(folder.GetListOfKeys()): # Loop on histograms (1 for BDT cut)
            if (subkey.ReadObj().InheritsFrom("TH1")):
                hist = subkey.ReadObj()
                if inputCfg['rebin'][iPt] > 1:
                    hist.Rebin(inputCfg['rebin'][iPt])
            else: 
                continue
            
            workspace, massFrame = fitInvMass(hist, inputCfg['mass']['min'],  inputCfg['mass']['max'], cutVars['massFitRange']['min'][iPt],  cutVars['massFitRange']['max'][iPt], inputCfg['pdgId'], ptAvg)  
            # retrive fit results from workspace
            mass = workspace.var("mass")
            nBkg = workspace.var("nBkg").getVal()
            nSig = workspace.var("nSig").getVal()
            mean = workspace.var("mean").getVal()
            sigma = workspace.var("sigma").getVal()
            totPDF = workspace.pdf("totPdf")
            signalPDF = workspace.pdf("signalPdf")
            bkgPDF = workspace.pdf("bkgPdf")
            data = workspace.data("data")
            # ndf = data.numEntries() - totPDF.getParameters(data).getSize()
            ndf =  hist.FindBin(cutVars['massFitRange']['max'][iPt]) - hist.FindBin(cutVars['massFitRange']['min'][iPt]) - totPDF.getParameters(data).getSize()
            chi2 = massFrame.chiSquare()
            # mass.setRange("fitRange", cutVars['massFitRange']['min'][iPt], cutVars['massFitRange']['max'][iPt])
            # chi2 = totPDF.createChi2(data, ROOT.RooFit.Range("fitRange")).getVal()
            # compute bkg 3 sigma and significance
            mass.setRange("intRange", mean - 3*sigma , mean + 3*sigma)
            int3Sigma = bkgPDF.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)), ROOT.RooFit.Range("intRange"))
            bkg3Sigma = nBkg*int3Sigma.getVal()
            significance = nSig / (np.sqrt(nSig + bkg3Sigma))
            meanList.append(mean)
            sMeanList.append(workspace.var("mean").getError())
            sigmaList.append(sigma)
            sSigmaList.append(workspace.var("sigma").getError())
            significanceList.append(significance)
            nSigList.append(nSig)
            sNSigList.append(workspace.var("nSig").getError())
            nBkgList.append(bkg3Sigma)
            sNBkgList.append(np.sqrt(bkg3Sigma))
            # Plotting options
            legend = ROOT.TLegend(0.4, 0.5, 0.8, 0.9)
            legend.SetTextSize(0.06)  # Adjust text size if needed
            legend.SetBorderSize(0)
            legend.SetFillStyle(0) 
            legend.AddEntry("","Fit results:", "")
            legend.AddEntry("",f"chi2/ndf = {(chi2*ndf):.1f}/{ndf:.0f}", "")
            # legend.AddEntry("",f"mean = {(mean*1000):.1f} MeV", "")
            # legend.AddEntry("",f"sigma = {(sigma*1000):.1f} MeV", "")
            legend.AddEntry("",f"nSig = {nSig:.0f}", "")
            # legend.AddEntry("",f"nBkg = {nBkg:.0f}", "")
            legend.AddEntry("",f"nBkg (3Sigma) = {bkg3Sigma:.0f}", "")
            legend.AddEntry("",f"significance (3Sigma) = {significance:.1f}", "")
            c.cd(iBDT + 1)
            pad = ROOT.gPad
            pad.SetLeftMargin(0.1)   
            pad.SetRightMargin(0.1)
            pad.SetTopMargin(0.1)
            pad.SetBottomMargin(0.1)
            massFrame.SetTitle(f"BKGscore < {bkgBDTlims[iBDT +1]:.3f}")
            massFrame.SetTitleSize(0.1)  # X-axis title size
            massFrame.GetXaxis().SetTitleSize(0.06)  # X-axis title size
            massFrame.GetYaxis().SetTitleSize(0.06)  # Y-axis title size
            massFrame.GetXaxis().SetLabelSize(0.06)  # X-axis label size
            massFrame.GetYaxis().SetLabelSize(0.06)  # Y-axis label size
            massFrame.DrawClone()
            legend.DrawClone("same")
        c.SaveAs(f"{outDir}/CutOpt_{cutVars['pt']['min'][iPt]:.0f}_{cutVars['pt']['max'][iPt]:.0f}.png")
        fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(16,9))
        plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.3, wspace=0.3)
        ax[0,0].errorbar(bkgBDTlims[1:], meanList,yerr=sMeanList, fmt='o', label='mean', markersize=5, color='black' )
        ax[0,0].set_xlabel("bkgBDT upper limit")
        ax[0,0].set_ylabel("mean (GeV)")
        ax[0,1].errorbar(bkgBDTlims[1:], sigmaList,yerr=sSigmaList, fmt='o', label='sigma', markersize=5, color='black' )
        ax[0,1].set_xlabel("bkgBDT upper limit")
        ax[0,1].set_ylabel("sigma (GeV)")
        ax[1,0].scatter(bkgBDTlims[1:], significanceList, label='significance', color='black' )
        ax[1,0].set_xlabel("bkgBDT upper limit")
        ax[1,0].set_ylabel("significance")
        ax[1,1].errorbar(bkgBDTlims[1:], nSigList,yerr=sNSigList, fmt='o', label='raw yield', markersize=5, color='black' )
        ax[1,1].set_xlabel("bkgBDT upper limit")
        ax[1,1].set_ylabel("RawYield")
        if inputCfg['pdgId'] == 411:
            ax[0,0].set_ylim(1.83, 1.90)
            ax[0,1].set_ylim(0, 0.1)
            ax[1,1].set_ylim(0, None)

        plt.legend()
        fig.savefig(f"{outDir}/summary_pt_{cutVars['pt']['min'][iPt]:.0f}_{cutVars['pt']['max'][iPt]:.0f}.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    parser.add_argument("--proj", help="perform only Sparse projection", action="store_true")
    parser.add_argument("--opt", help="perform only Optimization", action="store_true")
    args = parser.parse_args()

    if not args.opt:
        print("Start THn projection")
        ProjectSparse(args.cfgFileName)
        print("THn projection done!")
    
    if not args.proj:
        print("Start Fitting")
        FitForOptimization(args.cfgFileName)
        print("Fitting Done!")

    