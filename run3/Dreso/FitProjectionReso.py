from ROOT import *
import sys
import argparse
import yaml
import os
from ResoFitUtils import fit_invariant_mass_with_background, roofit_plot_with_matplotlib, simultaneous_fit_signal_background, template_fit_signal_background
import numpy as np


'''
Simple script to fit THnSparse projection of charmed resonances (for performance plot)
Author: Luca Aglietta
'''

def getChi2Ndf(hist, f, xmin, xmax):
    chi2 = 0
    binLow = hist.FindBin(xmin)
    binHigh = hist.FindBin(xmax)
    for b in range(binLow, binHigh):
        x = hist.GetBinCenter(b)
        y = f.Eval(x)
        n = hist.GetBinContent(b)
        chi2 += (n - y)*(n - y)/n
    ndf = binHigh - binLow - f.GetNpar()
    return chi2, ndf
# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    parser.add_argument('--outputDirectory','-d', default='.',
                    help='Directory for output plots. If ".", save in the same directory of the input file.')
    parser.add_argument('--plot', '-p', action='store_true',
                    help='Boolean to produce very basic plot')
    args = parser.parse_args()


    # acces config
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    conf = inputCfg['fit']
    inputFile = conf['filename']
    if not inputFile:
        inputFile = 'Proj_' + os.path.basename(inputCfg['proj']['filename'])
    histName = conf['histname']
    if not histName:
        histName = inputCfg['proj']['histname'] + '_proj_' + str(inputCfg['proj']['massAxis'])
    fitLims = conf['fitLims']
    pdgId = conf['pdgId']


    # open root file and get Sparse
    file = TFile.Open(inputFile, 'UPDATE')
    hist = file.Get(histName)

    # one line fitter
    workspace, mass_frame = fit_invariant_mass_with_background(hist, fitLims[0], fitLims[1],pdgId, False)
    
    # retrieve parameters from workspace and create standard Canvas
    mass = workspace.var("mass")
    data = workspace.data("data")

    alpha = workspace.var("alpha")
    beta = workspace.var("beta")
    gamma = workspace.var("gamma")
    l = workspace.var("l")
    mTh = workspace.var("mTh")
    nBkg = workspace.var("nBkg")
    bkgPdf = workspace.pdf("bkgPdf")

    meanR = workspace.var("meanR")
    sigmaR = workspace.var("sigmaR")
    widthR = workspace.var("widthR")
    meanR.setConstant(1)
    sigmaR.setConstant(1)
    widthR.setConstant(1)
    nSigR = workspace.var("nSigR")
    sigRPdf = workspace.pdf("sigRPdf")

    totPdf = workspace.pdf("totPdf")

    if pdgId == 435:
        meanPart = workspace.var('meanPart')
        sigmaPart = workspace.var('sigmaPart')
        nSigPart = workspace.var('nSigPart')
        meanPart.setConstant(1)
        sigmaPart.setConstant(1)
        sigPartPdf = workspace.pdf("sigPartPdf")
    
    # model = RooStats.ModelConfig()
    # model.SetWorkspace(workspace)
    # model.SetPdf("totPdf")
    # poi = RooArgSet(nSigR)
    # nullParams = poi.snapshot()
    # nullParams.setRealValue("nSigR",0.)

    # plc = RooStats.ProfileLikelihoodCalculator()

    # plc.SetData(data)
    # plc.SetModel(model)
    # plc.SetParameters(poi)
    # plc.SetNullParameters(nullParams)

    # #We get a HypoTestResult out of the calculator, and we can query it.
    # hypo_test_result = plc.GetHypoTest()
    # print ("-------------------------------------------------")
    # 
    # print ("The p-value for the null is ", hypo_test_result.NullPValue())
    # print ("Corresponding to a signifcance of ", hypo_test_result.Significance())
    # print ("-------------------------------------------------")
    # #input()
    # del plc

    def fTot (x, par):
            mass.setVal(x[0])
            val = totPdf.getVal(RooArgSet(mass))
            return par[0]*val
    totTF = TF1("totTF", fTot,  fitLims[0],  fitLims[1], 1)
    totTF.SetParameter(0,  hist.Integral(hist.FindBin(fitLims[0]),  hist.FindBin(fitLims[1]))* (hist.GetBinCenter(2) - hist.GetBinCenter(1)))
    totTF.SetNpx(100000)
    
    if pdgId ==10433:
        f = nSigR.getVal() / (nSigR.getVal() + nBkg.getVal())
    if pdgId ==435:
        f = nSigR.getVal() / (nSigR.getVal() + nBkg.getVal() + nSigPart.getVal())
        f2 = nSigPart.getVal() / (nSigR.getVal() + nBkg.getVal() + nSigPart.getVal())
    def fSig(x, par):
        mass.setVal(x[0])
        val = sigRPdf.getVal(RooArgSet(mass))
        return par[0]*val
    sigTF = TF1("sigTF", fSig,  fitLims[0],  fitLims[1], 1)
    sigTF.SetNpx(100000)

    def fBkg(x, par):
        mass.setVal(x[0])
        val = bkgPdf.getVal(RooArgSet(mass))
        return par[0]*val
    bkgTF = TF1("bkgTF", fBkg,  fitLims[0],  fitLims[1], 1)
    bkgTF.SetNpx(100000)

    if pdgId == 10433:
        sigTF.SetParameter(0, f * hist.Integral(hist.FindBin(fitLims[0]), hist.FindBin(fitLims[1])) * (hist.GetBinCenter(2) - hist.GetBinCenter(1)))
        bkgTF.SetParameter(0, (1-f) * hist.Integral(hist.FindBin(fitLims[0]), hist.FindBin(fitLims[1])) * (hist.GetBinCenter(2) - hist.GetBinCenter(1)))        

    if pdgId == 435:
        sigTF.SetParameter(0, f * hist.Integral(hist.FindBin(fitLims[0]), hist.FindBin(fitLims[1])) * (hist.GetBinCenter(2) - hist.GetBinCenter(1)))
        bkgTF.SetParameter(0, (1- f - f2) * hist.Integral(hist.FindBin(fitLims[0]), hist.FindBin(fitLims[1])) * (hist.GetBinCenter(2) - hist.GetBinCenter(1)))
        def fPart (x, par):
            mass.setVal(x[0])
            val = sigPartPdf.getVal(RooArgSet(mass))
            return par[0]*val
        partTF = TF1("partTF", fPart,  fitLims[0],  fitLims[1], 1)
        partTF.SetParameter(0,  f2 * hist.Integral(hist.FindBin(fitLims[0]),  hist.FindBin(fitLims[1]))* (hist.GetBinCenter(2) - hist.GetBinCenter(1)))
        partTF.SetNpx(100000)
    
    chi2, ndf = getChi2Ndf(hist, totTF, fitLims[0],  fitLims[1])
    hwhm = 0.5346 * widthR.getVal() + (0.2166 * widthR.getVal()**2 + 2*np.log(2)*sigmaR.getVal()**2)**0.5 #from JQSRT 17, P233
    # hwhm = widthR.getVal()
    sig2hwhm = sigTF.Integral(meanR.getVal() - 2* hwhm, meanR.getVal() + 2* hwhm) / (hist.GetBinCenter(2) - hist.GetBinCenter(1))
    print(sig2hwhm)
    bkg2HWHM = bkgTF.Integral(meanR.getVal() - 2* hwhm, meanR.getVal() + 2* hwhm) / (hist.GetBinCenter(2) - hist.GetBinCenter(1))
    print(bkg2HWHM)
    signif2Hwhm = sig2hwhm / sqrt(bkg2HWHM + sig2hwhm)
    print ("-------------------------------------------------")
    print(f'Fit Chi2/ndf = {chi2}/{ndf}')
    print (f"The Raw Yield is {nSigR.getVal():.0f}")
    print (f"The statistical significance in 2HWHM is: {signif2Hwhm:.1f}")
    print(f"mu = {meanR.getVal():.5f} pm {meanR.getError():.5f}, sigma = {sigmaR.getVal():.5f} pm {sigmaR.getError():.5f}")
    if pdgId == 435:
        print (f"Partly reco Ds1: {nSigPart.getVal()}")
    print ("-------------------------------------------------")


    gStyle.SetOptStat(0)
    plotLims = conf['plotLims']
    hist.GetXaxis().SetRangeUser(plotLims[0], plotLims[1])
    c1 = TCanvas("c1", "Projection with Range", 3200, 2400)
    hist.Draw()
    totTF.Draw('same')
    sigTF.Draw('same')
    bkgTF.Draw('same')
    if pdgId == 435:
        partTF.Draw('same')
    if args.plot:
        c1.SaveAs("fitProva.png")  # Save the plot as a PNG file
    file.cd()
    c1.Write()
    file.Close()

    