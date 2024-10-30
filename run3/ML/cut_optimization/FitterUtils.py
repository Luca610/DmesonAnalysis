'''
utility for fitting with Roofit D mesons particle spectra
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ROOT
from particle import Particle

def sigma_par(pT):
    return (0.00796 + 0.00176*pT) * 1000

def fitInvMass(histogram, xMin, xMax, xFitMin, xFitMax, pdgId, pt = 2):
    """
    Fit invariant mass histogram using ROOFit.
    
    Parameters:
    - histogram: ROOT.TH1, histogram to be fitted.
    - xMin: float, minimum x-value (mass)
    - xMax: float, maximum x-value (mass)
    - xFitMin: float, minimum x-value of fit range
    - xFitMax: float, maximum x-value of fit range
    - pdgId: id to distinguish between particles
    - pt: pt of bin for sigma parametrization
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - massFrame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    supported_ids = [411, 413] #D+, D*

    if pdgId not in supported_ids:
        raise KeyError(f"pdgId {pdgId} not supported!")
    
    mass = ROOT.RooRealVar("mass", "Invariant Mass", xMin, xMax)
    data = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), histogram)

    if pdgId == 413: #DStar
        # Crystal Bal signal and Threshold bkg
        m = (Particle.from_pdgid(pdgId).mass - Particle.from_pdgid(421).mass)*1e-3
        mPi = Particle.from_pdgid(211).mass*1e-3
        s = sigma_par(pt)/8
        mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m, m - 0.002, m + 0.002)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.001, 0., 0.05)
        alphaL = ROOT.RooRealVar("alphaL", "Left Boundry of Gaussian Core", 0.001, 0, 100)
        alphaR = ROOT.RooRealVar("alphaR", "Right Boundry of Gaussian Core", 0.001, 0, 100)
        nL = ROOT.RooRealVar("nL", "Slope of Left tail", 1, 0, 10)
        nR = ROOT.RooRealVar("nR", "Slope of Right tail", 1, 0, 10)
        signalPdf = ROOT.RooCrystalBall("signalPdf", "Gaussian Model", mass, mean, sigma, alphaL, nL, alphaR, nR)
        
        mTh = ROOT.RooRealVar("mTh", "Threshold mass", mPi)
        l = ROOT.RooRealVar("l", "exponent",0.5,  0.001, 1)        
        alpha = ROOT.RooRealVar("alpha", "Linear Coefficient",0.2,  -100, 0)
        beta = ROOT.RooRealVar("beta", "Quadratic Coefficient", 0.2, -100, 100)
        gamma = ROOT.RooRealVar("gamma", "Cubic Coefficient", 0.2, -1000, 1000)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, mTh, alpha, beta, gamma, l))
    
    elif pdgId == 411: #DPlus
        # Gaussian signal and pol3 bkg
        m = Particle.from_pdgid(pdgId).mass*1e-3
        s = sigma_par(pt)
        mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m, m - 0.05, m + 0.05)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.05)
        signalPdf = ROOT.RooGaussian("signalPdf", "Gaussian Model", mass, mean, sigma)
        a0 = ROOT.RooRealVar("a0", "Constant Coefficient", 100, 0, 1e7)
        a1 = ROOT.RooRealVar("a1", "Linear Coefficient", 0, -1e3, 1e3)
        a2 = ROOT.RooRealVar("a2", "Quadratic Coefficient", 0, -1e3, 1e3)
        a3 = ROOT.RooRealVar("a3", "Cubic Coefficient", 0, -1e3, 1e3)
        bkgPdf = ROOT.RooPolynomial("bkgPdf", "Polynomial Background", mass, ROOT.RooArgList(a0, a1, a2, a3))

    # Combine signal and background into a composite model
    nSig = ROOT.RooRealVar("nSig","Number of signal events",0.01*histogram.Integral(), 0 , histogram.Integral())
    nBkg = ROOT.RooRealVar("nBkg","Number of backgrouund events",0.9*histogram.Integral(), 0 , histogram.Integral())
    totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with nSig as parameter", ROOT.RooArgList(signalPdf, bkgPdf), ROOT.RooArgList(nSig, nBkg))
    print("fitting")
    # Fit the composite model to the data
    mass.setRange("fitRange", xFitMin, xFitMax)
    fitResult = totPdf.fitTo(data, ROOT.RooFit.Range("fitRange"))
    
    print("plotting")
    # Create a frame to draw the fit result and data
    massFrame = mass.frame()
    data.plotOn(massFrame)
    totPdf.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kBlue))
    # totPdf.plotOn(massFrame, ROOT.RooFit.Components("signalPdf"), ROOT.RooFit.LineColor(ROOT.kGreen))
    # totPdf.plotOn(massFrame, ROOT.RooFit.Components("bkgPdf"), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed))
    hresid = massFrame.residHist() 
    hpull = massFrame.pullHist()
    print("saving")
    # Save to workspace
    workspace = ROOT.RooWorkspace(f"w__{histogram.GetName()}", f"workspace_{histogram.GetName()}")
    getattr(workspace, 'import')(totPdf)
    getattr(workspace, 'import')(data)

    # Draw the frame
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    canvas.Divide(3,1)
    canvas.cd(1)
    massFrame.Draw()
    canvas.cd(2)
    hresid.Draw()
    canvas.cd(3)
    hpull.Draw()
    canvas.SaveAs("uglyplot.pdf")
    # print(massFrame.chiSquare())
    # input()

    print("done")
    return workspace, massFrame
