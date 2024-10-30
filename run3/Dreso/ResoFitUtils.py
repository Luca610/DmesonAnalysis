'''
utility for fitting with Roofit D mesons particle spectra
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ROOT
from particle import Particle


def fit_bkg_only(histogram, xMin, xMax, pdgId):
    supported_ids = [10433, 435, 411] #[Ds1, Ds2*, D+]

    if pdgId not in supported_ids:
        raise KeyError(f"pdgId {pdgId} not supported!")
    
    mass = ROOT.RooRealVar("mass", "Invariant Mass", xMin, xMax)
    data = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), histogram)


    if pdgId == 10433:
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(413).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
    elif pdgId == 435:
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(411).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
    # backgorund --> Threshold function
    mTh = ROOT.RooRealVar("mTh", "Threshold mass", (mD + mV0))
    l = ROOT.RooRealVar("l", "exponent",0.5,  0.001, 1000)        
    alpha = ROOT.RooRealVar("alpha", "Linear Coefficient",0.2,  -10000, 0)
    beta = ROOT.RooRealVar("beta", "Quadratic Coefficient", 0.2, -10000, 10000)
    gamma = ROOT.RooRealVar("gamma", "Cubic Coefficient", 0.2, -10000, 10000)
    threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
    bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, mTh, alpha, beta, gamma, l))
    
    fit_result = bkgPdf.fitTo(data, ROOT.RooFit.Save())
    workspace = ROOT.RooWorkspace(f"w__{histogram.GetName()}", f"workspace_{histogram.GetName()}")
    getattr(workspace, 'import')(bkgPdf)
    getattr(workspace, 'import')(data)
    getattr(workspace, 'import')(fit_result, "fitResults")

    mass_frame = mass.frame()
    data.plotOn(mass_frame)
    bkgPdf.plotOn(mass_frame)

    # Draw the frame
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    canvas.cd(1)
    mass_frame.Draw()
    canvas.SaveAs("uglyplotBKG.pdf")
    return workspace, mass_frame

def fit_invariant_mass_with_background(histogram, x_min, x_max, pdgId, delta_m):
    """
    Fit invariant mass histogram using ROOFit.
    
    Parameters:
    - histogram: ROOT.TH1, histogram to be fitted.
    - x_min: float, minimum x-value (mass) to consider in the fit.
    - x_max: float, maximum x-value (mass) to consider in the fit.
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - mass_frame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    supported_ids = [10433, 435, 411] #[Ds1, Ds2*, D+]

    if pdgId not in supported_ids:
        raise KeyError(f"pdgId {pdgId} not supported!")
    
    mass = ROOT.RooRealVar("mass", "Invariant Mass", x_min, x_max)
    data = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), histogram)


    if pdgId == 10433:
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(413).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
        m2700 = 2.714
        sm2700 = 0.005
        w2700 = 0.122
        sw2700 = 0.01
        m2860 = 2.859
        sm2860 = 0.027
        w2860 = 0.160 #One component with width 0.053 also observed
        sw2860 = 0.080
        m3040 = 3.044
        sm3040 = 0.03
        w3040 = 0.240
        sw3040 = 0.06

        # PDFs
        # First peak Ds1(2536) --> Voigtian
        meanR = ROOT.RooRealVar("meanR", "Mean of Gaussian", mR, mR - 0.01, mR + 0.01) #, mR - 0.01, mR + 0.01)
        sigmaR = ROOT.RooRealVar("sigmaR", "Width of Gaussian", 0.01, 0., 0.05)
        widthR = ROOT.RooRealVar("widthR", "Width of BW", wR)
        sigRPdf = ROOT.RooVoigtian("sigRPdf", "Voigtian Pdf", mass, meanR, widthR, sigmaR)
        nSigR = ROOT.RooRealVar("nSigR","Number of Ds1(2536) events",0.5*histogram.Integral(), 0 ,2*histogram.Integral())
        # Second broad peak Ds1*(2700)
        mean2700 = ROOT.RooRealVar("mean2700", "Mean of Gaussian", m2700, m2700 - 2*sm2700 - 0.01, m2700 + 2*sm2700 + 0.01) #, mR - 0.01, mR + 0.01)
        sigma2700 = ROOT.RooRealVar("sigma2700", "Width of Gaussian", 0.01, 0.005, 0.02)
        width2700 = ROOT.RooRealVar("width2700", "Width of BW", w2700, w2700 - 2*sw2700, w2700 + 2*sw2700)
        sig2700Pdf = ROOT.RooVoigtian("sig2700Pdf", "Voigtian Pdf", mass, mean2700, width2700, sigma2700)
        nSig2700 = ROOT.RooRealVar("nSig2700","Number of Ds1*(2700) events",0.01*histogram.Integral(), 0 ,0.1*histogram.Integral()) 
        # Third broad peak Ds1*(2860)
        mean2860 = ROOT.RooRealVar("mean2860", "Mean of Gaussian", m2860, m2860 - 2*sm2860, m2860 + 2*sm2860) #, mR - 0.01, mR + 0.01)
        sigma2860 = ROOT.RooRealVar("sigma2860", "Width of Gaussian", 0.01, 0.005, 0.02)
        width2860 = ROOT.RooRealVar("width2860", "Width of BW", w2860, w2860 - 2*sw2860, w2860 + 2*sw2860)
        sig2860Pdf = ROOT.RooVoigtian("sig2860Pdf", "Voigtian Pdf", mass, mean2860, width2860, sigma2860)
        nSig2860 = ROOT.RooRealVar("nSig2860","Number of Ds1*(2860) events",0.01*histogram.Integral(), 0 ,0.1*histogram.Integral()) 
        # Fourth broad peak Dsj(3040)
        mean3040 = ROOT.RooRealVar("mean3040", "Mean of Gaussian", m3040, m3040 - 2*sm3040, m3040 + 2*sm3040) #, mR - 0.01, mR + 0.01)
        sigma3040 = ROOT.RooRealVar("sigma3040", "Width of Gaussian", 0.01, 0.005, 0.02)
        width3040 = ROOT.RooRealVar("width3040", "Width of BW", w3040, w3040 - 2*sw3040, w3040 + 2*sw3040)
        sig3040Pdf = ROOT.RooVoigtian("sig3040Pdf", "Voigtian Pdf", mass, mean3040, width3040, sigma3040)
        nSig3040 = ROOT.RooRealVar("nSig3040","Number of Dsj(3040) events",0.01*histogram.Integral(), 0 ,0.1*histogram.Integral()) 
        # backgorund --> Threshold function
        mTh = ROOT.RooRealVar("mTh", "Threshold mass", (mD + mV0))
        l = ROOT.RooRealVar("l", "exponent",0.5,  0.001, 1)        
        alpha = ROOT.RooRealVar("alpha", "Linear Coefficient",0.2,  -100, 0)
        beta = ROOT.RooRealVar("beta", "Quadratic Coefficient", 0.2, -100, 100)
        gamma = ROOT.RooRealVar("gamma", "Cubic Coefficient", 0.2, -1000, 1000)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, mTh, alpha, beta, gamma, l))
        nBkg = ROOT.RooRealVar("nBkg","Number of background events",0.5*histogram.Integral(), 0 ,2*histogram.Integral())

        # totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(sigRPdf, bkgPdf), ROOT.RooArgList(nSigR, nBkg))
        totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(sigRPdf, bkgPdf), ROOT.RooArgList(nSigR, nBkg))

    
    elif pdgId == 435:
        # voigtian signal and pol1 bkg
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(411).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
        mPart = 2.395

        # Main peak Ds2Star(2573) --> Voigtian
        meanR = ROOT.RooRealVar("meanR", "Mean of Gaussian", mR,  mR - 0.01, mR + 0.01) #, mR - 0.01, mR + 0.01)
        sigmaR = ROOT.RooRealVar("sigmaR", "Width of Gaussian", 0.01, 0., 0.02)
        widthR = ROOT.RooRealVar("widthR", "Width of BW", wR)
        sigRPdf = ROOT.RooVoigtian("sigRPdf", "Voigtian Pdf", mass, meanR, widthR, sigmaR)
        nSigR = ROOT.RooRealVar("nSigR","Number of Ds1(2536) events",0.01*histogram.Integral(), 0 ,0.1*histogram.Integral())
        # Partial Ds1*(2536)
        meanPart = ROOT.RooRealVar("meanPart", "Mean of Gaussian", mPart, mPart - 0.005, mPart + 0.005) #, mR - 0.01, mR + 0.01)
        sigmaPart = ROOT.RooRealVar("sigmaPart", "Width of Gaussian", 0.01, 0.0, 0.02)
        sigPartPdf = ROOT.RooGaussian("sigPartPdf", "Gaussian Pdf", mass, meanPart,  sigmaPart)
        nSigPart = ROOT.RooRealVar("nSigPart","Number of Ds1(Part) events",0.01*histogram.Integral(), 0 ,0.1*histogram.Integral()) 
        # backgorund --> Threshold function
        mTh = ROOT.RooRealVar("mTh", "Threshold mass", (mD + mV0))
        l = ROOT.RooRealVar("l", "exponent",0.5,  0.001, 1)        
        alpha = ROOT.RooRealVar("alpha", "Linear Coefficient",0.2,  -10, 0)
        beta = ROOT.RooRealVar("beta", "Quadratic Coefficient", 0.2, -100, 100)
        gamma = ROOT.RooRealVar("gamma", "Cubic Coefficient", 0.2, -100, 100)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
        bkgPdf = ROOT.RooGenericPdf("bkgPdf", "Threshold Pdf", threshold_formula, ROOT.RooArgList(mass, mTh, alpha, beta, gamma, l))
        nBkg = ROOT.RooRealVar("nBkg","Number of background events",0.9*histogram.Integral(), 0 ,histogram.Integral())

        totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(sigRPdf, sigPartPdf, bkgPdf), ROOT.RooArgList(nSigR, nSigPart, nBkg))


    # elif pdgId == 411:
    #     # Gaussian signal and pol3 bkg
    #     m = Particle.from_pdgid(411).mass*1e-3
    #     mean = ROOT.RooRealVar("mean", "Mean of Gaussian", m, m - 0.02, m + 0.02)
    #     sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0., 0.02)
    #     signal_pdf = ROOT.RooGaussian("signal_pdf", "Gaussian Model", mass, mean, sigma)
    #     a0 = ROOT.RooRealVar("a0", "Constant Coefficient", 100, 0, 1e7)
    #     a1 = ROOT.RooRealVar("a1", "Linear Coefficient", 0, -1e3, 1e3)
    #     a2 = ROOT.RooRealVar("a2", "Quadratic Coefficient", 0, -1e3, 1e3)
    #     a3 = ROOT.RooRealVar("a3", "Cubic Coefficient", 0, -1e3, 1e3)
    #     bkg_pdf = ROOT.RooPolynomial("bkg_pdf", "Polynomial Background", mass, ROOT.RooArgList(a0, a1, a2, a3))
    #     width = ROOT.RooRealVar("width", "Width of BW", 0)

    # # Combine signal and background into a composite model
    # Nsig = ROOT.RooRealVar("Nsig","Number of signal events",0.01*histogram.Integral(), 0 ,0.3*histogram.Integral())
    # Nbkg = ROOT.RooRealVar("Nbkg","Number of backgrouund events",0.9*histogram.Integral(), 0 ,histogram.Integral())
    # totPdf = ROOT.RooAddPdf("totPdf", "Total PDF with Nsig as parameter", ROOT.RooArgList(signal_pdf, bkg_pdf), ROOT.RooArgList(Nsig, Nbkg))

    # Fit the composite model to the data
    fit_result = totPdf.fitTo(data, ROOT.RooFit.Save())

    # Create a frame to draw the fit result and data
    mass_frame = mass.frame()
    data.plotOn(mass_frame)
    totPdf.plotOn(mass_frame, ROOT.RooFit.Components("sigRPdf"), ROOT.RooFit.LineColor(ROOT.kRed))
    if pdgId == 10433:
        totPdf.plotOn(mass_frame, ROOT.RooFit.Components("sig2700Pdf"), ROOT.RooFit.LineColor(ROOT.kGreen))
    if pdgId == 435:
        totPdf.plotOn(mass_frame, ROOT.RooFit.Components("sigPartPdf"), ROOT.RooFit.LineColor(ROOT.kGreen))
    totPdf.plotOn(mass_frame, ROOT.RooFit.Components("bkgPdf"), ROOT.RooFit.LineStyle(ROOT.kDashed))
    totPdf.plotOn(mass_frame)
    

    hresid = mass_frame.residHist()
 
    # Construct a histogram with the pulls of the data w.r.t the curve
    hpull = mass_frame.pullHist()

    # Save to workspace
    workspace = ROOT.RooWorkspace(f"w__{histogram.GetName()}", f"workspace_{histogram.GetName()}")
    getattr(workspace, 'import')(totPdf)
    getattr(workspace, 'import')(data)
    getattr(workspace, 'import')(fit_result, "fitResults")  # Give a name to your fit results for later access

    # Draw the frame
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    canvas.Divide(3,1)
    canvas.cd(1)
    mass_frame.Draw()
    canvas.cd(2)
    hresid.Draw()
    canvas.cd(3)
    hpull.Draw()
    canvas.SaveAs("uglyplot.pdf")
    return workspace, mass_frame

def simultaneous_fit_signal_background(signalHist, bkgHist, xMin, xMax, pdgId):
    """
    Fit invariant mass histogram using ROOFit, Background estimated with simultaneous fit 
    from another histogram (sidebands, wrong sign, event mixing, rotation).
    
    Arguments:
    - signalHist: ROOT.TH1, invariant mass histogram to be fitted.
    - bkgHist: ROOT.TH1, BKG invariant mass distribution.
    - xMin: float, lower fit limit
    - xMax: float, upper fit limit
    - pdgId: int, PDG code of resonance to be fitted. Used to determine fitting functions and initialize parameters.
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - massFrame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """

    supportedIds = [10433, 435, 411] #[Ds1, Ds2*, D+]

    if pdgId not in supportedIds:
        raise KeyError(f"pdgId {pdgId} not supported!")
    
    mass = ROOT.RooRealVar("mass", "Invariant Mass", xMin, xMax)
    signalData = ROOT.RooDataHist("signalData", "Dataset with invariant mass", ROOT.RooArgList(mass), signalHist)
    bkgData = ROOT.RooDataHist("bkgData", "Dataset with invariant mass", ROOT.RooArgList(mass), bkgHist)
    if pdgId == 435:
        # physics constants
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(411).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
        mSig2 = 2.39

        # PDFs
        # signal Ds2*(2573) --> Voigtian
        mean = ROOT.RooRealVar("mean", "Mean of Gaussian", mR) #, mR - 0.01, mR + 0.01)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0.005, 0.02)
        width = ROOT.RooRealVar("width", "Width of BW", wR)
        signalPDF = ROOT.RooVoigtian("signalPDF", "Voigtian PDF", mass, mean, width, sigma)
        # Sig2 Incomplete Ds1 --> Gaussian
        meanSig2 = ROOT.RooRealVar("meanSig2", "Mean of Gaussian", mSig2, mSig2 - 0.01, mSig2 + 0.01)
        sigmaSig2 = ROOT.RooRealVar("sigmaSig2", "Width of Gaussian", 0.01, 1e-12, 0.01)
        sig2PDF = ROOT.RooGaussian("sig2PDF", "Gaussian Model", mass, meanSig2, sigmaSig2)
        # backgorund --> Threshold function
        mTh = ROOT.RooRealVar("mTh", "Threshold mass", (mD + mV0))
        l = ROOT.RooRealVar("l", "exponent",0.5,  0.001, 1)        
        alpha = ROOT.RooRealVar("alpha", "Linear Coefficient",0.2,  -10, 0)
        beta = ROOT.RooRealVar("beta", "Quadratic Coefficient", 0.2, -100, 100)
        gamma = ROOT.RooRealVar("gamma", "Cubic Coefficient", 0.2, -100, 100)
        threshold_formula = "(mass - mTh)^l * exp(alpha * (mass - mTh) + beta * (mass - mTh) * (mass - mTh)+ gamma * (mass - mTh)* (mass - mTh)* (mass - mTh))"
        # thresholdPDF = ROOT.RooGenericPdf("thresholdPDF", "resonance bkg PDF", "sqrt(abs(mass-mTh))*exp(alpha*(mass-mTh))", ROOT.RooArgSet(mass, mTh, alpha))

        thresholdPDF = ROOT.RooGenericPdf("thresholdPDF", "Threshold PDF", threshold_formula, ROOT.RooArgList(mass, mTh, alpha, beta, gamma, l))
        # bkg_pdf = ROOT.RooPolynomial("bkg_pdf", "Polynomial Background", mass, ROOT.RooArgList(a0, a1))
        # totalPDF (signal region only)
        nSig = ROOT.RooRealVar("nSig","Number of signal events",0.01*signalHist.Integral(), 0 ,0.1*signalHist.Integral())
        
        nSig2 = ROOT.RooRealVar("nSig2","Number of signal 2 events",0.01*signalHist.Integral(), 0 ,0.1*signalHist.Integral())
        nBkg = ROOT.RooRealVar("nBkg","Number of background events",0.9*signalHist.Integral(), 0 ,signalHist.Integral())
        # totPDF = ROOT.RooAddPdf("totPDF", "Total PDF", ROOT.RooArgList(signalPDF, sig2PDF, thresholdPDF), ROOT.RooArgList(nSig, nSig2, nBkg))        
        totPDF = ROOT.RooAddPdf("totPDF", "Total PDF", ROOT.RooArgList(signalPDF, thresholdPDF), ROOT.RooArgList(nSig,nBkg))        

        sigCat = ROOT.RooCategory("sigCat", "signal Categories")
        sigCat.defineType("signalRegion")
        sigCat.defineType("sidebandRegion")
        combinedData = ROOT.RooDataHist("combinedData", "combined data", ROOT.RooArgList(mass), ROOT.RooFit.Index(sigCat),
                    ROOT.RooFit.Import("signalRegion", signalData), ROOT.RooFit.Import("sidebandRegion", bkgData))
        # Create a simultaneous PDF using the category
        simPDF = ROOT.RooSimultaneous("simPDF", "simultaneous pdf", sigCat)
        simPDF.addPdf(totPDF, "signalRegion")
        simPDF.addPdf(thresholdPDF, "sidebandRegion")
        # Actual Fit
        fitResult = simPDF.fitTo(combinedData)
        # Save to workspace
        workspace = ROOT.RooWorkspace(f"w", f"workspace")
        getattr(workspace, 'import')(simPDF)
        getattr(workspace, 'import')(totPDF)
        getattr(workspace, 'import')(thresholdPDF)
        getattr(workspace, 'import')(bkgData)
        getattr(workspace, 'import')(signalData)
        # getattr(workspace, 'import')(fitResult, "fitResults")  # Give a name to your fit results for later access
        fOut = ROOT.TFile("prova.root","RECREATE")
        fOut.cd()
        workspace.Write()
        fOut.Close()
        #Do the plotting
        #together 
        massframe3 = mass.frame()
        massframe3.SetTitle("Combined histogram")
        combinedData.plotOn(massframe3)
        simPDF.plotOn(massframe3,ROOT.RooFit.ProjWData(sigCat,combinedData, True))
        #Separately for the two categories
        massframe1 = mass.frame()
        massframe1.SetTitle("Signal Region")
        combinedData.plotOn(massframe1, ROOT.RooFit.Cut("sigCat==sigCat::signalRegion"))
        simPDF.plotOn(massframe1, ROOT.RooFit.Slice(sigCat,"signalRegion"), ROOT.RooFit.ProjWData(sigCat,combinedData,True), ROOT.RooFit.Components("totPDF"), ROOT.RooFit.LineColor(ROOT.kRed))
        simPDF.plotOn(massframe1, ROOT.RooFit.Slice(sigCat,"signalRegion"), ROOT.RooFit.ProjWData(sigCat,combinedData,True), ROOT.RooFit.Components("sigPDF"), ROOT.RooFit.LineColor(ROOT.kGreen))
        # simPDF.plotOn(massframe1, ROOT.RooFit.Slice(sigCat,"signalRegion"), ROOT.RooFit.ProjWData(sigCat,combinedData,True), ROOT.RooFit.Components("sig2PDF"), ROOT.RooFit.LineColor(ROOT.kGreen))
        simPDF.plotOn(massframe1, ROOT.RooFit.Slice(sigCat,"signalRegion"), ROOT.RooFit.ProjWData(sigCat,combinedData,True), ROOT.RooFit.Components("thresholdPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed))
        massframe2 = mass.frame()
        massframe2.SetTitle("Side Band region")
        combinedData.plotOn(massframe2, ROOT.RooFit.Cut("sigCat==sigCat::sidebandRegion"))
        simPDF.plotOn(massframe2, ROOT.RooFit.Slice(sigCat,"sidebandRegion"), ROOT.RooFit.ProjWData(sigCat,combinedData,True), ROOT.RooFit.Components("thresholdPDF"))
        canvas = ROOT.TCanvas("cQa","QA",1650,900)
        canvas.Divide(3,1)
        canvas.cd(1)
        massframe1.Draw()
        canvas.cd(2)
        massframe2.Draw()
        canvas.cd(3)
        massframe3.Draw()
        canvas.SaveAs("symfit_prova.png")
        nSig.Print()
        nBkg.Print()
        

        


def template_fit_signal_background(signalHist, bkgHist, xMin, xMax, pdgId):
    """
    Fit invariant mass histogram using ROOFit, Background estimated with simultaneous fit 
    from another histogram (sidebands, wrong sign, event mixing, rotation).

    Arguments:
    - signalHist: ROOT.TH1, invariant mass histogram to be fitted.
    - bkgHist: ROOT.TH1, BKG invariant mass distribution.
    - xMin: float, lower fit limit
    - xMax: float, upper fit limit
    - pdgId: int, PDG code of resonance to be fitted. Used to determine fitting functions and initialize parameters.

    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - massFrame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    supportedIds = [10433, 435, 411] #[Ds1, Ds2*, D+]

    if pdgId not in supportedIds:
        raise KeyError(f"pdgId {pdgId} not supported!")

    mass = ROOT.RooRealVar("mass", "Invariant Mass", xMin, xMax)
    signalData = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), signalHist)
    bkgData = ROOT.RooDataHist("data", "Dataset with invariant mass", ROOT.RooArgList(mass), bkgHist)

    if pdgId == 435:
        # physics constants
        mR = Particle.from_pdgid(pdgId).mass*1e-3
        wR = Particle.from_pdgid(pdgId).width/2*1e-3
        mD = Particle.from_pdgid(411).mass*1e-3
        mV0 = Particle.from_pdgid(310).mass*1e-3
        mSig2 = 2.39

        # PDFs
        # signal Ds2*(2573) --> Voigtian
        mean = ROOT.RooRealVar("mean", "Mean of Gaussian", mR) #, mR - 0.01, mR + 0.01)
        sigma = ROOT.RooRealVar("sigma", "Width of Gaussian", 0.01, 0.005, 0.02)
        width = ROOT.RooRealVar("width", "Width of BW", wR)
        signalPDF = ROOT.RooVoigtian("signalPDF", "Voigtian PDF", mass, mean, width, sigma)
        bkgPDF = ROOT.RooHistPdf("bkgPDF", "Template PDF", ROOT.RooArgSet(mass), bkgData) #, 2 per interpolare
        nSig = ROOT.RooRealVar("nSig","Number of signal events",0.01*signalHist.Integral(), 0 ,0.1*signalHist.Integral())
        nBkg = ROOT.RooRealVar("nBkg","Number of background events",0.9*signalHist.Integral(), 0 ,signalHist.Integral())
        totPDF = ROOT.RooAddPdf("totPDF", "Total PDF", ROOT.RooArgList(signalPDF, bkgPDF), ROOT.RooArgList(nSig,nBkg))
        fit_result = totPDF.fitTo (signalData)
        # Create a frame to draw the fit result and data
        mass_frame = mass.frame()
        signalData.plotOn(mass_frame)
        totPDF.plotOn(mass_frame, ROOT.RooFit.Components("bkgPDF"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen)) 
        totPDF.plotOn(mass_frame, ROOT.RooFit.Components("signalPDF"), ROOT.RooFit.LineColor(ROOT.kRed))
        totPDF.plotOn(mass_frame)
        canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
        mass_frame.Draw()
        canvas.SaveAs("template_test.png")      


def roofit_plot_with_matplotlib(axis, histogram, workspace, x_min, x_max, fit_min, fit_max, plot):
    """
    Plot a TH1 and ROOFit fitted PDFs using matplotlib .
    
    Parameters:
    - axis: plt.axis, axis where to draw the plot.
    - histogram: ROOT.TH1, Invariant mass histogramselectedMasses = invMassReso[(ptReso>ptMin) & (ptReso<ptMax)]
    - x_min: float, minimum of mass plot range.
    - x_max: float, maximum of mass plot range.
    - fit_min: float, minimum of mass fit range (for normalization porpuses).
    - fit_max: float, maximum of mass fit range (for normalization porpuses).
    
    Returns:
    - workspace: ROOT.RooWorkspace, workspace containing the model and fit result.
    - mass_frame: ROOT.RooPlot, frame containing the fit and histogram for plotting.
    """
    
    # Access mass variable and data
    mass = workspace.var("mass")
    data = workspace.data("data")
    # frac = workspace.var("signal_fraction")
    nsig = workspace.var("Nsig")
    nbkg = workspace.var("Nbkg")
    fit_result = workspace.obj("fitResults")
   
    # Create the underlying histogram
    Nmin = histogram.FindBin(x_min)
    Nmax = histogram.FindBin(x_max)
    N = Nmax - Nmin
    N_fitmin = histogram.FindBin(fit_min)
    N_fitmax = histogram.FindBin(fit_max)
    N_fit = N_fitmax - N_fitmin
    N_plot = 250
    bin_center = np.zeros(N)
    bin_content = np.zeros(N)
    bin_error = np.zeros(N)
    # Accessing each data point
    for i in range(N):
        bin_center[i] = histogram.GetBinCenter(Nmin + i)
        bin_content[i] = histogram.GetBinContent(Nmin + i)
        bin_error[i] = histogram.GetBinError(Nmin + i)

    totPDF = workspace.pdf("totalPDF")
    signalPDF = workspace.pdf("signal_pdf")
    bkgPDF = workspace.pdf("bkg_pdf")
    y_totPDF = np.zeros(N_fit)
    y_signalPDF = np.zeros(N_fit)
    y_bkgPDF = np.zeros(N_fit)
    x_fit = np.linspace(fit_min, fit_max, N_fit)
    x_plot = np.linspace(fit_min, fit_max, N_plot)
    y_tot_plot = np.zeros(N_plot)
    y_signal_plot = np.zeros(N_plot)
    y_bkg_plot = np.zeros(N_plot)
    for i, x in enumerate(x_fit):
        mass.setVal(x)
        y_totPDF[i] = totPDF.getVal(ROOT.RooArgSet(mass)) 
        y_signalPDF[i] = signalPDF.getVal(ROOT.RooArgSet(mass)) 
        y_bkgPDF[i] = bkgPDF.getVal(ROOT.RooArgSet(mass)) 
    for i, x in enumerate(x_plot):
        mass.setVal(x)
        y_tot_plot[i] = totPDF.getVal(ROOT.RooArgSet(mass)) 
        y_signal_plot[i] = signalPDF.getVal(ROOT.RooArgSet(mass)) 
        y_bkg_plot[i] = bkgPDF.getVal(ROOT.RooArgSet(mass)) 

    # Retrieving useful information for legend 
    chi2 = totPDF.createChi2(data)
    mean = workspace.var("mean")
    sigma = workspace.var("sigma")
    width = workspace.var("width")
    
    # Compute S, B and Significance
    hwhm = 0.5346 * width.getVal() + (0.2166 * width.getVal()**2 + 2*np.log(2)*sigma.getVal()**2)**0.5 #from JQSRT 17, P233
    mass.setRange("intRange", mean.getVal() - 3*hwhm , mean.getVal() + 3*hwhm)
    bkg_3HWHM = bkgPDF.createIntegral(ROOT.RooArgSet(mass), ROOT.RooFit.NormSet(ROOT.RooArgSet(mass)), ROOT.RooFit.Range("intRange"))
    bkg = nbkg.getVal()*bkg_3HWHM.getVal()
    bkg_err = nbkg.getError()*bkg_3HWHM.getVal()
    sig = nsig.getVal()
    sig_err = nsig.getError()
    sig_plus_bkg = nsig.getVal()+nbkg.getVal()*bkg_3HWHM.getVal() 
    if sig <= 0 or bkg <= 0:
        return [0, 0, 0, 0]
    if chi2.getVal() < 2* (N_fit-fit_result.floatParsFinal().getSize()):
        significance = nsig.getVal()/np.sqrt(sig_plus_bkg)
        significance_err = significance*np.sqrt((sig_err**2 + bkg_err**2) / (4. * sig_plus_bkg**2) + (bkg/sig_plus_bkg) * sig_err**2 / sig**2)
    else:
        significance = 0
        significance_err = 0
        # sig = 0
        # sig_err = 0
    if plot:
        #Create legend text
        textstr1 = '\n'.join(( 
            '\u03BC = %.1f \u00B1 %.1f MeV/c$^2$' % (mean.getVal([0])*1000, mean.getError()*1000),
            '\u03C3 =  %.1f \u00B1 %.1f MeV/c$^2$' % (sigma.getVal([0])*1000, sigma.getError()*1000),
            '\u0393 =  %.1f \u00B1 %.1f MeV/c$^2$' % (2*width.getVal([0])*1000, width.getError()*1000),
            'HWHM = %.1f MeV/c$^2$'  % (hwhm * 1000),
            'S = %d \u00B1 %d' % (int(nsig.getVal()), int(nsig.getError())),
            'B(3 HWHM) = %d \u00B1 %d' % (int(nbkg.getVal()*bkg_3HWHM.getVal()), int((nbkg.getError()*bkg_3HWHM.getVal()))),
            'B = %d \u00B1 %d' % (int(nbkg.getVal()), int(nbkg.getError())),
                ))
        textstr2 = '\n'.join(( 
            '',
            '\u03C7$^2$/ndf = %.1f/%d' % (chi2.getVal(), N_fit-fit_result.floatParsFinal().getSize()),
            'Significance(3 HWHM) = %.1f \u00B1 %.1f' % (significance, significance_err),
                ))
        axis.text(0.015, 0.95, textstr2, transform=axis.transAxes, fontsize=15, verticalalignment='top')

        # Plotting data and fit PDFs
        axis.errorbar(bin_center, bin_content, yerr=bin_error, fmt='o', label='Data', markersize=5, color='black')
        
        norm = histogram.Integral(N_fitmin,N_fitmax)/np.sum(y_totPDF)
        f = nsig.getVal()/((nbkg.getVal() + nsig.getVal()))  

        axis.plot(x_plot, (1-f)*norm*y_bkg_plot, label='Bkg', linewidth=2, linestyle='--', color='red')
        axis.plot(x_plot, f*norm*y_signal_plot, label='Sig', linewidth=1, color='blue')
        axis.plot(x_plot, norm*y_tot_plot, label='Sig + Bkg', linewidth=2, color='blue')
        axis.fill_between(x_plot, f*norm*y_signal_plot, color='skyblue', alpha=0.4)
        axis.legend(fontsize=15)

        # Plot Options
        axis.set_title(f"{histogram.GetName()}", fontsize=18)
        # axis.set_xlabel(rf"{label_mass} (GeV/$c^2$)", fontsize=18)
        axis.set_xlabel(f"Invariant mass (GeV/c²)", fontsize=18)
        axis.set_ylabel(f'Counts/{((bin_center[2]-bin_center[1])*1000):.1f} MeV/c²', fontsize=18)
        axis.set_ylim(0., max(bin_content) * 1.5)
        props = dict(boxstyle='square', facecolor='white', alpha=0.8, edgecolor='black', zorder=10)
        axis.text(0.75, 0.3, textstr1, transform=axis.transAxes, fontsize=15, verticalalignment='top', color='black', bbox=props, backgroundcolor='white')

    results = [int(sig), int(sig_err), significance , significance_err ]
    return results