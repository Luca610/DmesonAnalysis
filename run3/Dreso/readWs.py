import ROOT

fInput = ROOT.TFile("prova.root")
fInput.cd()

workspace = fInput.Get("w")

workspace.Print()

mass = workspace.var("mass")
nBkg = workspace.var("nBkg")
nSig = workspace.var("nSig")

alpha = workspace.var("alpha")
beta = workspace.var("beta")
gamma = workspace.var("gamma")
l = workspace.var("l")
mTh = workspace.var("mTh")
mean = workspace.var("mean")
sigCat = workspace.var("sigCat")
sigma = workspace.var("sigma")
width = workspace.var("width")

totPDF = workspace.pdf("totPDF")
bkgPDF = workspace.pdf("thresholdPDF")

data = workspace.data("signalData")
dataBkg = workspace.data("bkgData")

# alpha.setConstant(1)
# beta.setConstant(1)
# gamma.setConstant(1)
# l.setConstant(1)
# mTh.setConstant(1)
# mean.setConstant(1)
# sigma.setConstant(0)
# width.setConstant(1)
# nBkg.setConstant(0)



# model = ROOT.RooStats.ModelConfig()
# model.SetWorkspace(workspace)
# model.SetPdf("totPDF")
# poi = ROOT.RooArgSet(nSig)
# nullParams = poi.snapshot()
# nullParams.setRealValue("nSig",0.)

# plc = ROOT.RooStats.ProfileLikelihoodCalculator()

# plc.SetData(data)
# plc.SetModel(model)
# plc.SetParameters(poi)
# plc.SetNullParameters(nullParams)

# #We get a HypoTestResult out of the calculator, and we can query it.
# hypo_test_result = plc.GetHypoTest()
# print ("-------------------------------------------------")
# print ("The p-value for the null is ", hypo_test_result.NullPValue())
# print ("Corresponding to a signifcance of ", hypo_test_result.Significance())
# print ("-------------------------------------------------")
# #input()
# del plc

# massplotSig = mass.frame()
# data.plotOn(massplotSig)
# totPDF.plotOn(massplotSig)
# chi2 = massplotSig.chiSquare()

# massplotBkg = mass.frame()
# dataBkg.plotOn(massplotBkg)
# bkgPDF.plotOn(massplotBkg)
# chi2Bkg = massplotBkg.chiSquare()

# canvas = ROOT.TCanvas()
# canvas.Divide(3,1)
# canvas.cd(1)
# massplotSig.Draw()
# canvas.cd(2)
# massplotBkg.Draw()

# nll = totPDF.createNLL(data)

# # Create the profile likelihood object
# profile_ll = ROOT.RooProfileLL("profile_ll", "Profile Likelihood of Nsig", nll, ROOT.RooArgSet(nSig))

# # Plot the profile likelihood
# frame = nSig.frame(ROOT.RooFit.Title("Profile Likelihood of Nsig"))
# profile_ll.plotOn(frame, ROOT.RooFit.ShiftToZero())  # ShiftToZero plots the profile likelihood relative to its minimum
# canvas.cd(3)
# frame.Draw()
# line_95 = ROOT.TLine(nSig.getMin(), 1.92, nSig.getMax(), 1.92)
# line_95.SetLineColor(ROOT.kRed)
# line_95.SetLineStyle(2)
# line_95.Draw("same")
# canvas.SaveAs("exercise_3.png")


histo = data.createHistogram("histo", mass)
histo.SetTitle("Data and PDF Comparison; mass; Events")
histo.Scale(1/histo.Integral())
histo.SetMarkerStyle(20)

f = nSig.getVal() / (nBkg.getVal() + nSig.getVal())
fSig = ROOT.RooRealVar("fSig", "fSig",f)
# Convert the PDF to a TF1 function
nSig.Print()
nBkg.Print()
# totPDF.fixCoefNormalization({nSig, nBkg})
pdf_function = totPDF.asTF(ROOT.RooArgList(mass),ROOT.RooArgList(alpha ,beta ,gamma ,l ,mTh ,mean ,sigma ,width, nSig, nBkg),ROOT.RooArgList(mass))

pdf_function.SetLineColor(ROOT.kRed)
pdf_function.SetLineWidth(2)
x_min = mass.getMin()
x_max = mass.getMax()



def normalized_formula(x, par):
    return par[0]* rooPdf.Eval(x[0])
normTF = ROOT.TF1(f"{rooTF.GetName()}_normalized", normalized_formula,xMin,xMax,rooTF.GetNpar()+1)
normTF.SetParameter(0, frac * binWidth)
return normTF

binWidth = histo.GetBinCenter(2) - histo.GetBinCenter(1)

# def normalized_function(x, par):
#     # par[0] is the normalization factor
#     normalization_factor = par[0]
#     # Evaluate the original function at x
#     original_value = pdf_function.Eval(x[0])
#     return normalization_factor * original_value
# # 
# # Convert the normalized function to a ROOT TF1
# npar = 1  # We have 1 parameter for normalization
# new_func = ROOT.TF1("newFunc", normalized_function, x_min, x_max, npar)
# new_func.SetParameter(0,(binWidth))

pars = ROOT.RooArgList(alpha ,beta ,gamma ,l ,mTh ,mean ,sigma ,width, nSig, nBkg)
new_func = normalized_TF1(totPDF, mass, pars, x_min, x_max, binWidth, 1)

sigPdf = workspace.pdf("signalPDF").asTF(ROOT.RooArgList(mass),ROOT.RooArgList(mean ,width, sigma),ROOT.RooArgList(mass))
def fSig (x, par):
    return par[0]*sigPdf.Eval(x[0])
sigTF = ROOT.TF1("sigTF", fSig, x_min, x_max, 1)
sigTF.SetParameter(0,  f* (histo.GetBinCenter(2) - histo.GetBinCenter(1)) )
histo.GetYaxis().SetRangeUser(0, 0.02)

bkgPdf = workspace.pdf("thresholdPDF").asTF(ROOT.RooArgList(mass),ROOT.RooArgList(alpha ,beta ,gamma ,l ,mTh ),ROOT.RooArgList(mass))
def fBkg (x, par):
    return par[0]*bkgPdf.Eval(x[0])
bkgTF = ROOT.TF1("bkgTF", fBkg, x_min, x_max, 1)
bkgTF.SetParameter(0,  (1-f)* (histo.GetBinCenter(2) - histo.GetBinCenter(1)) )

# Create a canvas to draw both the histogram and the function
canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)



# Draw the histogram
histo.Draw("E")  # "E" option draws error bars

# Draw the PDF as a TF1 on the same canvas
# pdf_function.Draw("")
new_func.Draw("same")
sigTF.Draw("same")
# bkgTF.Draw("same")

# Update the canvas to display the plot
canvas.Update()

# Optionally, save the plot to a file
canvas.SaveAs("data_pdf_comparison.png")

print (f)
print (sigPdf.Integral(x_min,x_max))
print(pdf_function.Integral(x_min,x_max))
print(bkgPdf.Integral(x_min,x_max))

# print(f"Chi2/ndof signal {chi2}, Chi2/ndof sidebands {chi2Bkg}")
