import argparse
import sys
import numpy as np
import yaml
sys.path.append('../../utils')
from DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error
import os
import pandas as pd
import uproot
import ROOT
import matplotlib.pyplot as plt
from ResoFitUtils import fit_invariant_mass_with_background, roofit_plot_with_matplotlib, simultaneous_fit_signal_background, template_fit_signal_background, fit_bkg_only
from particle import Particle

inputfile = '/home/luca/alice/Ds_reso/AnalysisResults.root'
folder =  'hf-candidate-creator-charm-reso-reduced'
histname = 'hMassDs1'

# Open the ROOT file
file = ROOT.TFile(inputfile, "READ")

# Navigate to the folder (directory) inside the file
dir = file.Get(folder)
if not dir:
    print(f"Error: Could not find directory {folder} in file {inputfile}")
    exit()

# Access the histogram inside the folder
hist = dir.Get(histname)
if not hist:
    print(f"Error: Could not find histogram {histname} in folder {folder}")
    exit()
hist.SetMarkerStyle(ROOT.kFullCircle)
hist.SetMarkerSize(0.5)
hist.SetTitle("DS1 mass")
hist.GetXaxis().SetTitleSize(0.045)
hist.GetXaxis().SetNdivisions(508)
hist.GetYaxis().SetNdivisions(505)
hist.GetXaxis().SetLabelSize(0.045)
hist.GetXaxis().SetDecimals()
hist.GetYaxis().SetTitleSize(0.045)
hist.GetYaxis().SetTitleOffset(1.4)
hist.GetYaxis().SetLabelSize(0.045)
hist.GetYaxis().SetDecimals(1)
hist.GetYaxis().SetRangeUser(1, hist.GetMaximum() * 1.4)
nBins = hist.GetNbinsX()
xmin = hist.GetXaxis().GetXmin()
xmax = hist.GetXaxis().GetXmax()
delta_x = mD = Particle.from_pdgid(413).mass*1e-3
new_hist = hist.ProjectionX()  # Clone the original histogram


# Loop over all the bins in the original histogram and shift the bin centers
# for bin_idx in range(1, hist.GetNbinsX() + 1):
#     # Get the content and error from the original histogram
#     content = hist.GetBinContent(bin_idx)
#     error = hist.GetBinError(bin_idx)
    
#     # Get the bin center and shift it
#     bin_center = hist.GetBinCenter(bin_idx)
#     new_bin_center = bin_center + delta_x
    
#     # Find the corresponding bin in the new histogram
#     new_bin_idx = new_hist.FindBin(new_bin_center)
    
#     # Assign the content and error to the new bin
#     new_hist.SetBinContent(new_bin_idx, content)
#     new_hist.SetBinError(new_bin_idx, error)
workspace, mass_frame = fit_invariant_mass_with_background(new_hist, 2.51, 2.63,10433, False)
canvas = ROOT.TCanvas("canv_masses", "", 500, 500)
Nsig = workspace.var("nSigR").getVal()
mass_frame.Draw()
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # x1, y1, x2, y2 (coordinates in NDC)
legend.SetHeader("Legend", "C")  # Optional header for the legend
legend.AddEntry(mass_frame.findObject("data"), "data", "p")  # Add PDF entry (line)
legend.AddEntry("", f"Yield = {Nsig:.0f}", "")  # Add data entry (points)

# Draw the legend
legend.Draw()

canvas.SaveAs('histDs1.pdf')
