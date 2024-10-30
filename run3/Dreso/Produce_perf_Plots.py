"""

"""
import argparse
import ROOT

def set_data_style(hist, had="dzero", isnorm=False):
    """

    """

    hist.SetDirectory(0)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlack)
    if isnorm:
        hist.GetYaxis().SetTitle(
            f"Normalised counts / {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
    else:
        hist.GetYaxis().SetTitle(f"Counts / {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")

    if "dzero" in had:
        hist.GetXaxis().SetTitle("#it{M}(K#pi) (GeV/#it{c}^{2})")
    elif "dplus" in had:
        hist.GetXaxis().SetTitle("#it{M}(K#pi#pi) (GeV/#it{c}^{2})")
    elif "ds" in had:
        hist.GetXaxis().SetTitle("#it{M}(KK#pi) (GeV/#it{c}^{2})")
    elif had == "lc":
        hist.GetXaxis().SetTitle("#it{M}(pK#pi) (GeV/#it{c}^{2})")

    hist.SetTitle("")
    hist.GetXaxis().SetTitleSize(0.045)
    hist.GetXaxis().SetNdivisions(508)
    hist.GetYaxis().SetNdivisions(505)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.4)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetDecimals(1)
    if "dzero" not in had:
        if "ds" in had and isnorm:
            hist.GetYaxis().SetRangeUser(hist.GetMinimum() / 60.0, hist.GetMaximum() * 1.5)
        elif "ds" in had and not isnorm:
            hist.GetYaxis().SetRangeUser(-5, hist.GetMaximum() * 1.5)
        else:
            hist.GetYaxis().SetRangeUser(hist.GetMinimum() / 1.2, hist.GetMaximum() * 1.2)
    else:
        hist.GetYaxis().SetRangeUser(0., hist.GetMaximum() * 1.3)


def get_spline_fromhist(hist, func_type="totfunc", norm_fact = 1.):
    """

    """
    hist.Scale(norm_fact)
    spline = ROOT.TSpline3(hist)
    spline.SetLineWidth(2)
    if func_type == "totfunc":
        spline.SetLineColor(ROOT.kBlue)
    elif func_type == "bkg":
        spline.SetLineColor(ROOT.kRed)
        spline.SetLineStyle(2)
    elif func_type == "signal":
        spline.SetLineColor(ROOT.kAzure + 4)

    return spline


def get_spline_fromfunc(func, func_type="totfunc", norm_fact = 1., x_min=None, x_max=None):
    """

    """

    if x_min is None:
        x_min = func.GetMaximumX()
    if x_max is None:
        x_max = func.GetMinimumX()

    hist = ROOT.TH1F(f"hist_from_{func.GetName()}", "", 1000, x_min, x_max)
    for ibin in range(1, 1001):
        hist.SetBinContent(ibin, func.Eval(hist.GetBinCenter(ibin)))

    return get_spline_fromhist(hist, func_type, norm_fact)


def plot(infile_name, had):
    """

    """

    ROOT.gStyle.SetPadRightMargin(0.035)
    ROOT.gStyle.SetPadTopMargin(0.065)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.TGaxis.SetMaxDigits(3)

    infile = ROOT.TFile.Open(infile_name)
    hist_mass, hist_mass_norm = None, None
    func_bkg, func_totfunc, func_bkg_norm, func_totfunc_norm = None, None, None, None
    if had == "lc":
        hist_mass = infile.Get("hdata_lc_pt6_8")
        set_data_style(hist_mass, had)
        norm = 1./hist_mass.Integral()
        hist_mass_norm = hist_mass.Clone()
        hist_mass_norm.Scale(norm)
        set_data_style(hist_mass_norm, had, True)
        func_bkg = get_spline_fromhist(infile.Get("bkg_0_lc_pt6_8"), "bkg")
        func_totfunc = get_spline_fromhist(infile.Get("total_func_lc_pt6_8"), "totfunc")
        func_bkg_norm = get_spline_fromhist(infile.Get("bkg_0_lc_pt6_8"), "bkg", norm)
        func_totfunc_norm = get_spline_fromhist(infile.Get("total_func_lc_pt6_8"), "totfunc", norm)
    elif "dplus" in had:
        canv = infile.Get("mass spectra")
        ipt = 6 if had == "dplushighpt" else 1
        for obj in canv.GetPad(ipt).GetListOfPrimitives():
            obj_name = obj.GetName()
            if obj_name == "fhistoInvMass":
                hist_mass = obj
                hist_mass.Scale(1.)
                norm = 1./hist_mass.Integral()
                set_data_style(hist_mass, had)
                hist_mass_norm = hist_mass.Clone()
                hist_mass_norm.Scale(norm)
                set_data_style(hist_mass_norm, had, True)
                func_bkg = get_spline_fromfunc(obj.GetFunction("funcbkgRecalc"), "bkg")
                func_bkg_norm = get_spline_fromfunc(obj.GetFunction("funcbkgRecalc"), "bkg", norm)
                obj.GetFunction("funcbkgRecalc").SetLineColor(ROOT.kWhite)
                obj.GetFunction("funcbkgFullRange").SetLineColor(ROOT.kWhite)
                obj.GetFunction("funcbkgRecalc").SetLineWidth(0)
                obj.GetFunction("funcbkgFullRange").SetLineWidth(0)
            elif obj_name == "funcmass":
                func_totfunc = get_spline_fromfunc(obj, "totfunc", 1., 1.71, 2.)
                func_totfunc_norm = get_spline_fromfunc(obj, "totfunc", norm, 1.71, 2.)
    elif "dzero" in had:
        canv = infile.Get("cout")
        for obj in canv.GetListOfPrimitives():
            obj_name = obj.GetName()
            if "hMass_bins" in obj_name:
                hist_mass = obj
                hist_mass.Scale(1.)
                norm = 1./hist_mass.Integral()
                set_data_style(hist_mass, had)
                hist_mass_norm = hist_mass.Clone()
                hist_mass_norm.Scale(norm)
                set_data_style(hist_mass_norm, had, True)
            elif obj_name == "func":
                func_totfunc = get_spline_fromfunc(obj, "totfunc", 1., 1.7, 2.05)
                func_totfunc_norm = get_spline_fromfunc(obj, "totfunc", norm, 1.7, 2.05)
            elif obj_name == "fbkg":
                func_bkg = get_spline_fromfunc(obj, "bkg", 1., 1.7, 2.05)
                func_bkg_norm = get_spline_fromfunc(obj, "bkg", norm, 1.7, 2.05)
    elif "ds" in had:
        canv = infile.Get("canvas_Pt8_12")
        hist_mass = canv.GetListOfPrimitives().FindObject('fHistoInvMass')
        norm = 1./hist_mass.Integral()
        set_data_style(hist_mass, had)
        hist_mass_norm = hist_mass.Clone()
        hist_mass_norm.Scale(norm)
        set_data_style(hist_mass_norm, had, True)
        func_totfunc = canv.GetListOfPrimitives().FindObject('funcmass')
        func_totfunc = get_spline_fromfunc(func_totfunc, "totfunc", 1., 1.7, 2.1)
        func_totfunc_norm = get_spline_fromfunc(func_totfunc, "totfunc", norm, 1.7, 2.1)
        func_bkg = canv.GetListOfPrimitives().FindObject('funcbkgrefit')
        func_bkg = get_spline_fromfunc(func_bkg, "bkg", 1., 1.7, 2.1)
        func_bkg_norm = get_spline_fromfunc(func_bkg, "bkg", norm, 1.7, 2.1)
                

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextSize(0.055)

    lat_small = ROOT.TLatex()
    lat_small.SetNDC()
    lat_small.SetTextFont(42)
    lat_small.SetTextColor(ROOT.kBlack)
    lat_small.SetTextSize(0.045)

    pt_labels = {
        "dzero": "0 < #it{p}_{T} < 1 GeV/#it{c}",
        "dzerohighpt": "7 < #it{p}_{T} < 8 GeV/#it{c}",
        "dplus": "3 < #it{p}_{T} < 4 GeV/#it{c}",
        "dplushighpt": "6 < #it{p}_{T} < 7 GeV/#it{c}",
        "ds": "6 < #it{p}_{T} < 8 GeV/#it{c}",
        "lc": "6 < #it{p}_{T} < 8 GeV/#it{c}"
    }

    labels = {
        "dzero": "D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.",
        "dzerohighpt": "D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.",
        "dplus": "D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+} and charge conj.",
        "dplushighpt": "D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+} and charge conj.",
        "ds": "D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#font[122]{-}}#pi^{+} and charge conj.",
        "lc": "#Lambda_{c}^{+} #rightarrow pK^{#minus}#pi^{+} and charge conj.",
    }

    leg = ROOT.TLegend(0.17, 0.49, 0.4, 0.64)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    #leg.SetHeader(labels[had])
    leg.AddEntry(hist_mass, "Data", "p")
    leg.AddEntry(func_bkg, "Background", "l")
    leg.AddEntry(func_totfunc, "Total fit function", "l")

    canv_masses = ROOT.TCanvas("canv_masses", "", 500, 500)
    if "dplus" in had:
        hist_mass.GetXaxis().SetRangeUser(1.75, 1.98)
    hist_mass.DrawCopy("e")
    func_bkg.Draw("lsame")
    func_totfunc.Draw("lsame")
    lat.DrawLatex(0.3, 0.86, "ALICE Performance")
    lat_small.DrawLatex(0.3, 0.8, "10#font[122]{-}30% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
    lat_small.DrawLatex(0.3, 0.74, pt_labels[had])
    lat_small.DrawLatex(0.3, 0.68, labels[had])
    leg.Draw()

    canv_masses.SaveAs(f"/home/fchinu/Run3/Ds_PbPb_5TeV/CutVariation/SignificanceScans/Figures/{had}_massfit_10_30_6_8.pdf")
    canv_masses.SaveAs(f"/home/fchinu/Run3/Ds_PbPb_5TeV/CutVariation/SignificanceScans/Figures/{had}_massfit_10_30_6_8.eps")

    canv_masses_norm = ROOT.TCanvas("canv_masses_norm", "", 500, 500)
    #if "dplus" in had:
    #    hist_mass_norm.GetXaxis().SetRangeUser(1.75, 1.98)
    canv_masses_norm.DrawFrame(1.71, 0, 2.1, 0.075, \
        f';#it{{M}}(KK#pi) (GeV/#it{{c}}^{{2}}); Normalised counts / {hist_mass.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}')
    hist_mass_norm.DrawCopy("esame")
    func_bkg_norm.Draw("lsame")
    func_totfunc_norm.Draw("lsame")
    lat.DrawLatex(0.17, 0.86, "This Thesis")
    lat_small.DrawLatex(0.17, 0.8, "10#font[122]{-}30% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.36 TeV")
    lat_small.DrawLatex(0.17, 0.74, pt_labels[had])
    lat_small.DrawLatex(0.17, 0.68, labels[had])
    mean = canv.GetListOfPrimitives().FindObject('funcmass').GetParameter(3)
    meanerr = canv.GetListOfPrimitives().FindObject('funcmass').GetParError(3)
    sigma = canv.GetListOfPrimitives().FindObject('funcmass').GetParameter(4)
    sigmaerr = canv.GetListOfPrimitives().FindObject('funcmass').GetParError(4)
    #lat_smaller = ROOT.TLatex()
    #lat_smaller.SetNDC()
    #lat_smaller.SetTextFont(42)
    #lat_smaller.SetTextColor(ROOT.kBlack)
    #lat_smaller.SetTextSize(0.04)
    #lat_smaller.DrawLatex(0.16, 0.15, f'#mu = ({mean*1000:.1f} #pm {meanerr*1000:.1f}) MeV/#it{{c}}^{{2}}')
    #lat_smaller.DrawLatex(0.16, 0.20, f'#sigma = ({sigma*1000:.1f} #pm {sigmaerr*1000:.1f}) MeV/#it{{c}}^{{2}}')

    leg.Draw()

    canv_masses_norm.SaveAs(f"/home/fchinu/Run3/Ds_PbPb_5TeV/CutVariation/SignificanceScans/Figures/{had}_massfit_norm_10_30_6_8.pdf")
    canv_masses_norm.SaveAs(f"/home/fchinu/Run3/Ds_PbPb_5TeV/CutVariation/SignificanceScans/Figures/{had}_massfit_norm_10_30_6_8.eps")
    canv_masses_norm.SaveAs(f"/home/fchinu/Run3/Ds_PbPb_5TeV/CutVariation/SignificanceScans/Figures/{had}_massfit_norm_10_30_6_8.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="mass.root", help="input root file", required=True)
    parser.add_argument("--particle", "-p", metavar="text",
                        default="dzero", help="particle species", required=False)
    args = parser.parse_args()

    plot(args.input_file, args.particle)
    input("Press ENTER to exit")
