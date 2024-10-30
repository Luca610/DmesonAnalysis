import argparse
import ROOT
from particle import Particle


def set_data_style(hist, had="dzero", isnorm=False):
    """

    """

    hist.SetDirectory(0)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(0.7)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kBlack)
    if isnorm:
        hist.GetYaxis().SetTitle(
            f"Normalised counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
    else:
        hist.GetYaxis().SetTitle(f"Counts per {hist.GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")

    if "ds1" in had:
        hist.GetXaxis().SetTitle("#it{M}(D^{*+}K_{S}^{0}) (GeV/#it{c}^{2})") 
    elif "ds2star" in had:
        hist.GetXaxis().SetTitle("#it{M}(D^{+}K_{S}^{0}) (GeV/#it{c}^{2})")

    hist.SetTitle("")
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


def get_spline_fromhist(hist, func_type="totfunc", norm_fact = 1.):
    """

    """
    hist.Scale(norm_fact)
    spline = ROOT.TSpline3(hist)
    spline.SetLineWidth(3)
    if func_type == "totfunc":
        spline.SetLineColor(ROOT.kBlue)
    elif func_type == "bkg":
        spline.SetLineColor(ROOT.kRed)
        spline.SetLineStyle(2)
    elif func_type == "signal":
        spline.SetLineColor(ROOT.kAzure + 4)
        spline.SetFillColor(ROOT.kAzure + 4)
        spline.SetFillStyle(3145)
    elif func_type == "signalPart":
        spline.SetLineColor(ROOT.kMagenta + 1)
        
    return spline


def get_spline_fromfunc(func, func_type="totfunc", norm_fact = 1., x_min=None, x_max=None):
    """

    """

    if x_min is None:
        x_min = func.GetMaximumX()
    if x_max is None:
        x_max = func.GetMinimumX()

    hist = ROOT.TH1F(f"hist_from_{func.GetName()}", "", 10000, x_min, x_max)
    for ibin in range(1, 10001):
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
    if had == "ds1":
        mR = Particle.from_pdgid(10433).mass*1e-3
        canv = infile.Get("c1")
        hist_mass = canv.GetListOfPrimitives().FindObject('hSparseMl_proj_3')
        norm = 1./hist_mass.Integral()
        set_data_style(hist_mass, had)
        hist_mass_norm = hist_mass.Clone()
        hist_mass_norm.Scale(norm)
        set_data_style(hist_mass_norm, had, True)
        func_totfunc = canv.GetListOfPrimitives().FindObject('totTF')
        func_totfunc = get_spline_fromfunc(func_totfunc, "totfunc", 1., 2.51, 2.95)
        func_totfunc_norm = get_spline_fromfunc(func_totfunc, "totfunc", norm, 2.51, 2.95)
        func_bkg = canv.GetListOfPrimitives().FindObject('bkgTF')
        func_bkg = get_spline_fromfunc(func_bkg, "bkg", 1., 2.51, 2.95)
        func_bkg_norm = get_spline_fromfunc(func_bkg, "bkg", norm, 2.51, 2.95)
        func_sig = canv.GetListOfPrimitives().FindObject('sigTF')
        func_sig = get_spline_fromfunc(func_sig, "signal", 1., mR - 0.012, mR + 0.012)
        func_sig_norm = get_spline_fromfunc(func_sig, "signal", norm, mR - 0.012, mR + 0.012)
    if had == "ds2star":
        mR = Particle.from_pdgid(435).mass*1e-3
        canv = infile.Get("c1")
        hist_mass = canv.GetListOfPrimitives().FindObject('hSparseMl_proj_3')
        norm = 1./hist_mass.Integral()
        set_data_style(hist_mass, had)
        hist_mass_norm = hist_mass.Clone()
        hist_mass_norm.Scale(norm)
        set_data_style(hist_mass_norm, had, True)
        func_totfunc = canv.GetListOfPrimitives().FindObject('totTF')
        func_totfunc = get_spline_fromfunc(func_totfunc, "totfunc", 1., 2.37, 2.73)
        func_totfunc_norm = get_spline_fromfunc(func_totfunc, "totfunc", norm, 2.37, 2.73)
        func_bkg = canv.GetListOfPrimitives().FindObject('bkgTF')
        func_bkg = get_spline_fromfunc(func_bkg, "bkg", 1., 2.37, 2.73)
        func_bkg_norm = get_spline_fromfunc(func_bkg, "bkg", norm, 2.37, 2.73)
        func_sig = canv.GetListOfPrimitives().FindObject('sigTF')
        func_sig = get_spline_fromfunc(func_sig, "signal", 1., mR - 0.04, mR + 0.04)
        func_sig_norm = get_spline_fromfunc(func_sig, "signal", norm, mR - 0.04, mR + 0.04)
        func_sig_part = canv.GetListOfPrimitives().FindObject('partTF')
        func_sig_part = get_spline_fromfunc(func_sig_part, "signalPart", 1., 2.37, 2.73)
        func_sig_part_norm = get_spline_fromfunc(func_sig_part, "signalPart", norm, 2.37, 2.73)
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
    # elif "ds" in had:
    #     canv = infile.Get("canvas_Pt8_12")
    #     hist_mass = canv.GetListOfPrimitives().FindObject('fHistoInvMass')
    #     norm = 1./hist_mass.Integral()
    #     set_data_style(hist_mass, had)
    #     hist_mass_norm = hist_mass.Clone()
    #     hist_mass_norm.Scale(norm)
    #     set_data_style(hist_mass_norm, had, True)
    #     func_totfunc = canv.GetListOfPrimitives().FindObject('funcmass')
    #     func_totfunc = get_spline_fromfunc(func_totfunc, "totfunc", 1., 1.7, 2.1)
    #     func_totfunc_norm = get_spline_fromfunc(func_totfunc, "totfunc", norm, 1.7, 2.1)
    #     func_bkg = canv.GetListOfPrimitives().FindObject('funcbkgrefit')
    #     func_bkg = get_spline_fromfunc(func_bkg, "bkg", 1., 1.7, 2.1)
    #     func_bkg_norm = get_spline_fromfunc(func_bkg, "bkg", norm, 1.7, 2.1)
                

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
        "ds1": "2 < #it{p}_{T} < 24 GeV/#it{c}",
        "ds2star": "5 < #it{p}_{T} < 24 GeV/#it{c}", 
        "dzerohighpt": "7 < #it{p}_{T} < 8 GeV/#it{c}",
        "dplus": "3 < #it{p}_{T} < 4 GeV/#it{c}",
        "dplushighpt": "6 < #it{p}_{T} < 7 GeV/#it{c}",
        "ds": "6 < #it{p}_{T} < 8 GeV/#it{c}",
        "lc": "6 < #it{p}_{T} < 8 GeV/#it{c}"
    }

    labels = {
        "ds1" : "D_{s1}^{+} #rightarrow D*^{+}K_{S}^{0} and charge conj.",
        "ds2star" : "D*_{s2}^{+} #rightarrow D^{+}K_{S}^{0} and charge conj.",
        "dzero": "D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.",
        "dzerohighpt": "D^{0} #rightarrow K^{#minus}#pi^{+} and charge conj.",
        "dplus": "D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+} and charge conj.",
        "dplushighpt": "D^{+} #rightarrow K^{#minus}#pi^{+}#pi^{+} and charge conj.",
        "ds": "D_{s}^{+} #rightarrow #phi#pi^{+} #rightarrow K^{+}K^{#font[122]{-}}#pi^{+} and charge conj.",
        "lc": "#Lambda_{c}^{+} #rightarrow pK^{#minus}#pi^{+} and charge conj.",
    }
    if had == "ds2star":
        leg = ROOT.TLegend(0.33, 0.25, 0.56, 0.5)
    else:
        leg = ROOT.TLegend(0.4, 0.25, 0.63, 0.5) 
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    #leg.SetHeader(labels[had])
    leg.AddEntry(func_sig, labels[had], "f")
    leg.AddEntry(hist_mass, "Data", "p")
    leg.AddEntry(func_bkg, "Background", "l")
    leg.AddEntry(func_totfunc, "Total fit", "l")
    if had == "ds2star":
        leg.AddEntry(func_sig_part, "Partly reconstructed D_{s1}^{+} #rightarrow D*^{+}K_{S}^{0}", "l")

    canv_masses = ROOT.TCanvas("canv_masses", "", 500, 500)
    if "ds1" in had:
        hist_mass.GetXaxis().SetRangeUser(2.51, 2.72)
    hist_mass.DrawCopy("e")
    func_bkg.Draw("lsame")
    func_totfunc.SetNpx(10000)
    func_totfunc.Draw("lsame")
    func_sig.SetNpx(10000)
    func_sig.Draw("FCsameR")
    if had == "ds2star":
        func_sig_part.Draw("FCsame")
    lat.DrawLatex(0.17, 0.86, "ALICE Performance")
    lat_small.DrawLatex(0.17, 0.8, "pp, #sqrt{#it{s}} = 13.6 TeV, #it{L}_{int} = 11 pb^{-1}")
    lat_small.DrawLatex(0.17, 0.74, pt_labels[had])
    # lat_small.DrawLatex(0.3, 0.68, labels[had])
    leg.Draw()

    canv_masses.SaveAs(f"{had}_massfit.pdf")
    canv_masses.SaveAs(f"{had}_massfit.png")
    canv_masses.SaveAs(f"{had}_massfit.eps")

    canv_masses_norm = ROOT.TCanvas("canv_masses_norm", "", 500, 500)
    if "ds1" in had:
        hist_mass.GetXaxis().SetRangeUser(2.51, 2.7)
        canv_masses_norm.DrawFrame(2.5, 0, 2.7, 0.02, \
        ';#it{M}(D*^{+}K_{S}^{0}) (GeV/#it{c}^{2}); Normalised counts per 2 MeV/#it{c}^{2}')
    if "ds2star" in had:
        canv_masses_norm.DrawFrame(2.37, 0, 2.73, 0.009, \
        ';#it{M}(D*^{+}K_{S}^{0}) (GeV/#it{c}^{2}); Normalised counts per 2 MeV/#it{c}^{2}')
    hist_mass_norm.DrawCopy("esame")
    func_bkg_norm.Draw("lsame")
    func_totfunc_norm.SetNpx(10000)
    func_totfunc_norm.Draw("lsame")
    func_sig_norm.Draw("FCsame")
    lat.DrawLatex(0.17, 0.86, "ALICE Performance")
    lat_small.DrawLatex(0.17, 0.8, "pp, #sqrt{#it{s}} = 13.6 TeV, #it{L}_{int} = 11 pb^{-1}")
    lat_small.DrawLatex(0.17, 0.74, pt_labels[had])
    # lat_small.DrawLatex(0.17, 0.68, labels[had])
    # mean = canv.GetListOfPrimitives().FindObject('funcmass').GetParameter(3)
    # meanerr = canv.GetListOfPrimitives().FindObject('funcmass').GetParError(3)
    # sigma = canv.GetListOfPrimitives().FindObject('funcmass').GetParameter(4)
    # sigmaerr = canv.GetListOfPrimitives().FindObject('funcmass').GetParError(4)
    #lat_smaller = ROOT.TLatex()
    #lat_smaller.SetNDC()
    #lat_smaller.SetTextFont(42)
    #lat_smaller.SetTextColor(ROOT.kBlack)
    #lat_smaller.SetTextSize(0.04)
    #lat_smaller.DrawLatex(0.16, 0.15, f'#mu = ({mean*1000:.1f} #pm {meanerr*1000:.1f}) MeV/#it{{c}}^{{2}}')
    #lat_smaller.DrawLatex(0.16, 0.20, f'#sigma = ({sigma*1000:.1f} #pm {sigmaerr*1000:.1f}) MeV/#it{{c}}^{{2}}')

    leg.Draw()

    canv_masses_norm.SaveAs(f"{had}_massfit_norm.pdf")
    canv_masses_norm.SaveAs(f"{had}_massfit_norm.eps")
    canv_masses_norm.SaveAs(f"{had}_massfit_norm.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--input_file", "-i", metavar="text",
                        default="mass.root", help="input root file", required=True)
    parser.add_argument("--particle", "-p", metavar="text",
                        default="ds1", help="particle species", required=False)
    args = parser.parse_args()

    plot(args.input_file, args.particle)