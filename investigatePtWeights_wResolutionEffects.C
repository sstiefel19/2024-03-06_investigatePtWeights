// purpose: find out how well pt weights work with resolution effect taken into account
#include "/analysisSoftware/SupportingMacros/GCo.h"
#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2024.h"
#include "/analysisSoftware/SupportingMacros/utils_TF1.h"

#include "computeResolutionFits.h"
#include "source/MCEffi.h"
#include "source/PtWeights.h"

#include <iostream>
#include <string.h>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"

TF1 &getMesonEfficiency(std::string fname, Double_t theXmax = 10.)
{

    TH1 &hEffi = *(TH1*)getObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dP/dptG", 1, 0., theXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();
    auto leg2 = new TLegend();
    auto drawAndAddToLegendF = [leg2](TF1 *f, EColor c = kBlue)
    {f->SetLineColor(c); f->Draw("same"); leg2->AddEntry(f,f->GetName()); };

    hEffi.Draw("same");
    leg2->Draw("same");

    TF1 &fEffi = *new TF1("fEffi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., theXmax);
    fEffi.SetParameters(1, 1, -1, 0.);

    hEffi.Fit(&fEffi, "N");
    drawAndAddToLegendF(&fEffi, kRed);
    return fEffi;
}

// set first parameter to "auto" for auto-concatenated name
TF1 &createGenDistFit(std::string const &theResultNameInfo,
                      TH1 &theTH1GenDist_dn_dptG,
                      bool theTH1IsInvariant,
                      TAxis const &thePtGaxis,
                      bool theMultiplyResultTF1ByX = false)
{
    // on top for readability
    bool lResultWillBeInvariant = theTH1IsInvariant && !theMultiplyResultTF1ByX;

    // must not be empty
    if (!theResultNameInfo.size())
    {
        printf("investigatePtWeights_wResolutionEffects::createGenDistFit():\n\t"
               "ERROR: theResultNameInfo can't be empty!\n"
               "\tChose explicit name or \'auto\' for auto-concatenated name (verbose).\n"
               "\tReturning dummy TF1.\n");
        return utils_TF1::GetDummyTF1(theTH1GenDist_dn_dptG.GetName(), true /*theAddDummyTag*/);
    }

    // compile result name
    std::string lResultFullName(
        theResultNameInfo != "auto" ? theResultNameInfo
                                    : Form("createGenDistFit_dn_dptG_from_%s_%s",
                                           theTH1GenDist_dn_dptG.GetName(),
                                           lResultWillBeInvariant ? "_inv" : ""));

    // create result TF1
    TF1 &lTF1GenDist_dn_dptG = *new TF1(lResultFullName.data(),
                                        theTH1IsInvariant ? "[0] + [1]/(x-[2])"
                                                          : "[0]",
                                        thePtGaxis.GetXmin(),
                                        thePtGaxis.GetXmax());
    // perform fit
    theTH1GenDist_dn_dptG.Fit(&lTF1GenDist_dn_dptG, "N");

    return (!theMultiplyResultTF1ByX)
               ? lTF1GenDist_dn_dptG
               : utils_TF1::MultiplyTF1ByX(lResultFullName,
                                           lTF1GenDist_dn_dptG);
}

void investigatePtWeights_wResolutionEffects()
{
    gROOT->Reset();

    //~ gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    // 0) definitions
    int ptBinStart = 1;
    int ptBinMax = 31;

    std::string meson("Pi0");
    std::string cent("10130053");
    std::string centAS("10130023");
    std::string fnameAS("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root");
    // std::string fnameAS("/trains/2024-01-26_LHC24a1_QA_noPtW/GCo_997_both.root");

    std::string lFnameInputEffiFit("input_for_effi-fit_101.root");
    std::string lFnameResFits(
        Form("%s_resolutionFits_%d-%d.root", centAS.data(), ptBinStart, ptBinMax));

    std::string fnameWeightsFile(Form(
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root"));

    GCo gAS({fnameAS, "GammaConvV1_997/", centAS, "_0d200009ab770c00amd0404000_0152101500000000"});
    TH2 &h2Resolution = *(TH2 *)gAS.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt");
    //~ TH2* h2Resolution = (TH2*)gMB.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt");

    auto setupWeightsInstance =
        [&meson, &cent, &centAS, &fnameWeightsFile](std::string const &theID,
                                                    TAxis const &thePtGaxis)
    {
        // get Infos for Weights instance
        TH1 &hGenDist_AS_inv = *(TH1 *)getObjectFromPathInFile(
            fnameWeightsFile, meson + "_LHC24a1_5TeV_" + centAS);
        TF1 &fGenData_dn_dptG_inv = *(TF1 *)getObjectFromPathInFile(
            fnameWeightsFile, meson + "_Data_5TeV_" + cent.substr(0, 6));

        return new PtWeights(theID,
                             hGenDist_AS_inv,
                             fGenData_dn_dptG_inv,
                             thePtGaxis);
    };

    // 2) get some parametrized effi
    TF1 &fEffiAtAll_dp_dptG = getMesonEfficiency(lFnameInputEffiFit);

    // 3) create resolution parametrizations
    TPairFitsWAxis lPair_vFits_ptG_i_dp_dr_Axis =
        computeResolutionFits(h2Resolution,
                              ptBinStart,
                              ptBinMax,
                              4 /*nRebin_r*/,
                              lFnameResFits,
                              true /*drawAllFitsOverlayed*/,
                              true /*plotSingles*/);

    TAxis &lPtGaxis = lPair_vFits_ptG_i_dp_dr_Axis.second;

    /* ptWeights from not inv form:
        0) get h_inv from weights file
        0.5) h_inv -> h -> f
        1) assume function f for mc generated particles
        2) sample a histo h from f

        3)calculate weights in the "normal" distributions, not the invariant one
        */

    PtWeights &lPtWeights = *setupWeightsInstance("lPtWeights",
                                                  lPtGaxis);

    // make a copy since we need to fit the histogram which stupidly changes the histogram itself
    TH1 &hGenDist_AS_inv = *cloneTH1(lPtWeights.GetTH1MCGen_dn_dptG());

    TF1 &lGenDist_AS_dn_dptG = createGenDistFit("auto",
                                                hGenDist_AS_inv,
                                                true /*theTH1IsInvariant*/,
                                                lPtGaxis,
                                                true /*theMultiplyResultTF1ByX*/);

    // more accurate way to get the genDist:
    // h_inv -> h -> f
    // TH1 &hGenDist_AS_dn_dptG = *multiplyTH1ByBinCenters(*hGenDist_AS_inv, "", "hGenDist_AS_dn_dptG");

    // TAxis &axisPtR = lPtGaxis;
    TAxis lAxisPtR(100, 0., 10.);

    // create MCEffi instances
    auto &lMCEffi_AS = *new MCEffi_wRes("lMCEffi_AS",                 //
                                        lGenDist_AS_dn_dptG,          // _fGenDist_dn_dptG
                                        fEffiAtAll_dp_dptG,           // _fEffi_dp_dptG
                                        lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                        lAxisPtR,                     // _axisPtR
                                        &lPtWeights);

    TF1 &fGenData_dn_dptG = lPtWeights.GetTF1TrgtDist_dn_dptG();
    auto &lMCEffi_D = *new MCEffi_wRes("lMCEffi_D",                  //
                                       fGenData_dn_dptG,         // _fGenDist_dn_dptG
                                       fEffiAtAll_dp_dptG,           // _fEffi_dp_dptG
                                       lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                       lAxisPtR);

    lMCEffi_AS.SampleMeasuredEffi_NW_2();
    lMCEffi_AS.PlotAll();
    lMCEffi_D.PlotAll();

    //~ auto &f2 = *lMCEffi_AS.fdN_dptR_NW.GetIntegrand().GetTF2();

    //~ double int2d = f2.Integral(.1, .2, 0., 9.9);
    //~ cout << int2d << endl;
    //~ return;

    //~ auto &f = *lMCEffi_AS.SampleMeasuredEffi_NW_2(kGreen);

    //~ auto &f = lMCEffi_AS.fdN_dptR_NW.GetTF1();

    //~ for (int i=1; i < axisPtR.GetNbins(); ++i){
    //~ double l = axisPtR.GetBinLowEdge(i);
    //~ double c = axisPtR.GetBinCenter(i);
    //~ double u = axisPtR.GetBinUpEdge(i);
    //~ printf("%s.Eval(bin %d) = %e %e %e\n", f.GetName(), i, f.Eval(l), f.Eval(c), f.Eval(u));
    //~ printf("%s.Integral() in bin %d = %f\n", f.GetName(), i, f.Integral(axisPtR.GetBinLowEdge(i), axisPtR.GetBinUpEdge(i)));
    //~ }

    // compareMeasuredEffis_TF1("", lMCEffi_AS, lMCEffi_D);
    //~ compareMeasuredEffis_TH1("", lMCEffi_AS, lMCEffi_D);
    // compareMeasuredEffis_TH1_New("", lMCEffi_AS, lMCEffi_D);

    // lMCEffi_AS.PlotAll();
    // lMCEffi_D.PlotAll();
}
