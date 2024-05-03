// purpose: find out how well pt weights work with resolution effect taken into account

#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2023.h"
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

TF1 *getMesonEfficiency(std::string fname, Double_t theXmax = 10.)
{

    TH1 *hEffi = (TH1 *)getObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dP/dptG", 1, 0., theXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();

    auto leg2 = new TLegend();

    auto drawAndAddToLegendF = [leg2](TF1 *f, EColor c = kBlue)
    {f->SetLineColor(c); f->Draw("same"); leg2->AddEntry(f,f->GetName()); };

    hEffi->Draw("same");
    leg2->Draw("same");

    TF1 *fEffi = new TF1("fEffi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., theXmax);
    fEffi->SetParameters(1, 1, -1, 0.);

    hEffi->Fit(fEffi, "N");
    drawAndAddToLegendF(fEffi, kRed);

    /*
    MyDerivative *fptr = new MyDerivative(fEffi, 1);  // create the user function class
    auto fEffi_d1 = fptr->GetD();

    MyDerivative *fptr2 = new MyDerivative(fEffi, 2);  // create the user function class
    auto fEffi_d2 = fptr2->GetD();

    drawAndAddToLegendF(fEffi_d1);
    drawAndAddToLegendF(fEffi_d2, kGreen);*/

    return fEffi;
}

void investigatePtWeights_wResolutionEffects()
{

    gROOT->Reset();

    //~ gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendTextSize(.03);

    std::string meson("Pi0");
    std::string cent("10130053");
    std::string centAS("10130023");
    std::string fnameAS("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root");
    // std::string fnameAS("/trains/2024-01-26_LHC24a1_QA_noPtW/GCo_997_both.root");

    std::string fnameWeightsFile(Form(
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root"));

    GCo gAS({fnameAS, "GammaConvV1_997/", centAS, "_0d200009ab770c00amd0404000_0152101500000000"});
    //~ TH2* h2Resolution = (TH2*)gMB.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt");
    TH2 &h2Resolution = *(TH2 *)gAS.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt");

    int ptBinStart = 1;
    int ptBinMax = 31;
    
    // get some parametrized effi
    TF1 *fEffiAtAll_dp_dptG = getMesonEfficiency("input_for_effi-fit_101.root");

    std::string lFnameResFits(
        Form("%s_resolutionFits_%d-%d.root", centAS.data(), ptBinStart, ptBinMax));

    // create resolution parametrizations
    TPairFitsWAxis lPair_vFits_ptG_i_dp_dr_Axis =
        computeResolutionFits(h2Resolution,
                              ptBinStart,
                              ptBinMax,
                              4 /*nRebin_r*/,
                              lFnameResFits,
                              true /*plotSingles*/);

    //~ drawAllFitsOnTop(lPair_vFits_ptG_i_dp_dr_Axis.first);
    TAxis &lPtGaxis = lPair_vFits_ptG_i_dp_dr_Axis.second;

    // get Infos for Weights instance
    TH1 *hGenDist_AS_inv = (TH1 *)getObjectFromPathInFile(
        fnameWeightsFile, meson + "_LHC24a1_5TeV_" + centAS);
    TF1 *fGenData_dn_dptG_inv = (TF1 *)getObjectFromPathInFile(
        fnameWeightsFile, std::string(meson + "_Data_5TeV_" + cent.substr(0, 6)));

    PtWeights *lPtWeights = new PtWeights("lPtWeights",
                                          *hGenDist_AS_inv,
                                          *fGenData_dn_dptG_inv,
                                          lPtGaxis);

    // h_inv -> f_inv -> f
    auto getGenDist = [&hGenDist_AS_inv, &lPtGaxis]()
    {
        TF1 *lGenDist_AS_dn_dptG_inv = new TF1("lGenDist_AS_dn_dptG_inv", "[0] + [1]/(x-[2])", lPtGaxis.GetXmin(), lPtGaxis.GetXmax());
        hGenDist_AS_inv->Fit(lGenDist_AS_dn_dptG_inv, "N");
        return multiplyTF1ByX(*lGenDist_AS_dn_dptG_inv, "lGenDist_AS_dn_dptG");
    };
    TF1 *lGenDist_AS_dn_dptG = getGenDist();

    // more accurate way to get the genDist:
    // h_inv -> h -> f
    TH1 &hGenDist_AS_dn_dptG = *multiplyTH1ByBinCenters(*hGenDist_AS_inv, "", "hGenDist_AS_dn_dptG");
     
    auto getGenDist_natural = [&hGenDist_AS_inv, &lPtGaxis]()
    {
    };

    // TAxis &axisPtR = lPtGaxis;
    TAxis lAxisPtR(100, 0., 10.);

    // create MCEffi instances
    auto &lMCEffi_AS = *new MCEffi_wRes("lMCEffi_AS",                 //
                                        *lGenDist_AS_dn_dptG,         // _fGenDist_dn_dptG
                                        *fEffiAtAll_dp_dptG,          // _fEffi_dp_dptG
                                        lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                        lAxisPtR,                     // _axisPtR
                                        lPtWeights);

    auto &lMCEffi_D = *new MCEffi_wRes("lMCEffi_D", //
                                       *fGenData_dn_dptG_inv, // _fGenDist_dn_dptG
                                       *fEffiAtAll_dp_dptG,   // _fEffi_dp_dptG
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
