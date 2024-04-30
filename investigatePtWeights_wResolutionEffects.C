// purpose: find out how well pt weights work with resolution effect taken into account

#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2023.h"
#include "computeResolutionFits.h"
#include "MCEffi.h"
#include "PtWeights.h"

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

// TCanvas &
// compareMeasuredEffis_TF1(std::string const &fID,
//                          MCEffi_wRes &lMCEffi_AS,
//                          MCEffi_wRes &lMCEffi_D)
// {

//     auto &fEffi_AS_NW = lMCEffi_AS.MeasuredEffiTF1_NW(kBlue);
//     auto &fEffi_AS_WW = *lMCEffi_AS.MeasuredEffiTF1_WW(kGreen);
//     auto &fEffi_D_NW = lMCEffi_D.MeasuredEffiTF1_NW(kRed);

//     auto &fEffi_AS_NW_over_D = *getTF1Division(fID + "_fEffi_AS_NW_over_D", &fEffi_AS_NW, &fEffi_D_NW);
//     auto &fEffi_AS_WW_over_D = *getTF1Division(fID + "_fEffi_AS_WW_over_D", &fEffi_AS_WW, &fEffi_D_NW);

//     auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
//     {
//         o.Draw("same");
//         if (leg)
//             leg->AddEntry(&o, o.GetName(), "l");
//     };
//     // setup canvas
//     TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TF1").data(),
//                                (fID + "compareMeasuredEffis_TF1").data(), 2000, 1000);

//     c1.Divide(1, 2);
//     c1.cd(1);
//     gPad->SetLogy();
//     TH2 *hd21 = new TH2F("hd121", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 4e-3);
//     hd21->Draw();

//     auto leg = new TLegend(.73, .64, .90, .90, "");
//     drawAndAdd(fEffi_AS_NW, leg);
//     drawAndAdd(fEffi_AS_WW, leg);
//     drawAndAdd(fEffi_D_NW, leg);
//     leg->Draw("same");

//     c1.cd(2);
//     fEffi_AS_WW_over_D.Draw();
//     // fEffi_AS_WW_over_D.Dump();

//     //~ gPad->SetLogy();
//     TH2 *hd22 = new TH2F("hd122", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
//     hd22->Draw();

//     auto leg2 = new TLegend();
//     drawAndAdd(fEffi_AS_NW_over_D);
//     drawAndAdd(fEffi_AS_WW_over_D);
//     leg2->Draw("same");

//     saveCanvasAs(c1);
//     return c1;
// }

// TCanvas &
// compareMeasuredEffis_TH1(std::string const &fID,
//                          MCEffi_wRes &lMCEffi_AS,
//                          MCEffi_wRes &lMCEffi_D)
// {

//     auto &hEffi_AS_NW = lMCEffi_AS.SampleMeasuredEffi_NW(kBlue);
//     auto &hEffi_AS_WW = *lMCEffi_AS.SampleMeasuredEffi_WW(kGreen);
//     auto &hEffi_D_NW = lMCEffi_D.SampleMeasuredEffi_NW(kRed);

//     auto &hEffi_AS_NW_over_D = *divideTH1ByTH1(hEffi_AS_NW, hEffi_D_NW, "", "hEffi_AS_NW_over_D");
//     auto &hEffi_AS_WW_over_D = *divideTH1ByTH1(hEffi_AS_WW, hEffi_D_NW, "", "hEffi_AS_WW_over_D");

//     auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
//     {
//         o.Draw("same");
//         if (leg)
//             leg->AddEntry(&o, o.GetName(), "l");
//     };
//     // setup canvas
//     TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TH1").data(),
//                                (fID + "compareMeasuredEffis_TH1").data(), 2000, 1000);

//     c1.Divide(1, 2);
//     c1.cd(1);
//     gPad->SetLogy();
//     TH2 *hd21 = new TH2F("hd21", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 8e-3);
//     hd21->Draw();

//     auto leg = new TLegend(.73, .64, .90, .90, "");
//     drawAndAdd(hEffi_AS_NW, leg);
//     drawAndAdd(hEffi_AS_WW, leg);
//     drawAndAdd(hEffi_D_NW, leg);
//     leg->Draw("same");

//     c1.cd(2);
//     //~ gPad->SetLogy();
//     TH2 *hd22 = new TH2F("hd22", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
//     hd22->Draw();

//     auto leg2 = new TLegend();
//     drawAndAdd(hEffi_AS_NW_over_D);
//     drawAndAdd(hEffi_AS_WW_over_D);
//     //~ leg2->Draw("same");

//     saveCanvasAs(c1);
//     return c1;
// }

// TCanvas &
// compareMeasuredEffis_TH1_New(std::string const &fID,
//                              MCEffi_wRes &lMCEffi_AS,
//                              MCEffi_wRes &lMCEffi_D)
// {

//     auto &hEffi_AS_NW = lMCEffi_AS.SampleMeasuredEffi_NW_2(kPink);
//     auto &hEffi_AS_WW = *lMCEffi_AS.SampleMeasuredEffi_WW_2(kCyan);
//     auto &hEffi_D_NW = lMCEffi_D.SampleMeasuredEffi_NW_2(kMagenta);

//     auto &hEffi_AS_NW_over_D = *divideTH1ByTH1(hEffi_AS_NW, hEffi_D_NW, "", "hEffi_AS_NW_over_D_2");
//     auto &hEffi_AS_WW_over_D = *divideTH1ByTH1(hEffi_AS_WW, hEffi_D_NW, "", "hEffi_AS_WW_over_D_2");

//     auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
//     {
//         o.Draw("same");
//         if (leg)
//             leg->AddEntry(&o, o.GetName(), "l");
//     };
//     // setup canvas
//     TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TH1_New").data(),
//                                (fID + "compareMeasuredEffis_TH1_New").data(), 2000, 1000);

//     c1.Divide(1, 2);
//     c1.cd(1);
//     gPad->SetLogy();
//     TH2 *hd21 = new TH2F("hd21", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 8e-3);
//     hd21->Draw();

//     auto leg = new TLegend(.73, .64, .90, .90, "");
//     drawAndAdd(hEffi_AS_NW, leg);
//     drawAndAdd(hEffi_AS_WW, leg);
//     drawAndAdd(hEffi_D_NW, leg);
//     leg->Draw("same");

//     c1.cd(2);
//     //~ gPad->SetLogy();
//     TH2 *hd22 = new TH2F("hd22", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
//     hd22->Draw();

//     auto leg2 = new TLegend();
//     drawAndAdd(hEffi_AS_NW_over_D);
//     drawAndAdd(hEffi_AS_WW_over_D);
//     //~ leg2->Draw("same");

//     saveCanvasAs(c1);
//     return c1;
// }

/* check effect of normalizing factor
 * fHistoTruePrimaryMotherInvMassMCPt[fiCut]->Fill(
*       Pi0Candidate->M(),
        lMotherMCPt,
        weighted*fWeightJetJetMC*weightMatBudget);*/

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

    auto getGenDist = [&hGenDist_AS_inv, &lPtGaxis]()
    {
        TF1 *lGenDist_AS_dn_dptG_inv = new TF1("lGenDist_AS_dn_dptG_inv", "[0] + [1]/(x-[2])", lPtGaxis.GetXmin(), lPtGaxis.GetXmax());
        hGenDist_AS_inv->Fit(lGenDist_AS_dn_dptG_inv, "N");
        return multiplyTF1ByX(*lGenDist_AS_dn_dptG_inv, "lGenDist_AS_dn_dptG");
    };
    TF1 *lGenDist_AS_dn_dptG = getGenDist();

    // TAxis &axisPtR = lPtGaxis;
    TAxis lAxisPtR(100, 0., 10.);

    cout << "before MCEffi_wRes\n";
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
