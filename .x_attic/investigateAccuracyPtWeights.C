/*purpose: find out how big the inaccuracy is that we introduce by using
 * histo->Eval(pt) in the conversion task when calculating the weights.
 * result: blue markers in bottom panel. They show how different I will
 *         I will meausure the efficicieny in both MCs even though they
 *         actually have the same.
 * This completely done in MCpt, so no resolution effect taken into account.*/

#include "MCEffi.h"

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

Double_t lXmax = 10.;

TF1 *getMesonEfficiency(std::string fname)
{

    TH1 *hEffi = (TH1 *)utils_files_strings::GetObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dN/dpT", 1, 0., lXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();

    auto leg2 = new TLegend();

    auto utils_plotting::DrawAndAddToLegendF = [leg2](TF1 *f, EColor c = kBlue)
    {f->SetLineColor(c); f->Draw("same"); leg2->AddEntry(f,f->GetName()); };

    hEffi->Draw("same");
    leg2->Draw("same");

    TF1 *fEffi = new TF1("fEffi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., lXmax);
    fEffi->SetParameters(1, 1, -1, 0.);

    hEffi->Fit(fEffi, "N");

    MyDerivative *fptr = new MyDerivative(fEffi, 1); // create the user function class
    auto fEffi_d1 = fptr->GetD();

    MyDerivative *fptr2 = new MyDerivative(fEffi, 2); // create the user function class
    auto fEffi_d2 = fptr2->GetD();

    utils_plotting::DrawAndAddToLegendF(fEffi, kRed);
    utils_plotting::DrawAndAddToLegendF(fEffi_d1);
    utils_plotting::DrawAndAddToLegendF(fEffi_d2, kGreen);

    return fEffi;
}

// change ratio for weights calculation to invariant functions
void investigateAccuracyPtWeights()
{

    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);

    Double_t lYmax = 3000;

    std::string meson("Pi0");
    std::string cent("101");

    // get some parametrized effi
    TF1 *myProviEffi = getMesonEfficiency("out_101.root");

    std::string fname(Form("/2023/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root"));
    std::string functionname(meson + "_Data_5TeV_" + cent + "300");
    std::string histoname(meson + "_LHC20e3a_5TeV_" + cent + "30053");
    std::string histonameA(meson + "_LHC24a1_5TeV_" + cent + "30023");

    TFile infile(fname.data(), "READ");
    TF1 *fD_inv = (TF1 *)infile.Get(functionname.data())->Clone((functionname + "_clone").data());
    TH1 *hMC_inv = (TH1 *)utils_files_strings::GetObjectFromPathInFile(infile, histoname, histoname + "_clone");
    TH1 *hMC_invA = (TH1 *)utils_files_strings::GetObjectFromPathInFile(infile, histonameA, histoname + "_cloneA");
    infile.Close();

    // setup canvas
    TCanvas *c1 = new TCanvas("c1", "c1", 2000, 1000);
    c1->Divide(1, 2);
    c1->cd(1);
    gPad->SetLogy();
    TH2 *hd21 = new TH2F("hd21", ";MC pT(GeV);dN/dpT", 1, 0., lXmax, 1, 1e-6, lYmax);
    hd21->Draw();
    TLegend *leg = new TLegend();
    leg->Draw("same");

    MCEffi lMCMB("MB_", hMC_inv, fD_inv, myProviEffi);
    MCEffi lMCAS("AS_", hMC_invA, fD_inv, myProviEffi);

    lMCMB.PlotAll(leg);

    c1->cd(2);
    TH2 *hd22 = new TH2F("hd22", ";MC pT(GeV);ratios", 1, 0., lXmax, 1, 0.8, 1.2);
    hd22->Draw();
    TLegend *leg2 = new TLegend();
    leg2->Draw("same");

    auto utils_plotting::DrawAndAddToLegend = [leg2](TH1 *h)
    {reduceErrorsByFactor(*h); h->Draw("same"); leg2->AddEntry(h,h->GetName()); };
    auto utils_plotting::DrawAndAddToLegendF = [leg2](TF1 *f)
    {f->Draw("same"); leg2->AddEntry(f,f->GetName()); };

    //~ TF1* fEffi_d2 = new TF1("fEffi_d2",[&](double*x, double *p){ return myProviEffi->Derivative2(*x); }, 0, lXmax, 0);
    //~ utils_plotting::DrawAndAddToLegendF(fEffi_d2);

    TH1 *hEffiWW_ASoverMB = utils_TH1::DivideTH1ByTH1(*lMCAS.GetMeasuredEffiWW(), *lMCMB.GetMeasuredEffiWW(), "", "hEffiWW_ASoverMB");
    utils_plotting::DrawAndAddToLegend(hEffiWW_ASoverMB);

    TH1 *hEffiNW_ASoverMB = utils_TH1::DivideTH1ByTH1(*lMCAS.GetMeasuredEffiWoW(), *lMCMB.GetMeasuredEffiWoW(), "", "hEffiNW_ASoverMB");

    TH1 *hASGenWWoverMBGenWW = utils_TH1::DivideTH1ByTH1(*lMCAS.hMCGenSampledWW, *lMCMB.hMCGenSampledWW, "", "hASGenWWoverMBGenWW");
    hASGenWWoverMBGenWW->SetLineColor(kRed);
    utils_plotting::DrawAndAddToLegend(hASGenWWoverMBGenWW);

    //~ TH1* hEffiWWoverNW = utils_TH1::DivideTH1ByTH1(*lMCMB.GetMeasuredEffiWW(), *lMCMB.GetMeasuredEffiWoW());

    c1->SaveAs(Form("investigateAccuracyPtWeights_%s_%s.png", meson.data(), cent.data()));
}
