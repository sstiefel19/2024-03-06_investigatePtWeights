#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2024.h"

#include <iostream>
#include <string.h>
#include <vector>

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TSystem.h"

/* currently:
    0) get h_inv from weights file
    0.5) h_inv -> h -> f
    1) assume function f for mc generated particles
    2) sample a histo h from f
    3) transform it into inv. yield: h_inv = 1/pT h
    4) use TH1::Interpolate at many points along h_inv
    5) compare to the same points on 1/f 

    idea:
    0) - 2) same as above
    calculate weights in the "normal" distributions, not the invariant one 
    */ 

void testTH1continuosEvaluationMethods() {
        
    std::string const lNameMacro("testTH1continuosEvaluationMethods");

    std::string centAS("10130023");
    std::string meson("Pi0");
    
    std::string fnameWeightsFile(Form(
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root"));

    // prepare canvas already here
    TCanvas &c1 = getCanvasWithTH2F(lNameMacro, 
                                    ";pT (GeV);[(1/pT)] dN/dpT", 
                                    0., 10., 
                                    0., 50., 
                                    false/*theLogY*/,
                                    2000,
                                    1000,
                                    nullptr, nullptr,
                                    .02 /*theNameTagTextSize*/,
                                    kGray, /*theNameTagTextColor*/
                                    1, 2, 1);
    auto leg = getLegend();
    
    // 0) get h_inv from weights file
    TH1 &hGenDist_AS_inv = *(TH1 *)getObjectFromPathInFile(
        fnameWeightsFile, meson + "_LHC24a1_5TeV_" + centAS);
    
    // h_inv -> h
    TH1 &hGenDist_AS_dn_dptG = *multiplyTH1ByBinCenters(hGenDist_AS_inv, "", "hGenDist_AS_dn_dptG");
    
    // h -> f
    TF1 &fGenDist_dn_dptG = *new TF1("fGenDist_dn_dptG", "[0]", 0., 10.);
    hGenDist_AS_dn_dptG.Fit(&fGenDist_dn_dptG, "N");

    // 2) sample a histo h from f
    TH1 &hGenDist_AS_dn_dptG_sampled = *getSampledH(fGenDist_dn_dptG, hGenDist_AS_inv);

    // 3) transform it into inv. yield: h_inv = 1/pT h
    TH1 &hGenDist_sampled_inv = *divideTH1ByBinCenters(hGenDist_AS_dn_dptG_sampled, "", "hGenDist_sampled_inv");

    drawAndAdd(hGenDist_AS_inv, "same l", kBlack, leg);
    drawAndAdd(hGenDist_AS_dn_dptG, "same l", kBlue, leg);
    drawAndAdd(fGenDist_dn_dptG, "same l", kRed, leg);
    drawAndAdd(hGenDist_AS_dn_dptG_sampled, "same l", kGreen, leg);
    drawAndAdd(hGenDist_sampled_inv, "same l", kOrange, leg);

    // 4) use TH1::Interpolate at many points along h_inv
    // use a fine binned histogram to move in small discrete steps along the x-axis
    TH1 &hExistingEval = *new TH1F("hExistingEval", ";pT (GeV);1/pT dN/dpT", 1000, 0., 10.);
    TH1 &hF_inv = *cloneTH1(hExistingEval, "", "hF_inv", "hF_inv");

    for (int i = 1; i <= hExistingEval.GetNbinsX(); i++) {
        float x = hExistingEval.GetBinCenter(i);
        hExistingEval.SetBinContent(i, hGenDist_sampled_inv.Interpolate(x));
        hF_inv.SetBinContent(i, x ? fGenDist_dn_dptG.Eval(x) / x : 0.);
    }
    drawAndAdd(hExistingEval, "same l", kMagenta, leg);
    drawAndAdd(hF_inv, "same l", kCyan, leg);

    // 5) compare to the same points on f
    c1.cd(2);
    TH2 &hd21 = *new TH2F("hd1", ";pTG (GeV);ratio", 1, 0., 10., 1, .8, 1.2);
    hd21.Draw();
    TH1 &hRatio = *divideTH1ByTH1(hExistingEval, hF_inv, "", "hRatio");
    hRatio.Draw("same");
    gPad->SetGridy();
    saveCanvasAs(c1, "pdf");
}
