// purpose: find out how well pt weights work with resolution effect taken into account
#include "/analysisSoftware/SupportingMacros/GCo.h"
#include "/analysisSoftware/SupportingMacros/utils_files_strings.h"
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

    TH1 &hEffi = *(TH1 *)utils_files_strings::GetObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dP/dptG", 1, 0., theXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();

    TF1 &fEffi = *new TF1("fEffi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., theXmax);
    fEffi.SetParameters(1, 1, -1, 0.);
    hEffi.Fit(&fEffi, "N");

    auto &leg = *new TLegend();
    utils_plotting::DrawAndAdd(hEffi, "same", kRed,
                               &leg, "", "l", true /* theDrawLegAlready */);
    utils_plotting::DrawAndAdd(fEffi, "same", kBlue,
                               &leg, "", "l");
    return fEffi;
}

// set first parameter to "auto" for auto-concatenated name
TF1 &fitGenDistHisto(std::string const &theResultNameInfo,
                     TH1 &theTH1GenDist_dn_dptG,
                     bool theTH1IsInvariant,
                     TAxis const &thePtGaxis,
                     bool theMultiplyResultTF1ByX,
                     bool &theResultIsInvariant_out)
{
    // on top for readability
    theResultIsInvariant_out = theTH1IsInvariant && !theMultiplyResultTF1ByX;

    // must not be empty
    if (!theResultNameInfo.size())
    {
        printf("investigatePtWeights_wResolutionEffects::createGenDistFit():\n\t"
               "ERROR: theResultNameInfo can't be empty!\n"
               "\tChose explicit name or \'auto\' for auto-concatenated name (verbose).\n"
               "\tReturning dummy TF1.\n");
        return utils_TF1::CreateDummyTF1(theTH1GenDist_dn_dptG.GetName(), true /*theAddDummyTag*/);
    }

    // compile full result name
    std::string lResultFullName(
        theResultNameInfo != "auto" ? theResultNameInfo
                                    : Form("createGenDistFit_dn_dptG_from_%s_%s",
                                           theTH1GenDist_dn_dptG.GetName(),
                                           theResultIsInvariant_out ? "_inv" : ""));

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
    // set some globals
    {
        gROOT->Reset();
        //~ gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);
        gStyle->SetLegendTextSize(.03);
    }

    // todo: put these definitions into a header file
    // 0) definitions
    int ptBinStart = 1;
    int ptBinMax = 31;

    bool lGenDistTH1IsInvariant = true;
    bool lMultiplyResultTF1ByX = true;

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

    // ==================== end 0) defintions ==========================

    auto setupWeightsInstance =
        [&meson, &cent, &centAS, &fnameWeightsFile](std::string const &theID,
                                                    bool theComputedInInvariantForm,
                                                    TAxis const &thePtGaxis)
    {
        // get Infos for Weights instance
        TH1 &hGenDist_AS_inv = *(TH1 *)utils_files_strings::GetObjectFromPathInFile(
            fnameWeightsFile, meson + "_LHC24a1_5TeV_" + centAS);
        TF1 &fGenData_dn_dptG_inv = *(TF1 *)utils_files_strings::GetObjectFromPathInFile(
            fnameWeightsFile, meson + "_Data_5TeV_" + cent.substr(0, 6));

        return new PtWeights(theID,
                             theComputedInInvariantForm,
                             hGenDist_AS_inv,
                             fGenData_dn_dptG_inv,
                             thePtGaxis);
    };

    // // detector parametrizations
    // {

    // 1) get some parametrized effi
    TF1 &fMCIntrinsicEffiAtAll_dp_dptG = getMesonEfficiency(lFnameInputEffiFit);

    // 2) create resolution parametrizations
    utils_fits::TPairFitsWAxis lPair_vFits_ptG_i_dp_dr_Axis =
        computeResolutionFits(h2Resolution,
                              ptBinStart,
                              ptBinMax,
                              4 /*nRebin_r*/,
                              lFnameResFits,
                              true /*drawAllFitsOverlayed*/,
                              true /*plotSingles*/);
    TAxis &lPtGaxis = lPair_vFits_ptG_i_dp_dr_Axis.second;
    //}

    /*
        more accurate way to get the genDist:
        h_inv -> h -> f
        TH1 &hGenDist_AS_dn_dptG = *utils_TH1::MultiplyTH1ByBinCenters(*hGenDist_AS_inv, "", "hGenDist_AS_dn_dptG");
    */
    bool lGenDistTF1IsInvariant = false; // this initialized value is not used
    PtWeights &lPtWeights = *setupWeightsInstance("lPtWeights",
                                                  lGenDistTF1IsInvariant,
                                                  lPtGaxis);

    // 3 obtain the genDist that was fed to lPtWeights
    // make a copy since we need to fit the histogram which stupidly changes the histogram itself
    TH1 &hGenDist_AS = *utils_utils::CloneTH1(lPtWeights.GetTH1MCGen_dn_dptG());

    // 4) fit the genDist
    TF1 &lGenDistTF1_dn_dptG_AS = fitGenDistHisto("auto",
                                                  hGenDist_AS,
                                                  lGenDistTH1IsInvariant /*theTH1IsInvariant*/,
                                                  lPtGaxis,
                                                  lMultiplyResultTF1ByX /*theMultiplyResultTF1ByX*/,
                                                  lGenDistTF1IsInvariant /*theResultIsInvariant_out*/);

    // 5 create MCEffi instances
    // TAxis &axisPtR = lPtGaxis;
    TAxis lAxisPtR(100, 0., 10.);
    auto &lMCEffi_AS = *new MCEffi("lMCEffi_AS",                  //
                                   lGenDistTF1_dn_dptG_AS,        // _fGenDist_dn_dptG
                                   fMCIntrinsicEffiAtAll_dp_dptG, // _fEffi_dp_dptG
                                   lPair_vFits_ptG_i_dp_dr_Axis,  // _vFits_ptG_i_dp_dr_wAxis
                                   lAxisPtR,                      // _axisPtR
                                   &lPtWeights);

    TF1 &fGenData_dn_dptG = lPtWeights.GetTF1TrgtDist_dn_dptG();
    auto &lMCEffi_D = *new MCEffi("lMCEffi_D",                   //
                                  fGenData_dn_dptG,              // _fGenDist_dn_dptG
                                  fMCIntrinsicEffiAtAll_dp_dptG, // _fEffi_dp_dptG
                                  lPair_vFits_ptG_i_dp_dr_Axis,  // _vFits_ptG_i_dp_dr_wAxis
                                  lAxisPtR);

    // 7) plot results
    lMCEffi_AS.SampleMeasuredEffi_MC_NW_2();
    lMCEffi_AS.PlotAll();
    lMCEffi_D.PlotAll();
}
