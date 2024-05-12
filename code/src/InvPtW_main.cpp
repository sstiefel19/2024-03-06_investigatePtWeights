// purpose: find out how well pt weights work with resolution effect taken into account
#include "../include/InvPtW_main.h"

#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"

#include "../include/dN_dptR_integrand.h"

#include "/2024/2024-03-06_investigatePtWeights/code/include/computeResolutionFits.h"
#include "../include/MCEffi.h"
#include "../include/PtWeights.h"

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

TF1 &InvPtW_main::GetMesonEfficiency(std::string fname, Double_t theXmax /*= 10.*/)
{

    TH1 &hEffi = *(TH1 *)getObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dP/dptG", 1, 0., theXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();

    TF1 &fEffi = *new TF1("fEffi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., theXmax);
    fEffi.SetParameters(1, 1, -1, 0.);
    hEffi.Fit(&fEffi, "N");

    auto &leg = *new TLegend();
    drawAndAdd(hEffi, "same", kRed,
               &leg, "", "l", true /* theDrawLegAlready */);
    drawAndAdd(fEffi, "same", kBlue,
               &leg, "", "l");
    std::cout << "fEffi: " << fEffi.GetTitle() << std::endl;
    return fEffi;
}

// set first parameter to "auto" for auto-concatenated name
TF1 &InvPtW_main::FitGenDistHisto(std::string const &theResultNameInfo,
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
        return utils_TF1::GetDummyTF1(theTH1GenDist_dn_dptG.GetName(), true /*theAddDummyTag*/);
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

// detector parametrizations
utils_fits::TPairFitsWAxis &InvPtW_main::FitDetector(const std::string &theFnameInputEffiFit,
                                                     const std::string &theFnameResFits,
                                                     TH2 &theH2Resolution,
                                                     int thePtBinStart,
                                                     int thePtBinMax,
                                                     int theNRebin_r,
                                                     bool theDrawAllFitsOverlayed,
                                                     bool thePlotSingles)
{
    lEffiAtAll_dp_dptG = &GetMesonEfficiency(theFnameInputEffiFit);
    if (!lEffiAtAll_dp_dptG)
    {
        printf("investigatePtWeights_wResolutionEffects::FitDetector():\n\t"
               "ERROR: lEffiAtAll_dp_dptG is nullptr!\n"
               "\tReturning dummy TPairFitsWAxis  NOT.\n");
    }
    std::cout << "ger2\n";

    // 2) create resolution parametrizations
    ComputeResolutionFits &lCRF = *new ComputeResolutionFits(
        theH2Resolution,
        thePtBinStart,
        thePtBinMax,
        theNRebin_r,
        theFnameResFits,
        theDrawAllFitsOverlayed,
        thePlotSingles);
    std::cout << "ger2\n";
    std::cout << "ger4\n";
    
    utils_fits::TPairFitsWAxis &lResult = lCRF.Compute(theH2Resolution,
                                                       thePtBinStart,
                                                       thePtBinMax,
                                                       theNRebin_r,
                                                       theFnameResFits,
                                                       true /*drawAllFitsOverlayed*/,
                                                       false /*plotSingles*/);
    std::cout << "ger3\n";
    return lResult;
}

InvPtW_main::InvPtW_main(std::string const &theMeson,
                         std::string const &theCent,
                         std::string const &theCentAS,
                         std::string const &theFnameAS,
                         std::string const &theFnameInputEffiFit,
                         std::string const &theFnameWeightsFile,
                         int thePtBinStart,
                         int thePtBinMax,
                         bool theGenDistTH1IsInvariant,
                         bool theMultiplyResultTF1ByX)
    : meson(theMeson),
      cent(theCent),
      centAS(theCentAS),
      fnameAS(theFnameAS),
      lFnameInputEffiFit(theFnameInputEffiFit),
      lFnameResFits(Form("input_root/%s_resolutionFits_%d-%d.root",
                         centAS.data(),
                         ptBinStart,
                         ptBinMax)),
      fnameWeightsFile(theFnameWeightsFile),
      ptBinStart(thePtBinStart),
      ptBinMax(thePtBinMax),
      lGenDistTH1IsInvariant(theGenDistTH1IsInvariant),
      lMultiplyResultTF1ByX(theMultiplyResultTF1ByX),
      gAS({fnameAS, "GammaConvV1_997/", centAS, "_0d200009ab770c00amd0404000_0152101500000000"}),
      h2Resolution(*(TH2F *)gAS.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt")),
      lEffiAtAll_dp_dptG(nullptr)
{
    printf("invPtWeights_class::invPtWeights_class(): created instance for %s %s %s %s %s %s %s %d %d %d %d\n",
           meson.data(),
           cent.data(),
           centAS.data(),
           fnameAS.data(),
           lFnameInputEffiFit.data(),
           lFnameResFits.data(),
           fnameWeightsFile.data(),
           ptBinStart,
           ptBinMax,
           lGenDistTH1IsInvariant,
           lMultiplyResultTF1ByX);
}

PtWeights &InvPtW_main::SetupWeightsInstance(std::string const &theID,
                                             bool theComputedInInvariantForm,
                                             TAxis const &thePtGaxis)
{
    TH1 &hGenDist_AS_inv = *(TH1 *)getObjectFromPathInFile(
        fnameWeightsFile, meson + "_LHC24a1_5TeV_" + centAS);
    TF1 &fGenData_dn_dptG_inv = *(TF1 *)getObjectFromPathInFile(
        fnameWeightsFile, meson + "_Data_5TeV_" + cent.substr(0, 6));

    return *new PtWeights(theID,
                          theComputedInInvariantForm,
                          hGenDist_AS_inv,
                          fGenData_dn_dptG_inv,
                          thePtGaxis);
}

int InvPtW_main::Main()
{
    utils_fits::TPairFitsWAxis &lPair_vFits_ptG_i_dp_dr_Axis = FitDetector(lFnameInputEffiFit,
                                                                           lFnameResFits,
                                                                           h2Resolution,
                                                                           ptBinStart,
                                                                           ptBinMax,
                                                                           4 /*nR*/,
                                                                           true /*theDrawAllFitsOverlayed*/,
                                                                           false /*thePlotSingles*/);

    /*
        more accurate way to get the genDist:
        h_inv -> h -> f
        TH1 &hGenDist_AS_dn_dptG = *multiplyTH1ByBinCenters(*hGenDist_AS_inv, "", "hGenDist_AS_dn_dptG");
    */
    bool lGenDistTF1IsInvariant = false; // this initialized value is not used
    TAxis &lPtGaxis = lPair_vFits_ptG_i_dp_dr_Axis.second;
    PtWeights &lPtWeights = SetupWeightsInstance("lPtWeights",
                                                 lGenDistTF1IsInvariant,
                                                 lPtGaxis);

    // 4) fit the genDist
    TF1 &lGenDistTF1_dn_dptG_AS = FitGenDistHisto("auto",
                                                  *cloneTH1(lPtWeights.GetTH1MCGen_dn_dptG()), // not sure if clone is necessary
                                                  lGenDistTH1IsInvariant /*theTH1IsInvariant*/,
                                                  lPtGaxis /*thePtGaxis*/,
                                                  lMultiplyResultTF1ByX /*theMultiplyResultTF1ByX*/,
                                                  lGenDistTF1IsInvariant /*theResultIsInvariant_out*/);

    // 5 create MCEffi instances
    // TAxis &axisPtR = lPtGaxis;
    TAxis lAxisPtR(100, 0., 10.);
    auto &lMCEffi_AS = *new MCEffi_wRes("lMCEffi_AS",                 //
                                        lGenDistTF1_dn_dptG_AS,       // _fGenDist_dn_dptG
                                        *lEffiAtAll_dp_dptG,          // _fEffi_dp_dptG
                                        lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                        lAxisPtR,                     // _axisPtR
                                        &lPtWeights);

    TF1 &fGenData_dn_dptG = lPtWeights.GetTF1TrgtDist_dn_dptG();
    auto &lMCEffi_D = *new MCEffi_wRes("lMCEffi_D",                  //
                                       fGenData_dn_dptG,             // _fGenDist_dn_dptG
                                       *lEffiAtAll_dp_dptG,          // _fEffi_dp_dptG
                                       lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                       lAxisPtR);

    // 7) plot results
    lMCEffi_AS.SampleMeasuredEffi_NW_2();
    lMCEffi_AS.PlotAll();
    lMCEffi_D.PlotAll();

    return 0;
}