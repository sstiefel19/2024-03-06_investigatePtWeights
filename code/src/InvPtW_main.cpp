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

InvPtW_main::InvPtW_main(std::string const &theID,
                         std::string const &theFnameInputEffiFit,
                         std::string const &theFnameWeightsFile,
                         std::string const &theEvCutWeightsFile,
                         std::string const &theMeson,
                         GCo const &theGCo_h2Resolution,
                         int theNR)
    : sID(theID),
      sFnameFitOverallEfficiency(theFnameInputEffiFit),
      sFnameWeightsFile(theFnameWeightsFile),
      sEvCutWeightsFile(theEvCutWeightsFile),
      sMeson(theMeson),
      gGCo_h2Resolution(theGCo_h2Resolution),

      // derive
      // 2) create resolution parametrizations
      cRF(ComputeResolutionFits(
          gGCo_h2Resolution,
          1,  // thePtBinStart
          31, // thePtBinMax
          theNR /*theNR*/,
          true, // theDrawAllFitsOverlayed
          false /*thePlotSingles*/)),

      fTargetGenData_dn_dptG_inv(*(TF1 *)getObjectFromPathInFile(
          sFnameWeightsFile,
          sMeson + "_Data_5TeV_" + sEvCutWeightsFile.substr(0, 6))),

      hGenDist_dn_dptG_inv(*(TH1D *)getObjectFromPathInFile(
          sFnameWeightsFile,
          sMeson + "_LHC24a1_5TeV_" + sEvCutWeightsFile)),

      // compute
      hGenDist_dn_dptG(*(TH1D *)multiplyTH1ByBinCenters(hGenDist_dn_dptG_inv,
                                                        "",
                                                        "hGenDist_dn_dptG")),
      fTargetGenData_dn_dptG(utils_TF1::MultiplyTF1ByX(sID + "_fTargetGenData_dn_dptG",
                                                       fTargetGenData_dn_dptG_inv)),

      // output info
      fEffiAtAll_dp_dptG(nullptr),
      // helper structures
      aAxisPtG(nullptr) // gets intialized in ParametrizeEfficiencies
{
    printf("invPtWeights_class::invPtWeights_class(): created instance:\n"
           "\tsID: %s\n"
           "\tsFnameFitOverallEfficiency: %s\n"
           "\tsFnameWeightsFile: %s\n"
           "\tsEvCutWeightsFile: %s\n"
           "\tsMeson: %s\n"
           "\tgGCo_h2Resolution: %s\n"
           "\tcRF: %s\n"
           "\tfTargetGenData_dn_dptG_inv: %s\n"
           "\thGenDist_dn_dptG_inv: %s\n"
           "\thGenDist_dn_dptG: %s\n"
           "\tfTargetGenData_dn_dptG: %s\n"
           "\tfEffiAtAll_dp_dptG: %s\n"
           "\taAxisPtG: %s\n",
           sID.data(),
           sFnameFitOverallEfficiency.data(),
           sFnameWeightsFile.data(),
           sEvCutWeightsFile.data(),
           sMeson.data(),
           gGCo_h2Resolution.fname.data(),
           cRF.GetID().data(),
           fTargetGenData_dn_dptG_inv.GetName(),
           hGenDist_dn_dptG_inv.GetName(),
           hGenDist_dn_dptG.GetName(),
           fTargetGenData_dn_dptG.GetName(),
           "nullptr",
           "nullptr");
}

// ===================== public member functions =================================
int InvPtW_main::Main(bool theCalcPtWInInvForm)
{
    // 1) fit overall efficiency
    // this call also initialized aAxisPtG
    int lNRebin_r = 4;
    utils_fits::TPairFitsWAxis &lPair_vFits_ptG_i_dp_dr_Axis =
        ParametrizeEfficiencies();

    // 2) fit the genDist
    TF1 &lGenDistTF1_dn_dptG_AS =
        FitMCGeneratedParticlesHisto("auto",
                                     theCalcPtWInInvForm ? hGenDist_dn_dptG_inv /*theTH1GenDist_dn_dptG*/
                                                         : hGenDist_dn_dptG,
                                     theCalcPtWInInvForm /*theTH1IsInvariant*/,
                                     theCalcPtWInInvForm /*theMultiplyResultTF1ByX*/);

    // 3) create PtWeights instance
    PtWeights &lPtWeights = CreatePtWeightsInstance(sID + "_lPtWeights",
                                                    theCalcPtWInInvForm);

    // 3 create MCEffi instances
    TAxis lAxisPtR(100, 0., 10.);
    auto &lMCEffi_AS = *new MCEffi(sID + "_lMCEffi_AS_" + lPtWeights.GetID(), //
                                   lGenDistTF1_dn_dptG_AS,                    // _fGenDist_dn_dptG
                                   *fEffiAtAll_dp_dptG,                       // _fEffi_dp_dptG
                                   lPair_vFits_ptG_i_dp_dr_Axis,              // _vFits_ptG_i_dp_dr_wAxis
                                   lAxisPtR,                                  // _axisPtR
                                   &lPtWeights);

    auto &lMCEffi_D = *new MCEffi(sID + "_lMCEffi_D",           //
                                  fTargetGenData_dn_dptG,       // _fGenDist_dn_dptG
                                  *fEffiAtAll_dp_dptG,          // _fEffi_dp_dptG
                                  lPair_vFits_ptG_i_dp_dr_Axis, // _vFits_ptG_i_dp_dr_wAxis
                                  lAxisPtR);

    // 4) plot results
    lMCEffi_AS.SampleMeasuredEffi_NW_2();
    lMCEffi_AS.PlotAll();
    lMCEffi_D.PlotAll();

    return 0;
}

// ===================== private member functions =========================
PtWeights &InvPtW_main::CreatePtWeightsInstance(std::string const &theID,
                                                bool theComputeInInvariantForm)
{
    std::string lFullName(theID + theComputeInInvariantForm ? "_inv" : "_special");
    if (!aAxisPtG)
    {
        printf("investigatePtWeights_wResolutionEffects::createPtWeightsInstance():\n\t"
               "ERROR: aAxisPtG is nullptr!\n"
               "\tReturning dummy PtWeights instance.\n");
        return *new PtWeights(lFullName);
    }

    return *new PtWeights(
        lFullName,
        theComputeInInvariantForm,
        theComputeInInvariantForm ? hGenDist_dn_dptG_inv
                                  : hGenDist_dn_dptG,
        theComputeInInvariantForm ? fTargetGenData_dn_dptG
                                  : fTargetGenData_dn_dptG_inv,
        *aAxisPtG);
}

// detector parametrizations
utils_fits::TPairFitsWAxis &InvPtW_main::ParametrizeEfficiencies()
{
    fEffiAtAll_dp_dptG = &GetMesonEfficiency(sFnameFitOverallEfficiency);
    if (!fEffiAtAll_dp_dptG)
    {
        printf("investigatePtWeights_wResolutionEffects::ParametrizeEfficiencies():\n\t"
               "ERROR: lEffiAtAll_dp_dptG is nullptr!\n"
               "\tReturning dummy TPairFitsWAxis  NOT.\n");
    }

    utils_fits::TPairFitsWAxis &lResult = cRF.Compute();
    aAxisPtG = &lResult.second;
    return lResult;
}

// set first parameter to "auto" for auto-concatenated name
TF1 &InvPtW_main::FitMCGeneratedParticlesHisto(std::string const &theResultNameInfo,
                                               TH1 &theTH1GenDist_dn_dptG,
                                               bool theTH1IsInvariant,
                                               bool theMultiplyResultTF1ByX)
{
    // on top for readability
    bool lResultIsInvariant_out = theTH1IsInvariant && !theMultiplyResultTF1ByX;

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
                                           lResultIsInvariant_out ? "_inv" : ""));

    // create result TF1
    TF1 &lTF1GenDist_dn_dptG = *new TF1(lResultFullName.data(),
                                        theTH1IsInvariant ? "[0] + [1]/(x-[2])"
                                                          : "[0]",
                                        aAxisPtG->GetXmin(),
                                        aAxisPtG->GetXmax());
    // perform fit
    theTH1GenDist_dn_dptG.Fit(&lTF1GenDist_dn_dptG, "N");

    return (!theMultiplyResultTF1ByX)
               ? lTF1GenDist_dn_dptG
               : utils_TF1::MultiplyTF1ByX(lResultFullName,
                                           lTF1GenDist_dn_dptG);
}

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
