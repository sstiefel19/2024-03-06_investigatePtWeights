// purpose: find out how well pt weights work with resolution effect taken into account
#include "../include/InvPtW.h"

#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_plotting.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TH1.h"

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

InvPtW::InvPtW(std::string const &theID,
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
      tGCo_forH2Resolution(theGCo_h2Resolution),

      // derive
      // 2) create resolution parametrizations
      vAllDrawableObjects(),
      tResolutionFits(ComputeResolutionFits(
          tGCo_forH2Resolution,
          1,  // thePtBinStart
          31, // thePtBinMax
          theNR /*theNR*/,
          true, // theDrawAllFitsOverlayed
          false /*thePlotSingles*/,
          &vAllDrawableObjects)),

      fTargetGenData_dn_dptG_inv(*(TF1 *)utils_files_strings::GetObjectFromPathInFile(
          sFnameWeightsFile,
          sMeson + "_Data_5TeV_" + sEvCutWeightsFile.substr(0, 6))),

      hGenDist_dn_dptG_inv(*(TH1D *)utils_files_strings::GetObjectFromPathInFile(
          sFnameWeightsFile,
          sMeson + "_LHC24a1_5TeV_" + sEvCutWeightsFile)),

      // compute
      hGenDist_dn_dptG(*(TH1D *)utils_TH1::MultiplyTH1ByBinCenters(hGenDist_dn_dptG_inv,
                                                                   "",
                                                                   "hGenDist_dn_dptG", "hGenDist_dn_dptG")),

      fTargetGenData_dn_dptG(utils_TF1::MultiplyTF1ByX(sID + "_fTargetGenData_dn_dptG",
                                                       fTargetGenData_dn_dptG_inv)),

      // get initialized in Initialize
      fEffiAtAll_dp_dptG(nullptr),
      tPair_vFits_ptG_i_dp_dr_Axis(nullptr),
      aAxisPtG(nullptr),
      fGenDistTF1_dn_dptG_AS(nullptr),
      fGenDistTF1_dn_dptG_AS_inv(nullptr),
      tMCEffi_D(nullptr),
      tMCEffi_AS_inv(nullptr),
      tMCEffi_AS_special(nullptr),
      bInitialized(false),
      vAllMCEffis()
{
    vAllDrawableObjects.insert(vAllDrawableObjects.end(),
                               {&fTargetGenData_dn_dptG_inv,
                                &fTargetGenData_dn_dptG,
                                &hGenDist_dn_dptG_inv,
                                &hGenDist_dn_dptG});

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
           "\tfTargetGenData_dn_dptG: %s\n",
           sID.data(),
           sFnameFitOverallEfficiency.data(),
           sFnameWeightsFile.data(),
           sEvCutWeightsFile.data(),
           sMeson.data(),
           tGCo_forH2Resolution.fname.data(),
           tResolutionFits.GetID().data(),
           fTargetGenData_dn_dptG_inv.GetName(),
           hGenDist_dn_dptG_inv.GetName(),
           hGenDist_dn_dptG.GetName(),
           fTargetGenData_dn_dptG.GetName());
}

// ===================== public member functions =================================
int InvPtW::Initialize()
{
    if (bInitialized)
    {
        printf("int InvPtW_main::Initialize(): instance %s: already initialized. returning.\n", sID.data());
    }
    printf("\n\n\nint InvPtW_main::Initialize(): instance %s: starting to initialize.\n", sID.data());

    // 1) fit overall efficiency
    fEffiAtAll_dp_dptG = &GetMesonEfficiency(sFnameFitOverallEfficiency);

    // 2) create resolution parametrizations
    tPair_vFits_ptG_i_dp_dr_Axis = &tResolutionFits.Compute();
    aAxisPtG = &tPair_vFits_ptG_i_dp_dr_Axis->second;

    // 3) fit the genDist, once for inv once for not inv
    fGenDistTF1_dn_dptG_AS =
        &FitMCGeneratedParticlesHisto("auto",
                                      hGenDist_dn_dptG, /*theTH1GenDist_dn_dptG_x*/
                                      false /*theTH1IsInvariant*/);

    fGenDistTF1_dn_dptG_AS_inv =
        &FitMCGeneratedParticlesHisto("auto",
                                      hGenDist_dn_dptG_inv /*theTH1GenDist_dn_dptG_x*/,
                                      true /*theTH1IsInvariant*/);

    // 4) create MCEffi instances
    tMCEffi_D = new MCEffi(sID + "_lMCEffi_D",
                           fTargetGenData_dn_dptG,
                           *fEffiAtAll_dp_dptG,
                           *tPair_vFits_ptG_i_dp_dr_Axis,
                           *aAxisPtG,
                           nullptr,
                           &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_D);

    tMCEffi_AS_inv = new MCEffi(sID + "_lMCEffi_AS",
                                *fGenDistTF1_dn_dptG_AS_inv,
                                *fEffiAtAll_dp_dptG,
                                *tPair_vFits_ptG_i_dp_dr_Axis,
                                *aAxisPtG,
                                &CreatePtWeightsInstance(sID + "_lPtWeights",
                                                         true),
                                &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_AS_inv);

    tMCEffi_AS_special = new MCEffi(sID + "_lMCEffi_AS_special",
                                    *fGenDistTF1_dn_dptG_AS,
                                    *fEffiAtAll_dp_dptG,
                                    *tPair_vFits_ptG_i_dp_dr_Axis,
                                    *aAxisPtG,
                                    &CreatePtWeightsInstance(sID + "_lPtWeights",
                                                             false),
                                    &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_AS_special);

    printf("int InvPtW_main::Initialize(): instance %s: done initializing.\n\n\n", sID.data());
    bInitialized = true;
    return 1;
}

TCanvas &InvPtW::CompareAllMeasuredEfficiencies(TLegend *theLeg /*= nullptr*/)
{
    std::string lMethodName(sID + "_CompareAllMeasuredEfficiencies");
    TCanvas &lCanvas = utils_plotting::GetCanvasWithTH2F(lMethodName + "_canvas",
                                                         lMethodName /*theTitleH*/,
                                                         0., 10.5,
                                                         1.e-6, 1.e-2,
                                                         true /*theLogY*/,
                                                         2000, 1000,
                                                         &lMethodName /*theTitleC*/,
                                                         &lMethodName,
                                                         0.03 /*theNameTagTextSize*/,
                                                         kGray /*theNameTagTextColor*/,
                                                         1, 2, /*theNx, theNy*/
                                                         1 /*theDrawFirstOnI*/);
    TLegend *lLeg = theLeg ? theLeg : new TLegend(.73, .64, .90, .90, "");
    utils_plotting::DrawAndAdd(GetMCEffi_D().GetMeasuredEffi_NW_clone(), "same", kGreen,
                               lLeg, "", "l", true /* theDrawLegAlready */);
    utils_plotting::DrawAndAdd(GetMCEffi_AS_inv().GetMeasuredEffi_NW_clone(), "same", kBlue,
                               lLeg, "", "l");
    utils_plotting::DrawAndAdd(*GetMCEffi_AS_inv().GetMeasuredEffi_WW_clone(), "same", kRed,
                               lLeg, "", "l");
    utils_plotting::DrawAndAdd(*GetMCEffi_AS_special().GetMeasuredEffi_WW_clone(), "same", kMagenta,
                               lLeg, "", "l");

    // utils_plotting::SaveCanvasAs(lCanvas);
    return lCanvas;
}

int InvPtW::Main()
{
    if (!bInitialized)
    {
        printf("int InvPtW_main::Main(): Instance%s::Main() called for the first time. Initializing now.\n", sID.data());
        Initialize();
    }

    PlotAll();

    return 0;
}

void InvPtW::PlotAll()
{
    if (!bInitialized)
    {
        printf("int InvPtW_main::PlotAll(): instance %s: not initialized. Returning.\n", sID.data());
        return;
    }

    for (auto &lMCEffi : vAllMCEffis)
    {
        lMCEffi->PlotAll();
    }
}

// ===================== private member functions =========================
PtWeights &InvPtW::CreatePtWeightsInstance(std::string const &theID,
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

// set first parameter to "auto" for auto-concatenated name
TF1 &InvPtW::FitMCGeneratedParticlesHisto(std::string const &theResultNameInfo,
                                          TH1 &theTH1GenDist_dn_dptG_x, // x = inv or not
                                          bool theTH1IsInvariant)
{
    // must not be empty
    if (!theResultNameInfo.size())
    {
        printf("investigatePtWeights_wResolutionEffects::createGenDistFit():\n\t"
               "ERROR: theResultNameInfo can't be empty!\n"
               "\tChose explicit name or \'auto\' for auto-concatenated name (verbose).\n"
               "\tReturning dummy TF1.\n");
        return utils_TF1::CreateDummyTF1(theTH1GenDist_dn_dptG_x.GetName(), true /*theAddDummyTag*/);
    }

    // compile full result name
    std::string lResultFullName(
        theResultNameInfo != "auto" ? theResultNameInfo
                                    : Form("createGenDistFit_dn_dptG_from_%s_%s",
                                           theTH1GenDist_dn_dptG_x.GetName(),
                                           theTH1IsInvariant ? "_inv" : ""));

    // create result TF1
    TF1 &lTF1GenDist_dn_dptG_x = *new TF1(lResultFullName.data(),
                                          theTH1IsInvariant ? "[0] + [1]/(x-[2])"
                                                            : "[0]",
                                          aAxisPtG->GetXmin(),
                                          aAxisPtG->GetXmax());
    // perform fit
    theTH1GenDist_dn_dptG_x.Fit(&lTF1GenDist_dn_dptG_x, "N");
    return lTF1GenDist_dn_dptG_x;
}

TF1 &InvPtW::GetMesonEfficiency(std::string fname, Double_t theXmax /*= 10.*/)
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
    std::cout << "fEffi: " << fEffi.GetTitle() << std::endl;
    return fEffi;
}
