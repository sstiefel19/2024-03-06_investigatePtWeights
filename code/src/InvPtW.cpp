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

      hMCGenDist_dn_dptG_inv(*(TH1D *)utils_files_strings::GetObjectFromPathInFile(
          sFnameWeightsFile,
          sMeson + "_LHC24a1_5TeV_" + sEvCutWeightsFile)),

      // compute
      hMCGenDist_dn_dptG(*(TH1D *)utils_TH1::MultiplyTH1ByBinCenters(hMCGenDist_dn_dptG_inv,
                                                                     "",
                                                                     "hMCGenDist_dn_dptG", "hMCGenDist_dn_dptG")),

      fTargetGenData_dn_dptG(utils_TF1::MultiplyTF1ByX(sID + "_fTargetGenData_dn_dptG",
                                                       fTargetGenData_dn_dptG_inv)),

      // get initialized in Initialize
      fMCIntrinsicEffiAtAll_dp_dptG(nullptr),
      tPair_vFits_ptG_i_dp_dr_Axis(nullptr),
      aAxisPtG(nullptr),
      fMCGenDistTF1_dn_dptG_AS(nullptr),
      fMCGenDistTF1_dn_dptG_AS_inv(nullptr),
      tMCEffi_D(nullptr),
      tMCEffi_AS_inv_ptW(nullptr),
      tMCEffi_AS_special_ptW(nullptr),
      bInitialized(false),
      vAllMCEffis()
{
    vAllDrawableObjects.insert(vAllDrawableObjects.end(),
                               {&fTargetGenData_dn_dptG_inv,
                                &fTargetGenData_dn_dptG,
                                &hMCGenDist_dn_dptG_inv,
                                &hMCGenDist_dn_dptG});

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
           hMCGenDist_dn_dptG_inv.GetName(),
           hMCGenDist_dn_dptG.GetName(),
           fTargetGenData_dn_dptG.GetName());
}

//  ===================== public member functions =================================
TCanvas &InvPtW::CompareGeneratedSpectra(TLegend &theLeg)
{
    std::string const &lObjVarNameBase("fGen_dn_dptG");
    return CompareObservables_generic(lObjVarNameBase,
                                      "all" /*theSelectWhich*/,
                                      Form("%s;ptG (GeV);dN/dptG (1./GeV)",
                                           lObjVarNameBase.data()),
                                      theLeg);
}

TCanvas &InvPtW::CompareMeasuredEfficiencies(TLegend &theLeg)
{
    std::string const &lObjVarNameBase("MeasuredEffi");
    return CompareObservables_generic(lObjVarNameBase,
                                      "all", /*theSelectWhich*/
                                      Form("%s;ptR (GeV);efficiency",
                                           lObjVarNameBase.data()),
                                      theLeg);
}

int InvPtW::Initialize()
{
    if (bInitialized)
    {
        printf("int InvPtW_main::Initialize(): INFO: instance %s: already initialized. returning.\n", sID.data());
    }
    printf("\n\n\nint InvPtW_main::Initialize(): INFO: instance %s: starting to initialize.\n", sID.data());

    // 1) fit overall efficiency
    fMCIntrinsicEffiAtAll_dp_dptG = &FitTrueDetectorIntrinsicMesonEfficiency(sFnameFitOverallEfficiency);

    // 2) create resolution parametrizations
    tPair_vFits_ptG_i_dp_dr_Axis = &tResolutionFits.Compute();
    aAxisPtG = &tPair_vFits_ptG_i_dp_dr_Axis->second;

    // 3) fit the genDist, once for inv once for not inv
    fMCGenDistTF1_dn_dptG_AS =
        &FitMCGeneratedParticlesHisto("auto",
                                      hMCGenDist_dn_dptG, /*theTH1GenDist_dn_dptG_x*/
                                      false /*theTH1IsInvariant*/);

    fMCGenDistTF1_dn_dptG_AS_inv =
        &FitMCGeneratedParticlesHisto("auto",
                                      hMCGenDist_dn_dptG_inv /*theTH1GenDist_dn_dptG_x*/,
                                      true /*theTH1IsInvariant*/);

    // 4) create MCEffi instances
    tMCEffi_D = new MCEffi(sID + "_lMCEffi_D",
                           fTargetGenData_dn_dptG,
                           *fMCIntrinsicEffiAtAll_dp_dptG,
                           *tPair_vFits_ptG_i_dp_dr_Axis,
                           *aAxisPtG,
                           nullptr,
                           &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_D);

    tMCEffi_AS_inv_ptW = new MCEffi(sID + "_lMCEffi_AS_inv",
                                    *fMCGenDistTF1_dn_dptG_AS,
                                    *fMCIntrinsicEffiAtAll_dp_dptG,
                                    *tPair_vFits_ptG_i_dp_dr_Axis,
                                    *aAxisPtG,
                                    &CreatePtWeightsInstance(sID + "_lPtWeights",
                                                             true),
                                    &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_AS_inv_ptW);

    tMCEffi_AS_special_ptW = new MCEffi(sID + "_lMCEffi_AS_special",
                                        *fMCGenDistTF1_dn_dptG_AS,
                                        *fMCIntrinsicEffiAtAll_dp_dptG,
                                        *tPair_vFits_ptG_i_dp_dr_Axis,
                                        *aAxisPtG,
                                        &CreatePtWeightsInstance(sID + "_lPtWeights",
                                                                 false),
                                        &vAllDrawableObjects);
    vAllMCEffis.push_back(tMCEffi_AS_special_ptW);

    printf("int InvPtW_main::Initialize(): instance %s: done initializing.\n\n\n", sID.data());
    bInitialized = true;
    return 1;
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
TCanvas &InvPtW::CompareObservables_generic(std::string const &theObservableNameBase,
                                            std::string const &theSelectWhich,
                                            std::string const &theTitleH,
                                            TLegend &theLeg, float theLegTextSize /*= 0.03*/,
                                            float theXmin, float theXmax /*= 0., 10.5 */,
                                            float theYmin, float theYmax /*= 1.e-6, 1.e+4 */,
                                            bool theLogY /*= true*/)
{
    std::string lMethodName(sID +
                            "_CompareObservables_" +
                            theObservableNameBase +
                            "_" + theSelectWhich);

    TCanvas &lCanvas = utils_plotting::GetCanvasWithTH2F(lMethodName + "_canvas" /*theNameC*/,
                                                         theTitleH /*theTitleH*/,
                                                         theXmin, theXmax /*theXmin, theXmax*/,
                                                         theYmin, theYmax /*theYmin, theYmax*/,
                                                         theLogY,
                                                         2000, 1000 /*theWidth, theHeight*/,
                                                         &lMethodName /*theTitleC*/,
                                                         &lMethodName /*theNameH*/,
                                                         0.03 /*theNameTagTextSize*/,
                                                         kGray /*theNameTagTextColor*/,
                                                         1, 2, /*theNx, theNy*/
                                                         1 /*theDrawFirstOnI*/);

    // prepare vector of objects to draw
    /* create these one heap because if not they might cause problems
       with drawing because of lifetime */
    auto &lVectorBundles = *new std::vector<utils_plotting::DrawAndAddBundle>();

    auto &lIterationOuter =
        theSelectWhich == "all" ? *new std::vector<std::string>{"NW", "WW"}
                                : *new std::vector<std::string>{theSelectWhich};

    // outer loop over NW, WW
    for (auto &iOuter : lIterationOuter)
    {
        // inner loop over all MCEffi instances
        for (auto const &iMCEffi : vAllMCEffis)
        {
            if (iMCEffi)
            {
                printf("int InvPtW_main::CompareObservables_generic(): INFO: instance %s\n"
                       "\tiOuter = %s, iMCEffi = %s ... -> ",
                       sID.data(), iOuter.data(), iMCEffi->GetID().data());

                if ((iOuter == "WW") && !iMCEffi->CanRunWithPtWeights())
                {
                    printf("skipping.\n");
                    continue;
                }

                // prepare for constructor of DrawAndAddBundle
                printf("preparing for constructor of DrawAndAddBundle...\n");
                std::string &lFullObservableName =
                    *new std::string(theObservableNameBase + "_" + iOuter);
                std::string lMCID(iMCEffi->GetID());
                std::string &iLegLable = *new std::string(
                    lMCID.erase(0, std::string("lInvPtW_main_lMCEffi_").size()) +
                    "_" +
                    lFullObservableName);

                float lLineWidth = iLegLable.find("D") != std::string::npos ? 1.5
                                   : iOuter == "WW"                         ? 2.
                                                                            : 1.;

                std::string lLegDrawOption(""); // empty enables auto leg draw option

                // create DrawAndAddBundle
                utils_plotting::DrawAndAddBundle lBundle(
                    *iMCEffi->GetObservableObject(lFullObservableName),
                    "same",
                    kGreen + lVectorBundles.size() * 2,
                    lLineWidth,
                    &theLeg,
                    // lInvPtW_main_lMCEffi_ AS_spec
                    iLegLable,
                    lLegDrawOption,
                    theLegTextSize,
                    true /* theDrawLegAlready */);

                // add to vector
                lVectorBundles.push_back(lBundle);
            }
        }
    }

    for (auto iBundle : lVectorBundles)
    {
        utils_plotting::DrawAndAdd(iBundle);
    }

    // utils_plotting::SaveCanvasAs(lCanvas);
    return lCanvas;
}

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
        theComputeInInvariantForm ? hMCGenDist_dn_dptG_inv
                                  : hMCGenDist_dn_dptG,
        theComputeInInvariantForm ? fTargetGenData_dn_dptG
                                  : fTargetGenData_dn_dptG_inv,
        *aAxisPtG);
}

// set first parameter to "auto" for auto-concatenated name
TF1 &InvPtW::FitMCGeneratedParticlesHisto(std::string const &theResultNameInfo,
                                          TH1 &theTH1GenDist_dn_dptG_x, // x = inv or not
                                          bool theTH1IsInvariant) const
{
    // must not be empty
    if (!theResultNameInfo.size())
    {
        printf("InvPtW::FitMCGeneratedParticlesHisto():\n\t"
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

TF1 &InvPtW::FitTrueDetectorIntrinsicMesonEfficiency(std::string fname, Double_t theXmax /*= 10.*/)
{
    TH1 &lTH1Effi = *(TH1 *)utils_files_strings::GetObjectFromPathInFile(fname, "hMyTrueEffiA");
    TCanvas *c3 = new TCanvas("c3", "c3", 2000, 1000);
    TH2 *hd21 = new TH2F("hd1", ";MC pT(GeV);dP/dptG", 1, 0., theXmax, 1, 1e-6, 2e-3);
    hd21->Draw();
    //~ gPad->SetLogy();

    TF1 &lTF1Effi = *new TF1("lTF1Effi", " ( 1e-5 + [0]*x + [1]*x*x + [2]*x*x*x + [3]*x*x*x*x)", 0., theXmax);
    lTF1Effi.SetParameters(1, 1, -1, 0.);
    lTH1Effi.Fit(&lTF1Effi, "N");

    auto &leg = *new TLegend();
    utils_plotting::DrawAndAdd(lTH1Effi, "same", kRed, 1. /* theObjLineWidth */,
                               &leg, "", "l", true /* theDrawLegAlready */);
    utils_plotting::DrawAndAdd(lTF1Effi, "same", kBlue, 1. /* theObjLineWidth */,
                               &leg, "", "l");
    return lTF1Effi;
}
