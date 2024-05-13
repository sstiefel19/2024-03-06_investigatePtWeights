// purpose: find out how well pt weights work with resolution effect taken into account
#pragma once
#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"

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

class InvPtW_main
{
public:
    InvPtW_main(std::string const &theID,
                std::string const &theFnameInputEffiFit,
                std::string const &theFnameWeightsFile,
                std::string const &theFnameAS,
                std::string const &theMeson,
                std::string const &theEvCut,
                std::string const &theEvCutAS,
                std::string const &theMainDirAS,
                std::string const &thePhotMesCutNo);

    int Main(bool theUseInvariantForm);

private:
    // detector parametrizations
    utils_fits::TPairFitsWAxis &FitDetector(int theNRebin_r,
                                            bool theDrawAllFitsOverlayed,
                                            bool thePlotSingles);

    // set first parameter to "auto" for auto-concatenated name
    TF1 &FitMCGeneratedParticlesHisto(std::string const &theResultNameInfo,
                                      TH1 &theTH1GenDist_dn_dptG,
                                      bool theTH1IsInvariant,
                                      bool theMultiplyResultTF1ByX);

    // detector parametrizations
    TF1 &GetMesonEfficiency(std::string fname, Double_t theXmax = 10.);

    // pt weights related
    PtWeights &CreatePtWeightsInstance(std::string const &theID,
                                       bool theComputeInInvariantForm);

    // ==================== defining & intrinsic properties ====================
    // fnames for detector info and weights
    std::string const id;
    std::string sFnameInputEffiFit;
    std::string sFnameWeightsFile;

    // input physiscs configs filenames
    std::string sFnameAS;
    std::string sMeson;
    std::string sEvCut;
    std::string sEvCutAS;
    GCo const gConvV1_AS;

    // for the fitting
    int iPtBinStart = 1;
    int iPtBinMax = 31;
    std::string sFnameResFits;

    // ======================= derived properties ==============================

    // holding detector reponse ptR vs. ptG
    TH2F h2Resolution;

    // defining data and MC distributions, will be extracted from above files
    TF1 &fTargetGenData_dn_dptG_inv;
    TH1 &hGenDist_dn_dptG_inv;

    // from input derived information
    TH1D hGenDist_dn_dptG;
    TF1 fTargetGenData_dn_dptG;

    // measured output variables
    TF1 *fEffiAtAll_dp_dptG;

    // ======================= helper functions ===============================
    // helper members for accessing GammaConvV1 output files
    TAxis *aAxisPtG; // gets intialized in FitDetector
};
