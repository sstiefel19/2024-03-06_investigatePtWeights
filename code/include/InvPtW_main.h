// purpose: find out how well pt weights work with resolution effect taken into account
#pragma once
#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
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
                std::string const &theMeson,
                std::string const &theCent,
                std::string const &theCentAS,
                std::string const &theFnameAS);

    TF1 &GetMesonEfficiency(std::string fname, Double_t theXmax = 10.);

    // set first parameter to "auto" for auto-concatenated name
    TF1 &FitGenDistHisto(std::string const &theResultNameInfo,
                         TH1 &theTH1GenDist_dn_dptG,
                         bool theTH1IsInvariant,
                         bool theMultiplyResultTF1ByX,
                         bool &theResultIsInvariant_out);

    // detector parametrizations
    utils_fits::TPairFitsWAxis &FitDetector(int theNRebin_r,
                                            bool theDrawAllFitsOverlayed,
                                            bool thePlotSingles);

    int Main(bool theUseInvariantForm);

private:
    PtWeights &CreatePtWeightsInstance(std::string const &theID,
                                       bool theComputeInInvariantForm);

    // ==================== defining & intrinsic properties ====================
    // fnames for detector info and weights
    std::string sFnameInputEffiFit;
    std::string sFnameWeightsFile;

    // for the fitting
    int iPtBinStart = 1;
    int iPtBinMax = 31;

    // input filenames
    std::string sMeson;
    std::string sCent;
    std::string sCentAS;
    std::string sFnameAS;
    // std::string fnameAS("/trains/2024-01-26_LHC24a1_QA_noPtW/GCo_997_both.root");

    // ======================= derived properties ==============================
    std::string sFnameResFits;

    // holding detector reponse ptR vs. ptG
    TH2  *h2Resolution;

    // defining data and MC distributions, will be extracted from above files
    TH1 &hGenDist_dn_dptG_inv;
    TF1 &fTargetGenData_dn_dptG_inv;

    // from input derived information
    TH1  *hGenDist_dn_dptG;
    TF1  *fTargetGenData_dn_dptG;

    // measured output variables
    TF1 *fEffiAtAll_dp_dptG;

    // ======================= helper functions ===============================
    // helper members for accessing GammaConvV1 output files
    std::string const id;
    GCo *gConvV1_AS;
    TAxis *aAxisPtG; // gets intialized in FitDetector

};
