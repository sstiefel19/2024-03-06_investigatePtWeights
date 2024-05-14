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
                std::string const &theEvCutWeightsFile,
                std::string const &theMeson,
                GCo const &theGCo_h2Resolution,
                int theNR);

    int Main(bool theUseInvariantForm);

private:
    // member functions
    utils_fits::TPairFitsWAxis &ParametrizeEfficiencies();

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

    // ==================== member variables ====================
    // fnames for detector info and weights
    // intrinsic properties
    std::string const sID;
    std::string sFnameFitOverallEfficiency;
    std::string sFnameWeightsFile;
    std::string sEvCutWeightsFile;
    std::string sMeson;
    GCo const gGCo_h2Resolution; // from above
    ComputeResolutionFits cRF;

    // defining data and MC distributions, will be extracted from above files
    TF1 fTargetGenData_dn_dptG_inv; // from sFnameWeightsFile
    TH1D hGenDist_dn_dptG_inv;      // from sFnameWeightsFile

    // from input derived information
    TH1D hGenDist_dn_dptG;
    TF1 fTargetGenData_dn_dptG;

    // measured output variables
    TF1 *fEffiAtAll_dp_dptG;

    // ======================= helper functions ===============================
    // helper members for accessing GammaConvV1 output files
    TAxis *aAxisPtG; // gets intialized in ParametrizeEfficiencies
};
