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
    InvPtW_main(std::string const &theMeson,
                std::string const &theCent,
                std::string const &theCentAS,
                std::string const &theFnameAS,
                std::string const &theFnameInputEffiFit,
                std::string const &theFnameWeightsFile,
                int thePtBinStart,
                int thePtBinMax,
                bool theGenDistTH1IsInvariant,
                bool theMultiplyResultTF1ByX);

    TF1 &GetMesonEfficiency(std::string fname, Double_t theXmax = 10.);

    // set first parameter to "auto" for auto-concatenated name
    TF1 &FitGenDistHisto(std::string const &theResultNameInfo,
                         TH1 &theTH1GenDist_dn_dptG,
                         bool theTH1IsInvariant,
                         TAxis const &thePtGaxis,
                         bool theMultiplyResultTF1ByX,
                         bool &theResultIsInvariant_out);

    // detector parametrizations
    utils_fits::TPairFitsWAxis FitDetector(const std::string &theFnameInputEffiFit,
                               const std::string &theFnameResFits,
                               TH2 &theH2Resolution,
                               int thePtBinStart,
                               int thePtBinMax,
                               int theNRebin_r,
                               bool theDrawAllFitsOverlayed,
                               bool thePlotSingles,
                               TF1 *&theEffiAtAll_dp_dptG_out);
    
    PtWeights &SetupWeightsInstance(std::string const &theID,
                                    bool theComputedInInvariantForm,
                                    TAxis const &thePtGaxis);

    int Main();

private:
    int ptBinStart;
    int ptBinMax;

    bool lGenDistTH1IsInvariant;
    bool lMultiplyResultTF1ByX;

    std::string meson;
    std::string cent;
    std::string centAS;
    std::string fnameAS;
    // std::string fnameAS("/trains/2024-01-26_LHC24a1_QA_noPtW/GCo_997_both.root");

    std::string lFnameInputEffiFit;
    std::string lFnameResFits;
    std::string fnameWeightsFile;
    GCo gAS;
    TH2F h2Resolution;
};
