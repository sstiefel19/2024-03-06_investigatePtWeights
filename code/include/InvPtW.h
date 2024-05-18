// purpose: find out how well pt weights work with resolution effect taken into account
#pragma once
#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
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

class MCEffi;

class InvPtW
{
public:
    InvPtW(std::string const &theID,
           std::string const &theFnameInputEffiFit,
           std::string const &theFnameWeightsFile,
           std::string const &theEvCutWeightsFile,
           std::string const &theMeson,
           GCo const &theGCo_h2Resolution,
           int theNR);

    int Initialize();
    int Main();

    TCanvas &CompareMeasuredEfficiencies(TLegend &theLeg);
    TCanvas &CompareGeneratedSpectra(TLegend &theLeg);

    // getters
    TAxis const *GetAxisPtG() const { return aAxisPtG; }
    TF1 const *GetEffiAtAll_dp_dptG() const { return fMCIntrinsicEffiAtAll_dp_dptG; }
    GCo const &GetGCo_h2Resolution() const { return tGCo_forH2Resolution; }

    ComputeResolutionFits const &GetResFitsInstance() const { return tResolutionFits; }
    utils_fits::TPairFitsWAxis const *GetPair_vFits_ptG_i_dp_dr_Axis() const { return tPair_vFits_ptG_i_dp_dr_Axis; }

    std::string const &GetEvCutWeightsFile() const { return sEvCutWeightsFile; }
    std::string const &GetFnameFitOverallEfficiency() const { return sFnameFitOverallEfficiency; }
    std::string const &GetFnameWeightsFile() const { return sFnameWeightsFile; }
    std::string const &GetMeson() const { return sMeson; }

    TH1D const &GetGenDist_dn_dptG_inv() const { return hMCGenDist_dn_dptG_inv; }
    TH1D const &GetGenDist_dn_dptG() const { return hMCGenDist_dn_dptG; }
    TF1 const *GetGenDistTF1_dn_dptG_AS() const { return fMCGenDistTF1_dn_dptG_AS; }
    TF1 const *GetGenDistTF1_dn_dptG_AS_inv() const { return fMCGenDistTF1_dn_dptG_AS_inv; }
    TF1 const &GetTargetGenData_dn_dptG_inv() const { return fTargetGenData_dn_dptG_inv; }
    TF1 const &GetTargetGenData_dn_dptG() const { return fTargetGenData_dn_dptG; }

    MCEffi &GetMCEffi_D() { return *tMCEffi_D; }
    MCEffi &GetMCEffi_AS_inv() { return *tMCEffi_AS_inv_ptW; }
    MCEffi &GetMCEffi_AS_special() { return *tMCEffi_AS_special_ptW; }

    std::string const &GetID() const { return sID; }
    bool IsInitialized() const { return bInitialized; }

    void PlotAll();

private:
    // ==================== member variables ====================
    // intrinsic properties
    // fnames for detector info and weights
    std::string const sID;
    std::string sFnameFitOverallEfficiency;
    std::string sFnameWeightsFile;
    std::string sEvCutWeightsFile;
    std::string sMeson;

    GCo const tGCo_forH2Resolution; // from above
    ComputeResolutionFits tResolutionFits;

    // defining data and MC distributions, will be extracted from above files
    TF1 fTargetGenData_dn_dptG_inv; // from sFnameWeightsFile
    TH1D hMCGenDist_dn_dptG_inv;    // from sFnameWeightsFile

    // derived information
    TH1D hMCGenDist_dn_dptG;
    TF1 fTargetGenData_dn_dptG;

    // members that can't be initialized right away
    TF1 *fMCIntrinsicEffiAtAll_dp_dptG;
    utils_fits::TPairFitsWAxis *tPair_vFits_ptG_i_dp_dr_Axis;
    TAxis *aAxisPtG; // gets intialized in ParametrizeEfficiencies

    TF1 *fMCGenDistTF1_dn_dptG_AS;
    TF1 *fMCGenDistTF1_dn_dptG_AS_inv;

    MCEffi *tMCEffi_D;              // Ngen shape as in data, no pt-weights
    MCEffi *tMCEffi_AS_inv_ptW;     // Ngen shape as in AS MC, can run with and without pt-weights
    MCEffi *tMCEffi_AS_special_ptW; // Ngen shape as in AS MC, can run with and without pt-weights

    // accounting
    bool bInitialized;
    std::vector<MCEffi *> vAllMCEffis;
    std::vector<TObject *> vAllDrawableObjects;

    // ====================== private member functions =========================
    TCanvas &CompareObservables_generic(std::string const &theObservableNameBase,
                                        std::string const &theSelectWhich,
                                        std::string const &theTitleH,
                                        TLegend &theLeg, float theLegTextSize = 0.03,
                                        float theXmin = 0., float theXmax = 10.5,
                                        float theYmin = 1.e-6, float theYmax = 1.e+4,
                                        bool theLogY = true);

    // set first parameter to "auto" for auto-concatenated name
    TF1 &FitMCGeneratedParticlesHisto(std::string const &theResultNameInfo,
                                      TH1 &theTH1GenDist_dn_dptG_x, // x = inv or not
                                      bool theTH1IsInvariant);

    // detector parametrizations
    TF1 &GetMesonEfficiency(std::string fname, Double_t theXmax = 10.);

    // pt weights related
    PtWeights &CreatePtWeightsInstance(std::string const &theID,
                                       bool theComputeInInvariantForm);
};
