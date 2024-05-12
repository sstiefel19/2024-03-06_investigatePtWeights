#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2024.h"
#include "source/InvPtW_main.h"
#include "TROOT.h"

void run_InvPtW_main()
{
    bool u0 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/GCo.cpp+");
    bool u2 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/utils_sstiefel_2024.cpp+");
    bool u1 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/utils_TF1.cpp+");
    bool b0 = gROOT->ProcessLine(".L computeResolutionFits.h+");
    bool b1 = gROOT->ProcessLine(".L source/PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L source/dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L source/dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L source/MCEffi.cpp+");
    bool b5 = gROOT->ProcessLine(".L source/InvPtW_main.cpp+");

    InvPtW_main &lInvPtW_main = *new InvPtW_main(
        "Pi0",
        "10130053",
        "10130023",
        "/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root",
        "input_for_effi-fit_101.root",
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root",
        1,
        31,
        true,
        true);
    lInvPtW_main.Main();

}