#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/2024/2024-03-06_investigatePtWeights/code/include/InvPtW_main.h"
#include "TROOT.h"

void run_InvPtW_main()
{
    // gROOT->ProcessLine(".L compileAllLibs.C");

    InvPtW_main &lInvPtW_main = *new InvPtW_main(
        "lInvPtW_main",
        "input_root/input_for_effi-fit_101.root",
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root",
        "/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root",
        "Pi0",
        "10130053",
        "10130023",
        "GammaConvV1_997/",
        "_0d200009ab770c00amd0404000_0152101500000000");

    lInvPtW_main.Main(true /* theUseInvariantForm */);
}