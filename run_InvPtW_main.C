#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/2024/2024-03-06_investigatePtWeights/code/include/InvPtW_main.h"
#include "TROOT.h"

void run_InvPtW_main()
{
    // gROOT->ProcessLine(".L compileAllLibs.C");

    int lNR = 4;
    GCo const lGCo_AS("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root",
                      "GammaConvV1_997/",
                      "10130023",
                      "_0d200009ab770c00amd0404000_0152101500000000",
                      false /* _keepFileOpen */);

    InvPtW_main &lInvPtW_main = *new InvPtW_main(
        "lInvPtW_main",
        "input_root/input_for_effi-fit_101.root",
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root",
        "10130023",
        "Pi0",
        lGCo_AS,
        lNR);

    lInvPtW_main.Main(true /* theUseInvariantForm */);

    lInvPtW_main.Main(false /* theUseInvariantForm */);
}