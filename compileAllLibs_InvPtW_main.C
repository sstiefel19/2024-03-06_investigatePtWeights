#include "TROOT.h"

void compileAllLibs_InvPtW_main()
{
    bool u1 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_files_strings.cpp+");
    bool u2 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/GCo.cpp+");
    bool u4 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_plotting.cpp+");
    bool u3 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_fits.cpp+");
    bool u5 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_computational.cpp+");
    bool u6 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_TF1.cpp+");
    bool u7 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_TH1.cpp+");

    bool b0 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/computeResolutionFits.cpp+");
    bool b1 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/MCEffi.cpp+");
    bool b5 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/InvPtW.cpp+");
}