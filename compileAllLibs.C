#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/2024/2024-03-06_investigatePtWeights/code/include/InvPtW_main.h"
#include "TROOT.h"

void compileAllLibs()
{   
    bool u0 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_sstiefel_2024.cpp+");
    bool u1 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_fits.cpp+");
    bool u2 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/GCo.cpp+");
    bool u3 = gROOT->ProcessLine(".L /analysisSoftware/utils_sstiefel_2024/src/utils_TF1.cpp+");
    bool b0 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/computeResolutionFits.cpp+");
    bool b1 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/MCEffi.cpp+");
    bool b5 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/code/src/InvPtW_main.cpp+");
}