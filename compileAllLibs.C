#include "TROOT.h"

void compileAllLibs()
{   
    bool u2 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/utils_sstiefel_2024.cpp+");
    bool u1 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/utils_TF1.cpp+");
    bool u0 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/GCo.cpp+");
    bool b1 = gROOT->ProcessLine(".L source/PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L source/dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L source/dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L source/MCEffi.cpp+");
}