#include "TROOT.h"

void run()
{   
    bool b1 = gROOT->ProcessLine(".L source/PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L source/dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L source/dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L source/MCEffi.cpp+");

    gROOT->ProcessLine(".x investigatePtWeights_wResolutionEffects.C");
}