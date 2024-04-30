#include "TROOT.h"

void run()
{   
    bool b1 = gROOT->ProcessLine(".L PtWeights.cpp+");
    bool b2 = gROOT->ProcessLine(".L dN_dptR_integrand.cpp+");
    bool b3 = gROOT->ProcessLine(".L dN_dptR.cpp+");
    bool b4 = gROOT->ProcessLine(".L MCEffi.cpp+");

    gROOT->ProcessLine(".x investigatePtWeights_wResolutionEffects.C");
    // if (b1 && b2 && b3 && b4)
    // {
    //     gROOT->ProcessLine(".x investigatePtWeights_wResolutionEffects.C");
    // }
    // else 
    // {
    //     printf("Error loading the libraries\n");
    // }
}