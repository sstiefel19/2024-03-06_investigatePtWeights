#include "TROOT.h"

void run_testTH1continuosEvaluationMethods()
{   
    bool u2 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/utils_sstiefel_2024.cpp+");
    bool u0 = gROOT->ProcessLine(".L /analysisSoftware/SupportingMacros/GCo.cpp+");

    gROOT->ProcessLine(".x testTH1continuosEvaluationMethods.C");
}