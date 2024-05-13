#include "TROOT.h"

int start_compile_run_main()
{
    bool b0 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/compileAllLibs.C");
    bool b1 = gROOT->ProcessLine(".L /2024/2024-03-06_investigatePtWeights/run_InvPtW_main.C");

    return b0 && b1;
}