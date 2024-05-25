#include "TROOT.h"
#include "TSystem.h"

void both_compile_and_run_main()
{
    gSystem->Exec("touch log_both_compile_and_run_main.txt");
    bool b0 = gSystem->Exec("root -L /2024/2024-03-06_investigatePtWeights/compileAllLibs_InvPtW_main.C >> log_both_compile_and_run_main.txt");
    // bool b1 = gSystem->Exec("root -x /2024/2024-03-06_investigatePtWeights/run_InvPtW_main.C >> log_both_compile_and_run_main.txt");
    bool b1 = gROOT->ProcessLine(".x /2024/2024-03-06_investigatePtWeights/run_InvPtW_main.C");
}
