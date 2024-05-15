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

    lInvPtW_main.Initialize();
    // lInvPtW_main.Main();

    MCEffi &lEffi_data = lInvPtW_main.GetMCEffi_D();
    // MCEffi  &lEffi_AS_inv = lInvPtW_main.GetMCEffi_AS_inv();
    // MCEffi &lEffi_AS_special = lInvPtW_main.GetMCEffi_AS_special();

    lEffi_data.PlotAll();

    /*
        what do I wanna know?
        in data it holds:
            - Ngen IS correct (by definition)
            - Nrec_i = effi^reco_i * Ngen
            -> Ngen_i = Nrec_i / effi^reco_i          (A)

        overall goal: measure effi^reco_i in MC such  (A) holds

        problem: Nrec_i depends on pT-shape of Ngen because feed-ins != feed-outs
        solution: use pt weights to simulate Ngen shape of data in MC

        goals:
        1) measure effi^reco_i in MC
            a) with Ngen as in data                                     - x
            b) with Ngen flat as in AS MC without pt-weights            - x
            c) with Ngen flat as in AS MC and with pt-weights_inv       - x
            d) with Ngen flat as in AS MC and with pt-weights_special   - y

            * x: Main(true), y: Main(false)
    */

    // notes:
    /*
        lInvPtW_main: constructor sets all information such that an instance can be called with both weights options
    */
}