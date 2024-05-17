#include "/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_plotting.h"
#include "/2024/2024-03-06_investigatePtWeights/code/include/InvPtW.h"
#include "TROOT.h"

void run_InvPtW_main()
{
    // gROOT->ProcessLine(".L compileAllLibs.C");

    printf("run_InvPtW_main.C::run_InvPtW_main(): Start\n");
    int lNR = 4;
    GCo const &lGCo_AS = *new GCo("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_997.root",
                                  "GammaConvV1_997/",
                                  "10130023",
                                  "_0d200009ab770c00amd0404000_0152101500000000",
                                  true /* _keepFileOpen */);

    printf("run_InvPtW_main.C::run_InvPtW_main(): calling InvPtW_main::InvPtW()\n");
    InvPtW &lInvPtW_main = *new InvPtW(
        "lInvPtW_main",
        "input_root/input_for_effi-fit_101.root",
        "/2024/2024-01-29_determinePtWeights/newUploadedFiles/MCSpectraInputPbPb_Stephan_it0b_with24a1.root",
        "10130023",
        "Pi0",
        lGCo_AS,
        lNR);

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
    // Initialize && PlotAll
    lInvPtW_main.Main();
    utils_plotting::SaveCanvasAs(lInvPtW_main.CompareMeasuredEfficiencies());
}