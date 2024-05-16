#pragma once
#include "../include/ComputeResolutionFits.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

ComputeResolutionFits::ComputeResolutionFits(
    GCo const &theGCo_h2Res,
    int thePtBinStart,
    int thePtBinMax,
    int theNR,
    bool theDrawAllFitsOverlayed,
    bool thePlotSingles,
    std::vector<TObject *> *theVAllDrawableObjects /*= nullptr*/)
    : gCo_h2Res(theGCo_h2Res),
      vAllDrawableObjects(
          theVAllDrawableObjects ? theVAllDrawableObjects : new std::vector<TObject *>()),
      h2Resolution(*(TH2F *)gCo_h2Res.GetFromTrue("ESD_TruePrimaryPi0_MCPt_ResolPt")),
      iPtBinStart(thePtBinStart),
      iPtBinMax(thePtBinMax),
      iNR(theNR),
      sFnameSaveFits(Form("input_root/%s_resolutionFits_%d-%d.root",
                          gCo_h2Res.evCut.data(),
                          thePtBinStart,
                          thePtBinMax)),
      theDrawAllFitsOverlayed(theDrawAllFitsOverlayed),
      thePlotSingles(thePlotSingles)
{
    printf("ComputeResolutionFits::ComputeResolutionFits(): created instance:\n"
           "\th2Resolution: %s iPtBinStart: %d, iPtBinMax: %d\n",
           h2Resolution.GetName(), iPtBinStart, iPtBinMax);
}

utils_fits::TPairFitsWAxis &
ComputeResolutionFits::Compute()

{
    printf("ComputeResolutionFits::Compute(): iPtBinStart: %d, iPtBinMax: %d\n", iPtBinStart, iPtBinMax);

    std::vector<TH1 *> &vHists_ptG_i_dn_dr = *new std::vector<TH1 *>(iPtBinMax + 1, static_cast<TH1 *>(nullptr));
    std::vector<TF1 *> &vFits_ptG_i_dp_dr = *new std::vector<TF1 *>(iPtBinMax + 1, static_cast<TF1 *>(nullptr));

    // TAxis          &lPtGAxis = *h2Resolution.GetXaxis();
    TAxis &lPtGAxis = *copyTAxisUpToPt(*h2Resolution.GetXaxis(), 9.9);
    utils_fits::TPairFitsWAxis &lResult = *new utils_fits::TPairFitsWAxis{vFits_ptG_i_dp_dr, lPtGAxis};

    // plot th2
    TCanvas &cR1 = *new TCanvas("computeResolutionFits_TH2", "computeResolutionFits_TH2", 2000, 1000);
    h2Resolution.GetXaxis()->SetRangeUser(0., 10.);
    h2Resolution.Draw("colz");
    gPad->SetLogz();
    saveCanvasAs(cR1);

    std::cout << "ComputeResolutionFits::Compute(): iPtBinStart: " << iPtBinStart << ", iPtBinMax: " << iPtBinMax << " " << sFnameSaveFits.data() << std::endl;

    // first see if there is a file
    bool lFile = sFnameSaveFits.size();
    if (lFile)
    {
        TFile &infile = *new TFile(sFnameSaveFits.data(), "READ");
        if (infile.IsOpen())
        {
            printf("a file with name %s exists alerady. Trying to obtain fits from there..\n",
                   sFnameSaveFits.data());
            if (GetFitsFromFile(infile, iPtBinStart, iPtBinMax, lResult))
            {
                printf(" worked :)\n");
                return lResult;
            }
            printf("did not work. recompute.\n");
        }
        infile.Close();
    }

    TCanvas &cR = *new TCanvas("computeResolutionFits_mos", "computeResolutionFits_mos", 2000, 1000);
    cR.Divide(5, 6);

    if (lFile)
    {
        TFile *outfile = new TFile(sFnameSaveFits.data(), "RECREATE");
    }

    // make projections to get resolution in each bin of ptG
    Double_t *lLastBinFoundParams = new Double_t[7];
    for (int i = iPtBinStart; i <= iPtBinMax; ++i)
    {
        cout << i << endl;
        Double_t ptMin = lPtGAxis.GetBinLowEdge(i);
        Double_t ptMax = lPtGAxis.GetBinLowEdge(i + 1);

        // project on r
        TH1 &h1 = *h2Resolution.ProjectionY(Form("hRes_PtBin%d_%.1f-%.1f_GeV", i, ptMin, ptMax), i, i);
        h1.SetTitle(h1.GetName());
        vHists_ptG_i_dn_dr[i] = &h1;

        // rebin in r
        int lNR = (i == 3) ? 8 : iNR;
        TH1 &h1r = *h1.Rebin(lNR, Form("%s_reb%d", h1.GetName(), lNR));
        h1r.SetTitle(h1r.GetName());

        // do the fit
        auto lPairBefNormNorm = utils_fits::GetCrystallBallFit(
            h1r,
            i == 3, /*theGaussOnly*/
            -.6, .2, -1., 1.,
            i == 7 ? lLastBinFoundParams : nullptr);

        TF1 &fitBN = lPairBefNormNorm.first;
        TF1 &fitN = lPairBefNormNorm.second;
        vFits_ptG_i_dp_dr[i] = &fitN;

        fitBN.GetParameters(lLastBinFoundParams);

        // drawing
        {
            auto drawAll = [&](TLegend *leg)
            {
                drawAndAdd(h1r, "", kBlue, leg, "histo");
                drawAndAdd(fitBN, "same", kRed, leg, "fit");
                drawAndAdd(fitN, "same", kGreen, leg, "fit norm.");
                leg->Draw("same");
            };

            cR.cd(i);
            gPad->SetLogy();
            drawAll(new TLegend(.14, .62, .31, .88, ""));
            if (thePlotSingles)
            {
                auto c = new TCanvas(h1r.GetName(), h1r.GetName(), 1000, 1000);
                gPad->SetLogy();
                drawAll(new TLegend(.14, .62, .31, .88, ""));
            }
        }
        if (lFile)
        {
            h1.Write();
            h1r.Write();
            fitBN.Write();
            fitN.Write(Form("fitN_ptGbin_%d", i));
        }
    }
    saveCanvasAs(cR);
    if (lFile)
    {
        lPtGAxis.Write("lPtGAxis");
        cR.Write();
    }
    if (theDrawAllFitsOverlayed)
    {
        utils_fits::DrawAllFitsOnTop(vFits_ptG_i_dp_dr);
    }
    return lResult;
}

bool ComputeResolutionFits::
    GetFitsFromFile(TFile &theFile,
                    int iPtBinStart,
                    int iPtBinMax,
                    utils_fits::TPairFitsWAxis &thePair)
{
    if (!theFile.IsOpen())
    {
        printf("file not open\n");
        return false;
    }
    thePair.second = *(TAxis *)theFile.Get("lPtGAxis");

    for (int i = iPtBinStart; i <= iPtBinMax; ++i)
    {

        const char *fitname = Form("fitN_ptGbin_%d", i);
        TF1 *fit = (TF1 *)theFile.Get(fitname);
        if (fit)
        {
            printf("found %s \n", fitname);
        }
        thePair.first[i] = fit;
    }
    return true;
}
