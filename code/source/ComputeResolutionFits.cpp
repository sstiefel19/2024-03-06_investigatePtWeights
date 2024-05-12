#pragma once
#include "ComputeResolutionFits.h"
#include "source/dN_dptR_integrand.h" // needed for TPairFitsWAxis
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

ComputeResolutionFits::ComputeResolutionFits(TH2 &theH2Resolution,
                                             int theBinStart,
                                             int thePtBinMax,
                                             int theNR,
                                             std::string const &theFnameSave,
                                             bool theDrawAllFitsOverlayed,
                                             bool thePlotSingles)
    : binStart(theBinStart),
      ptBinMax(thePtBinMax),
      nR(theNR),
      fnameSave(theFnameSave),
      theDrawAllFitsOverlayed(theDrawAllFitsOverlayed),
      thePlotSingles(thePlotSingles)
{
    ComputeResolutionFits(theH2Resolution, binStart, ptBinMax, nR, fnameSave, theDrawAllFitsOverlayed, thePlotSingles);
}

static TPairFitsWAxis &
ComputeResolutionFits::Compute(TH2 &theH2Resolution,
                               int binStart,
                               int ptBinMax,
                               int nR = 4,
                               std::string const &fnameSave = "computeResolutionFits.root",
                               bool theDrawAllFitsOverlayed = true,
                               bool thePlotSingles = false)

{
    std::vector<TH1 *> &vHists_ptG_i_dn_dr = *new std::vector<TH1 *>(ptBinMax + 1, static_cast<TH1 *>(nullptr));
    std::vector<TF1 *> &vFits_ptG_i_dp_dr = *new std::vector<TF1 *>(ptBinMax + 1, static_cast<TF1 *>(nullptr));

    // TAxis          &lPtGAxis = *theH2Resolution.GetXaxis();
    TAxis &lPtGAxis = *copyTAxisUpToPt(*theH2Resolution.GetXaxis(), 9.9);
    TPairFitsWAxis &lResult = *new TPairFitsWAxis{vFits_ptG_i_dp_dr, lPtGAxis};

    // plot th2
    TCanvas &cR1 = *new TCanvas("computeResolutionFits_TH2", "computeResolutionFits_TH2", 2000, 1000);
    theH2Resolution.GetXaxis()->SetRangeUser(0., 10.);
    theH2Resolution.Draw("colz");
    gPad->SetLogz();
    saveCanvasAs(cR1);

    // first see if there is a file
    bool lFile = fnameSave.size();
    if (lFile)
    {
        TFile &infile = *new TFile(fnameSave.data(), "READ");
        if (infile.IsOpen())
        {
            printf("a file with name %s exists alerady. Trying to obtain fits from there..\n",
                   fnameSave.data());
            if (getFitsFromFile(infile, binStart, ptBinMax, lResult))
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
        TFile *outfile = new TFile(fnameSave.data(), "RECREATE");
    }

    // make projections to get resolution in each bin of ptG
    Double_t *lLastBinFoundParams = new Double_t[7];
    for (int i = binStart; i <= ptBinMax; ++i)
    {
        cout << i << endl;
        Double_t ptMin = lPtGAxis.GetBinLowEdge(i);
        Double_t ptMax = lPtGAxis.GetBinLowEdge(i + 1);

        // project on r
        TH1 &h1 = *theH2Resolution.ProjectionY(Form("hRes_PtBin%d_%.1f-%.1f_GeV", i, ptMin, ptMax), i, i);
        h1.SetTitle(h1.GetName());
        vHists_ptG_i_dn_dr[i] = &h1;

        // rebin in r
        int lNR = (i == 3) ? 8 : nR;
        TH1 &h1r = *h1.Rebin(lNR, Form("%s_reb%d", h1.GetName(), lNR));
        h1r.SetTitle(h1r.GetName());

        // do the fit
        auto lPairBefNormNorm =
            getCrystallBallFit(h1r,
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
        drawAllFitsOnTop(vFits_ptG_i_dp_dr);
    }
    return lResult;
}

static bool ComputeResolutionFits::
    GetFitsFromFile(TFile &theFile,
                    int binStart,
                    int ptBinMax,
                    TPairFitsWAxis &thePair)
{
    if (!theFile.IsOpen())
    {
        printf("file not open\n");
        return false;
    }
    thePair.second = *(TAxis *)theFile.Get("lPtGAxis");

    for (int i = binStart; i <= ptBinMax; ++i)
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

