#include "source/dN_dptR_integrand.h" // needed for utils_fits::TPairFitsWAxis

// extended crystal Ball function
/*
 consists on a gaussian 'core', extended continuously by two power law tails on left and right side
 parameters:
 mean and sigma for the gaussian core
 alpha1 the point (in standard deviations from the mean) at which the left side tail starts
 n1 the exponent of the left side tail
 alpha2 and n2: the same for the right side
 */
Double_t CrystallBall2(Double_t x, Double_t mean, Double_t sigma, Double_t alpha1, Double_t n1, Double_t alpha2, Double_t n2)
{
    auto SQUARE = [](Double_t x)
    { return x * x; };

    Double_t t = (x - mean) / sigma;
    if (t < -alpha1)
    {
        Double_t a = TMath::Power(n1 / alpha1, n1) * TMath::Exp(-SQUARE(alpha1) / 2);
        Double_t b = n1 / alpha1 - alpha1;
        return a / TMath::Power(b - t, n1);
    }
    else if (t > alpha2)
    {

        Double_t a = TMath::Power(n2 / alpha2, n2) * TMath::Exp(-SQUARE(alpha2) / 2);
        Double_t b = n2 / alpha2 - alpha2;
        return a / TMath::Power(b + t, n2);
    }
    else
        return TMath::Exp(-SQUARE(t) / 2);
}

Double_t CrystallBall(Double_t *x, Double_t *par)
{
    return CrystallBall2(*x, par[0], par[1], par[2], par[3], par[4], par[5]) * par[6];
}

void drawAllFitsOnTop(std::vector<TF1 *> &theFits)
{

    auto c = new TCanvas("drawAllFitsOnTop", "drawAllFitsOnTop", 2000, 1000);
    TH2 *hd21 = new TH2F("fits", ";r (GeV);dp/dr (1./GeV)", 1, -1., 1., 1, 1e-3, 4e1);
    hd21->Draw();
    gPad->SetLogy();

    auto leg = new TLegend(.8, .1, .9, .9);
    printf("drawFits(): theFits.size() = %zu\n", theFits.size());
    for (size_t i = 1; i < theFits.size(); ++i)
    {
        auto f = theFits.at(i);
        if (f)
        {
            f->SetLineColor(i);
            leg->AddEntry(f, Form("fit%zu", i), "l");
            f->Draw("same");
        }
    }
    leg->Draw("same");
    utils_plotting::SaveCanvasAs(*c);
}

// fit, normalized fit
std::pair<TF1 &, TF1 &>
getCrystallBallFit(TH1 &theH,
                   bool theGaussOnly = false,
                   Double_t fitXmin = -.6,
                   Double_t fitXmax = .2,
                   Double_t xMin = -1.,
                   Double_t xMax = 1.,
                   Double_t *theInitialPars = nullptr)
{

    const char *lNameFitNorm = Form("%s_fit_Norm", theH.GetName());
    const char *lNameFitBefNorm = Form("%s_fit_BefNorm", theH.GetName());
    if (!theH.GetEntries())
    {
        return std::pair<TF1 &, TF1 &>{*new TF1(lNameFitBefNorm, "0.5", xMin, xMax),
                                       *new TF1(lNameFitNorm, "0.5", xMin, xMax)};
    }
    auto &f_BefNorm = *new TF1(lNameFitBefNorm, CrystallBall, xMin, xMax, 7);

    Double_t lXYmax = theH.GetBinCenter(theH.GetMaximumBin());
    Double_t lYmax = theH.GetBinContent(theH.GetMaximumBin());

    std::vector<Double_t> lInitialPars({lXYmax,  // 0 mean gauss
                                        .1,      // 1 sig gauss
                                        1.,      // 2 alpha1 the point (in standard deviations from the mean) at which the left side tail starts
                                        2.,      // 3 n1 the exponent of the left side tail
                                        1.,      // 4 alpha rhs
                                        2.,      // 5 n1 rhs
                                        lYmax}); // 6 amplitude
    f_BefNorm.SetParameters(theInitialPars ? theInitialPars : lInitialPars.data());

    // special cases
    if (theGaussOnly)
    {
        f_BefNorm.FixParameter(2, 1000.);
        f_BefNorm.FixParameter(3, 0.);
        f_BefNorm.FixParameter(4, 1000.);
        f_BefNorm.FixParameter(5, 0.);
    }
    theH.Fit(&f_BefNorm, "NQ", "", fitXmin, fitXmax);
    Double_t lIntegral = f_BefNorm.Integral(xMin, xMax);
    Double_t lNewAmp = f_BefNorm.GetParameter(6) / lIntegral;

    TF1 &f_Norm = *new TF1(lNameFitNorm, CrystallBall, xMin, xMax, 7);
    f_Norm.SetParameters(f_BefNorm.GetParameters());
    f_Norm.SetParameter(6, lNewAmp);
    Double_t lIntegralN = f_Norm.Integral(xMin, xMax);
    if (abs(lIntegralN - 1.) > 1e-5)
    {
        printf("getCrystallBallFit(): WARNING: Integral of %s should be 1 but is %f\n",
               lNameFitNorm, lIntegralN);
    }
    return std::pair<TF1 &, TF1 &>{f_BefNorm, f_Norm};
}

bool getFitsFromFile(TFile &theFile,
                     int binStart,
                     int ptBinMax,
                     utils_fits::TPairFitsWAxis &thePair)
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

utils_fits::TPairFitsWAxis &
computeResolutionFits(TH2 &theH2Resolution,
                      int binStart,
                      int ptBinMax,
                      int nR = 4,
                      std::string const &sFnameSaveFits = "computeResolutionFits.root",
                      bool theDrawAllFitsOverlayed = true,
                      bool thePlotSingles = false)
{
    std::vector<TH1 *> &vHists_ptG_i_dn_dr = *new std::vector<TH1 *>(ptBinMax + 1, static_cast<TH1 *>(nullptr));
    std::vector<TF1 *> &vFits_ptG_i_dp_dr = *new std::vector<TF1 *>(ptBinMax + 1, static_cast<TF1 *>(nullptr));

    // TAxis          &lPtGAxis = *theH2Resolution.GetXaxis();
    TAxis &lPtGAxis = *utils_computational::CopyTAxisUpToPt(*theH2Resolution.GetXaxis(), 9.9);
    utils_fits::TPairFitsWAxis &lResult = *new utils_fits::TPairFitsWAxis{vFits_ptG_i_dp_dr, lPtGAxis};

    // plot th2
    TCanvas &cR1 = *new TCanvas("computeResolutionFits_TH2", "computeResolutionFits_TH2", 2000, 1000);
    theH2Resolution.GetXaxis()->SetRangeUser(0., 10.);
    theH2Resolution.Draw("colz");
    gPad->SetLogz();
    utils_plotting::SaveCanvasAs(cR1);

    // first see if there is a file
    bool lFile = sFnameSaveFits.size();
    if (lFile)
    {
        TFile &infile = *new TFile(sFnameSaveFits.data(), "READ");
        if (infile.IsOpen())
        {
            printf("a file with name %s exists alerady. Trying to obtain fits from there..\n",
                   sFnameSaveFits.data());
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
        TFile *outfile = new TFile(sFnameSaveFits.data(), "RECREATE");
    }

    // make projections to get resolution in each bin of ptG
    Double_t *lLastBinFoundParams = new Double_t[7];
    for (int i = binStart; i <= ptBinMax; ++i)
    {
        std::cout << i << std::endl;
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
                utils_plotting::DrawAndAdd(h1r, "", kBlue, leg, "histo");
                utils_plotting::DrawAndAdd(fitBN, "same", kRed, leg, "fit");
                utils_plotting::DrawAndAdd(fitN, "same", kGreen, leg, "fit norm.");
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
    utils_plotting::SaveCanvasAs(cR);
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
