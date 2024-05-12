#pragma once;

class TF1;
class TH1;
class TFile;
class TPairFitsWAxis;

class<utility>;
class<vector>;

class utils_fits
{
public:
    static double CrystallBall(double *x, double *par);

    static double CrystallBall2(double x,
                                double mean,
                                double sigma,
                                double alpha1,
                                double n1,
                                double alpha2,
                                double n2);
    static void DrawAllFitsOnTop(std::vector<TF1 *> &theFits);

    static std::pair<TF1 &, TF1 &>
    GetCrystallBallFit(TH1 &theH,
                       bool theGaussOnly,
                       double fitXmin,
                       double fitXmax,
                       double xMin,
                       double xMax,
                       double *theInitialPars);

    static bool GetFitsFromFile(TFile &theFile,
                                int binStart,
                                int ptBinMax,
                                TPairFitsWAxis &thePair);
};
