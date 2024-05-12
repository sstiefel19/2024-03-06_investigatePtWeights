#pragma once

#include "/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"

class TAxis;
class TCanvas;
class TH2;
class TF1;
class TFile;
class TLegend;

#include <string>

class ComputeResolutionFits
{
public:
    ComputeResolutionFits(TH2 &theH2Resolution,
                          int theBinStart,
                          int thePtBinMax,
                          int theNR,
                          std::string const &theFnameSave,
                          bool theDrawAllFitsOverlayed,
                          bool thePlotSingles);

    utils_fits::TPairFitsWAxis &
    Compute(TH2 &theH2Resolution,
            int binStart,
            int ptBinMax,
            int nR,
            std::string const &fnameSave,
            bool theDrawAllFitsOverlayed = true,
            bool thePlotSingles = false);

    bool GetFitsFromFile(TFile &theFile,
                         int binStart,
                         int ptBinMax,
                         utils_fits::TPairFitsWAxis &thePair);

private:
    TH2 &h2Resolution;
    int binStart;
    int ptBinMax;
    int nR;
    std::string const fnameSave;
    bool theDrawAllFitsOverlayed;
    bool thePlotSingles;
};