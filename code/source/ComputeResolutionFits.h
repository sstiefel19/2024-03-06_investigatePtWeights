#pragma once
#include "source/dN_dptR_integrand.h" // needed for TPairFitsWAxis
#include "TCanvas.h"
#include "TLegend.h"

#include <string>

class TH2;
class TF1;
class TFile;
class TPairFitsWAxis;


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

    

    private:
    int binStart;
    int ptBinMax;
    int nR;
    std::string const fnameSave;
    bool theDrawAllFitsOverlayed;
    bool thePlotSingles;
};