#pragma once

#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"
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
    ComputeResolutionFits(GCo const &theGCo_h2Res,
                          int theBinStart,
                          int thePtBinMax,
                          int theNR,
                          bool theDrawAllFitsOverlayed,
                          bool thePlotSingles,
                          std::vector<TObject *> *theVAllDrawableObjects);

    utils_fits::TPairFitsWAxis &
    Compute();

    bool GetFitsFromFile(TFile &theFile,
                         int binStart,
                         int ptBinMax,
                         utils_fits::TPairFitsWAxis &thePair);

    std::string const &GetID() const { return sID; }

private:
    std::string const sID;
    std::vector<TObject *> *vAllDrawableObjects;
    GCo const &gCo_h2Res; // from above
    TH2F h2Resolution;
    int const iPtBinStart = 1;
    int const iPtBinMax = 31;
    int iNR;

    // derived
    std::string const sFnameSaveFits;
    bool theDrawAllFitsOverlayed;
    bool thePlotSingles;
};