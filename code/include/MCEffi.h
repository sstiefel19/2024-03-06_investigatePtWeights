// encapsulates a set of:
/*
// intrinsic data members
    1. std::string sID;
    2. TF1 &fMCGenDist_dn_dptG; // the natural distribution of generated particles in this MC.
    TAxis axisPtR;
    // defining since it holds the pt-weights instance
    3. dN_dptR *tdN_dptR_WW_opt; // to simulate the numbers of reconstructed particles with pt-weights

// expressive members
    4. dN_dptR tdN_dptR_NW; // to simulate the numbers of reconstructed particles without pt-weights
    10. TH1 *hMeasuredEffi_NW_2;
    11. TH1 *hMeasuredEffi_WW_2_opt; // optional

*/

#pragma once

#include "TAxis.h"
#include "TObject.h"
#include "dN_dptR.h"
#include "PtWeights.h"

#include <string>
#include <list>
#include <vector>

class TLegend;

class MCEffi
{
public:
    MCEffi(std::string const &_id,
           TF1 &_fGenDist_dn_dptG,
           TF1 &_fEffi_dp_dptG,
           utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
           TAxis &_axisPtR,
           PtWeights *_tPtWeights,
           std::vector<TObject *> *theVAllDrawableObjects);

public:
    std::string const &GetID() const { return sID; }
    TH1 const &GetMeasuredEffi_NW();
    TH1 const *GetMeasuredEffi_WW();
    TH1  &GetMeasuredEffi_NW_clone();
    TH1  *GetMeasuredEffi_WW_clone();

    void PlotAll(TLegend *theLeg = nullptr);

private:
    // intrinsic data members
    std::string sID;
    std::vector<TObject *> *vAllDrawableObjects;
    TAxis axisPtR;
    // defining since it holds the pt-weights instance
    dN_dptR tdN_dptR_NW;      // to simulate the numbers of reconstructed particles without pt-weights
    dN_dptR *tdN_dptR_WW_opt; // to simulate the numbers of reconstructed particles with pt-weights

    // expressive members
    TH1 *hMeasuredEffi_NW_2;
    TH1 *hMeasuredEffi_WW_2_opt; // optional

    // private member functions
    // helper functions - move out?
    TH1 &SampleMeasuredEffi_generic_1D(std::string const &theSuffix,
                                       TF1 const &theNumF,
                                       TF1 const &theDenF,
                                       TAxis const &theAxis) const;

    TH1 &SampleMeasuredEffi_generic_2D(std::string const &theName,
                                       TF2 const &theNumTF2_dN_dptG_dptR,
                                       TF1 const &theDenTF1_dN_dptR) const;

    TH1 &SampleMeasuredEffi_NW_2();
    TH1 *SampleMeasuredEffi_WW_2();
};
