#pragma once

#include "TAxis.h"
#include "TList.h"
#include "dN_dptR.h"
#include "PtWeights.h"
#include <string>

class TLegend;

class MCEffi
{
public:
    MCEffi(std::string const &_id,
           TF1 &_fGenDist_dn_dptG,
           TF1 &_fEffi_dp_dptG,
           utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
           TAxis &_axisPtR,
           PtWeights *_tPtWeights = nullptr);

    // integrate TF2
    TH1 &SampleMeasuredEffi_NW_2(Color_t theLineColor = kBlue);
    TH1 *SampleMeasuredEffi_WW_2(Color_t theLineColor = kBlue);

    std::string const &GetID() const { return id; }

    void PlotAll(TLegend *theLeg = nullptr);

private:
    // private member functions
    TH1 &SampleTH1FromTF2(TF2 &theTF2_d_dptG_dptR, std::string name = "");

    TH1 &SampleMeasuredEffi_generic(std::string const &theSuffix,
                                    TF1 &theNumF,
                                    TF1 &theDenF,
                                    TAxis const &theAxis) const;

    // defining data members
    std::string id;
    TF1 &fGenDist_dn_dptG; // todo make const ?!
    dN_dptR tdN_dptR_NW;
    TAxis axisPtR;

    // optional members
    PtWeights *tPtWeights;
    dN_dptR *tdN_dptR_WW;

    // helper members
    TF1 *fGenDist_dn_dptG_WW;
    TH1 *hMCGen_NW;
    TH1 *hMCRec_NW;
    TH1 *hEffiRec_NW;

    TH1 *hMCGen_WW;
    TH1 *hMCRec_WW;
    TH1 *hEffiRec_WW;

    TList fListAllInit;
};
