#pragma once

#include "TAxis.h"
#include "TList.h"
#include "dN_dptR.h"
#include "PtWeights.h"
#include <string>

class TLegend;


class MCEffi_wRes
{
public:
    MCEffi_wRes(std::string const &_id,
                TF1 &_fGenDist_dn_dptG,
                TF1 &_fEffi_dp_dptG,
                TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                TAxis &_axisPtR,
                PtWeights *_tPtWeights = nullptr);
                
    std::string const &GetID() const { return id; }

    // TF1 &MeasuredEffiTF1_NW(Color_t theLineColor = kBlue);
    // TF1 *MeasuredEffiTF1_WW(Color_t theLineColor = kRed);

    // integrate TF1s
    // TH1 &SampleMeasuredEffi_NW(Color_t theLineColor = kBlue);
    // TH1 *SampleMeasuredEffi_WW(Color_t theLineColor = kRed);

    // integrate TF2
    TH1 &SampleMeasuredEffi_NW_2(Color_t theLineColor = kBlue);
    TH1 *SampleMeasuredEffi_WW_2(Color_t theLineColor = kBlue);

    void PlotAll(TLegend *theLeg = nullptr);

    // public data members
    TF1 &fGenDist_dn_dptG; // the assumed distribution of generated particles in MC (will be flat for AS MCs )
    dN_dptR tdN_dptR_NW;

private:
    // defining data members
    std::string id;
    PtWeights *tPtWeights;
    dN_dptR *tdN_dptR_WW;
    TAxis axisPtR;

    // expressive data members
    TH1 &SampleTH1FromTF2(TF2 &theTF2_d_dptG_dptR, std::string name = "");

    TH1 &SampleMeasuredEffi_generic(std::string const &theSuffix,
                                    TF1  &theNumF,
                                    TF1  &theDenF,
                                    TAxis const &theAxis) const;


    TF1 *fGenDist_dn_dptG_WW;
    TH1 *hMCGen_NW;
    TH1 *hMCRec_NW;
    TH1 *hEffiRec_NW;

    TH1 *hMCGen_WW;
    TH1 *hMCRec_WW;
    TH1 *hEffiRec_WW;

    TList fListAllInit;
};
