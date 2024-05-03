#pragma once

#include "TAxis.h"
class TF1;
class TF2;
class TH1;

class PtWeights
{
public:
    PtWeights(std::string const &_fID,
              TH1 &_hMCGen_dn_dptG_inv,
              TF1 &_fTrgtDist_dn_dptG_inv,
              TAxis &_axisPtG);
              
    std::string const &GetID() const { return id; }    
    double GetPtWeight(double *ptG, double *) const;
    TF1 &GetTF1();
    
private:
    std::string id;
    TH1 &hMCGen_dn_dptG_inv;
    TF1 &fTrgtDist_dn_dptG_inv;
    TAxis axisPtG;
    TF1 &fTF1;
};
