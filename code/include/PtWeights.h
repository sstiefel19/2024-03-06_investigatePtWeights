#pragma once

#include "TAxis.h"
class TF1;
class TF2;
class TH1;

class PtWeights
{
public:
    PtWeights(std::string const &_fID,
              bool _bComputeInInvariantForm,
              TH1 const &_hMCGen_dn_dptG,
              TF1 const &_fTrgtDist_dn_dptG,
              TAxis const &_axisPtG);

    TAxis const &GetAxisPtG() const { return axisPtG; }
    TF1 &GetTF1TrgtDist_dn_dptG() const { return fTrgtDist_dn_dptG; }
    TF1 const &GetTF1TrgtDist_dn_dptG_const() const { return fTrgtDist_dn_dptG; }
    std::string const &GetID() const { return id; }
    TH1 const &GetTH1MCGen_dn_dptG() const { return hMCGen_dn_dptG; }

    double GetPtWeight(double *ptG, double *) const;
    TF1 &GetTF1();

private:
    std::string id;
    bool bComputeInInvariantForm;

    // create on heap so I can work with TH1 plus I don need to worry about lifetimes
    TH1 &hMCGen_dn_dptG;
    TF1 &fTrgtDist_dn_dptG;
    TAxis axisPtG;
    TF1 &fTF1;
};
