#pragma once

#include "../include/dN_dptR_integrand.h"

#include <string>

class TF1;
class TF2;
class PtWeights;

class dN_dptR
{
public:
    dN_dptR(std::string const &_id,
            TF1 &_fGen_dn_dptG,
            TF1 &_fEffi_dp_dptG,
            utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
            PtWeights *_fPtWeights = nullptr);

    // ptG will be integrated out
    // can't be const since tIntegrand.fTF1.par[0] will be set to *ptR
    double Evaluate_1D(double *ptR, double *);

    std::string const &GetID() const { return id; }    
    dN_dptR_integrand  &GetIntegrand();    
    TF1 &GetNewTF1(double ptR) const;
    TF2 &GetNewTF2() const;

    ~dN_dptR();

private:
    std::string id;
    dN_dptR_integrand tIntegrand;
    PtWeights *tPtWeights;
};
