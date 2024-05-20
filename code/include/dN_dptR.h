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
            PtWeights *_fPtWeights = nullptr,
            std::vector<TObject *> *theVAllDrawableObjects = nullptr);

    ~dN_dptR();

    // ptG will be integrated out
    // can't be const since tIntegrand.fTF1.par[0] will be set to *ptR
    double Evaluate_1D(double *ptR, double *);

    TF1 const &GetGenDist_dn_dptG_NW() const { return GetIntegrand().GetGenDist_dn_dptG_NW(); }
    TF1 &GetGenDist_dn_dptG_NW_clone() const { return GetIntegrand().GetGenDist_dn_dptG_NW_clone(); }

    TF1 const *GetGenDist_dn_dptG_WW_opt() const { return fGen_dn_dptG_WW_opt; }
    TF1 *GetGenDist_dn_dptG_WW_opt_clone() const
    {
        return fGen_dn_dptG_WW_opt ? (TF1 *)fGen_dn_dptG_WW_opt->Clone(Form("%s_clone",
                                                                            fGen_dn_dptG_WW_opt->GetName()))
                                   : (TF1 *)nullptr;
    }

    std::string const &GetID() const { return id; }
    dN_dptR_integrand const &GetIntegrand() const { return tIntegrand; }

    TF1 const &GetTF1WithLastSetPtR() { return tIntegrand.GetTF1WithLastSetPtR(); }
    TF2 const &GetTF2_dN_dptG_dptR() const { return tIntegrand.GetTF2_dN_dptG_dptR(); }

    // TF2 &GetTF2_dN_dptG_dptR_clone() const { return tIntegrand.GetTF2_dN_dptG_dptR(); }

private:
    std::string id;
    std::vector<TObject *> *vAllDrawableObjects;
    dN_dptR_integrand tIntegrand;
    PtWeights *tPtWeights_opt;
    TF1 *fGen_dn_dptG_WW_opt;
};
