#pragma once

// for utils_fits::TPairFitsWAxis
#include "ComputeResolutionFits.h"

#include "TF1.h"
#include "TF2.h"

class TAxis;
class PtWeights;

class dN_dptR_integrand
{
public:
    dN_dptR_integrand(std::string const &_id,
                      TF1 &_fGen_dn_dptG,
                      TF1 &_fEffi_dp_dptG,
                      utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                      PtWeights *_tPtWeights = nullptr,
                      std::vector<TObject *> *theVAllDrawableObjects = nullptr);

    ~dN_dptR_integrand();

    TF1 const &GetGenDist_dn_dptG_NW() const;
    TF1 &GetGenDist_dn_dptG_NW_clone() const;
    std::string const &GetID() const;

    TF1 &GetTF1Reference();
    TF1 const &GetTF1WithLastSetPtR() const;
    TF2 const &GetTF2_dN_dptG_dptR() const;

    void SetTF1ParameterPtR(double ptR);

private:
    // defining data members
    std::string id;
    std::vector<TObject *> *vAllDrawableObjects;
    TF1 fGen_dn_dptG_NW;                   // the distribution of generated particles
    TF1 fEffi_dp_dptG;                     // prob that a meson generated with ptG gets reconstructed (no matter with which ptR)
    std::vector<TF1 *> &vFits_ptG_i_dp_dr; // for every bin i in ptG: prob dist that the meson gets smeared with r =  (ptR - ptG) / ptG . Normalized to 1
    TAxis const axisPtG;                   // the binning in ptG used in vFits_ptG_i_dp_dr
    PtWeights *tPtWeights_opt;

    // expressive data members
    // holds a TF1 in which ptG is a variable and ptR a parameter.
    // It is used for dN_dptR::Evaluate_1D in which ptG is integrated out
    TF1 &fTF1; // create on heap so the actual TF1 object can outlive *this* instance

    // holds a TF2 with variables ptG, ptR. It is used for dN_dptR::Evaluate_2D
    TF2 &fTF2; // create on heap so the actual TF1 object can outlive *this* instance

    // x = ptG, pars[0] = ptR. needed for fTF1
    double Evaluate_1D(double *ptG, double *ptR) const;

    // xy[0] = ptG, xy[1] = ptR. needed for fTF2
    double Evaluate_2D(double *xy, double *) const;

    // creates a new TF1 instance in which ptR is a parameter
    TF1 &GetNewTF1(double ptR) const;

    // creates a new TF2 instance with ptG, ptR as variables
    TF2 &GetNewTF2_dN_dptG_dptR() const;
};
