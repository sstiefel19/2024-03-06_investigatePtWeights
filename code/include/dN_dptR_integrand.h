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

    TF1 const &GetGenDist_dn_dptG() const { return fGen_dn_dptG; }
    TF1 &GetGenDist_dn_dptG_clone() const
    {
        return *(TF1 *)fGen_dn_dptG.Clone(Form("%s_clone", fGen_dn_dptG.GetName()));
    }
    std::string const &GetID() const { return id; }

    TF1 &GetTF1Reference() { return fTF1; }
    TF1 const &GetTF1WithLastSetPtR() const { return fTF1; }
    TF2 const &GetTF2_dN_dptG_dptR() const { return fTF2; }

    void SetTF1ParameterPtR(double ptR) { fTF1.SetParameter(0, ptR); }

private:
    // defining data members
    std::string id;
    std::vector<TObject *> *vAllDrawableObjects;
    TF1 &fGen_dn_dptG;                     // the distribution of generated particles
    TF1 &fEffi_dp_dptG;                    // prob that a meson generated with ptG gets reconstructed (no matter with which ptR)
    std::vector<TF1 *> &vFits_ptG_i_dp_dr; // for every bin i in ptG: prob dist that the meson gets smeared with r =  (ptR - ptG) / ptG . Normalized to 1
    TAxis &axisPtG;                        // the binning in ptG used in vFits_ptG_i_dp_dr
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
