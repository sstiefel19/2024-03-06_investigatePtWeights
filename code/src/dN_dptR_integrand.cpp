#include "../include/dN_dptR_integrand.h"
#include "../include/PtWeights.h"

#include "TF1.h"
#include "TF2.h"

dN_dptR_integrand::dN_dptR_integrand(std::string const &_id,
                                     TF1 &_fGen_dn_dptG,
                                     TF1 &_fEffi_dp_dptG,
                                     utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                                     PtWeights *_tPtWeights /*= nullptr*/)
    : id(_id),
      fGen_dn_dptG(_fGen_dn_dptG),
      fEffi_dp_dptG(_fEffi_dp_dptG),
      vFits_ptG_i_dp_dr(_vFits_ptG_i_dp_dr_wAxis.first),
      axisPtG(_vFits_ptG_i_dp_dr_wAxis.second),
      tPtWeights(_tPtWeights),
      fTF1(nullptr),
      fTF2(nullptr)
{
    fTF1 = &GetNewTF1(0.);
    fTF2 = &GetNewTF2();

    printf("dN_dptR_integrand::dN_dptR_integrand(): created instance %s %s ptWeights.\n",
           id.data(), tPtWeights ? "with" : "without");
    printf("dN_dptR_integrand::dN_dptR_integrand(): id: %s, axisPtG: xmin: %f, xmax: %f, nBins: %d\n",
           id.data(), axisPtG.GetXmin(), axisPtG.GetXmax(), axisPtG.GetNbins());
}

dN_dptR_integrand::~dN_dptR_integrand()
{
    printf("dN_dptR_integrand::~dN_dptR_integrand(): destructor for %s called.\n", id.data());
}

// x = ptG, pars[0] = ptR. needed for fTF1
double dN_dptR_integrand::Evaluate_1D(double *ptG, double *ptR) const
{
    double result = 0.;

    double r = *ptG ? (*ptR - *ptG) / *ptG : 0.;

    if (r > -1. && r < 1.)
    {
        size_t iPtGBin = axisPtG.FindBin(*ptG);
        bool inRange = iPtGBin < vFits_ptG_i_dp_dr.size();

        auto f_dp_dr = inRange ? vFits_ptG_i_dp_dr.at(iPtGBin) : static_cast<TF1 *>(nullptr);

        double p_GtoR = f_dp_dr ? f_dp_dr->Eval(r) : 0.;
        double ptWeight = !p_GtoR ? 0. : tPtWeights ? tPtWeights->GetPtWeight(ptG, nullptr)
                                                    : 1.;
        result = fGen_dn_dptG.Eval(*ptG) * fEffi_dp_dptG.Eval(*ptG) * p_GtoR * ptWeight;

        if (result)
        {
            //~ printf("dN_dptR_integrand::Evaluate(): ptG = %.1f ptR = %.1f r = %.1f i = %zu inRange %d p_GtoR = %f result = %f\n", *ptG, ptR, r, iPtGBin, inRange, p_GtoR, result);
        }
    }
    return result;
}

// xy[0] = ptG, xy[1] = ptR. needed for fTF2
double dN_dptR_integrand::Evaluate_2D(double *xy, double *) const
{
    double result = 0.;

    double ptG = xy[0];
    double ptR = xy[1];

    double r = ptG ? (ptR - ptG) / ptG : 0.;

    if (r > -1. && r < 1.)
    {
        size_t iPtGBin = axisPtG.FindBin(ptG);
        bool inRange = iPtGBin < vFits_ptG_i_dp_dr.size();

        auto f_dp_dr = inRange ? vFits_ptG_i_dp_dr.at(iPtGBin) : static_cast<TF1 *>(nullptr);

        double p_GtoR = f_dp_dr ? f_dp_dr->Eval(r) : 0.;
        double ptWeight = !p_GtoR ? 0. : tPtWeights ? tPtWeights->GetPtWeight(&ptG, nullptr)
                                                    : 1.;
        result = fGen_dn_dptG.Eval(ptG) * fEffi_dp_dptG.Eval(ptG) * p_GtoR * ptWeight;

        if (result)
        {
            //~ printf("dN_dptR_integrand::Evaluate(): ptG = %.1f ptR = %.1f r = %.1f i = %zu inRange %d p_GtoR = %f result = %f\n", ptG, ptR, r, iPtGBin, inRange, p_GtoR, result);
        }
    }
    return result;
}

// creates a new TF1 instance in which ptR is a parameter
TF1 &dN_dptR_integrand::GetNewTF1(double ptR) const
{
    auto &f = *new TF1((id + "_TF1").data(), this,
                       &dN_dptR_integrand::Evaluate_1D,
                       axisPtG.GetXmin(), axisPtG.GetXmax(), 1,
                       "dN_dptR_integrand", "Evaluate");
    f.SetParameter(0, ptR);
    return f;
}

TF1 &dN_dptR_integrand::GetTF1()
{
    return fTF1 ? *fTF1 : GetNewTF1(0.);
}

// creates a new TF2 instance with ptG, ptR as variables
TF2 &dN_dptR_integrand::GetNewTF2() const
{
    auto &f = *new TF2((id + "_TF2").data(), this,
                       &dN_dptR_integrand::Evaluate_2D,
                       axisPtG.GetXmin(), axisPtG.GetXmax(),
                       axisPtG.GetXmin(), axisPtG.GetXmax(),
                       0,
                       "dN_dptR_integrand", "Evaluate_2D");
    return f;
}

TF2 &dN_dptR_integrand::GetTF2()
{
    return fTF2 ? *fTF2 : GetNewTF2();
}
