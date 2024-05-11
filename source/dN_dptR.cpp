#include "dN_dptR.h"

#include "TF1.h"
#include "TF2.h"
#include "PtWeights.h"

dN_dptR::dN_dptR(std::string const &_id,
                 TF1 &_fGen_dn_dptG,
                 TF1 &_fEffi_dp_dptG,
                 TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                 PtWeights *_tPtWeights /*= nullptr*/)
    : id(_id),
      tIntegrand(_id + "_integrand",
                 _fGen_dn_dptG,
                 _fEffi_dp_dptG,
                 _vFits_ptG_i_dp_dr_wAxis,
                 _tPtWeights),
      tPtWeights(_tPtWeights)
{
    printf("dN_dptR::dN_dptR(): created instance %s %s ptWeights.\n",
           id.data(), tPtWeights ? "with" : "without");
}

double dN_dptR::Evaluate_1D(double *ptR, double *)
{
    TF1 &lTF1 = tIntegrand.GetTF1();
    lTF1.SetParameter(0, *ptR);

    // calculate necessary ptG range.
    /* r = (ptR - ptG) / ptG
     * ptG = ptR / (r + 1)
     * dp/dr = 0 for r outside [-1., 1.] */
    double lIntMin = TMath::Max(lTF1.GetXmin(), 0.5 * *ptR);
    return lTF1.Integral(lIntMin, lTF1.GetXmax());
}

TF1 &dN_dptR::GetNewTF1(double ptR) const
{
    return tIntegrand.GetNewTF1(ptR);
}

TF2 &dN_dptR::GetNewTF2() const
{
    return tIntegrand.GetNewTF2();
}

dN_dptR_integrand  &dN_dptR::GetIntegrand() 
{
    return tIntegrand;
}

dN_dptR::~dN_dptR()
{
    printf("dN_dptR::~dN_dptR(): destructor for %s called.\n", id.data());
}