#include "../include/dN_dptR.h"
#include "../include/dN_dptR_integrand.h"
#include "../include/PtWeights.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"

#include "TF1.h"
#include "TF2.h"

dN_dptR::dN_dptR(std::string const &_id,
                 TF1 &_fGen_dn_dptG,
                 TF1 &_fEffi_dp_dptG,
                 utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                 PtWeights *_tPtWeights /*= nullptr*/,
                 std::vector<TObject *> *theVAllDrawableObjects /*= nullptr*/)
    : id(_id),
      vAllDrawableObjects(
          theVAllDrawableObjects ? theVAllDrawableObjects : new std::vector<TObject *>()),
      tIntegrand(_id + "_integrand",
                 _fGen_dn_dptG,
                 _fEffi_dp_dptG,
                 _vFits_ptG_i_dp_dr_wAxis,
                 _tPtWeights),
      tPtWeights_opt(_tPtWeights),
      fGen_dn_dptG_WW_opt(tPtWeights_opt ? &utils_TF1::TF1Product(id + "_fGen_dn_dptG_WW_opt",
                                                                  GetGenDist_dn_dptG_NW_clone(),
                                                                  tPtWeights_opt->GetTF1())
                                         : nullptr)
{
    printf("dN_dptR::dN_dptR(): created instance %s %s ptWeights.\n",
           id.data(), tPtWeights_opt ? "with" : "without");
}

dN_dptR::~dN_dptR()
{
    printf("dN_dptR::~dN_dptR(): destructor for %s called.\n", id.data());
}

double dN_dptR::Evaluate_1D(double *ptR, double *)
{
    // calculate necessary ptG range.
    /* r = (ptR - ptG) / ptG
     * ptG = ptR / (r + 1)
     * dp/dr = 0 for r outside [-1., 1.] */
    tIntegrand.SetTF1ParameterPtR(*ptR);
    TF1 &lTF1 = tIntegrand.GetTF1Reference();
    double lIntMin = TMath::Max(lTF1.GetXmin(), 0.5 * *ptR);
    return lTF1.Integral(lIntMin, lTF1.GetXmax());
}
