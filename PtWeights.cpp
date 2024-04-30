#include "PtWeights.h"

#include "TF1.h"
#include "TF2.h"
#include "TH1.h"

PtWeights::PtWeights(std::string const &_id,
                     TH1 &_hMCGen_dn_dptG_inv,
                     TF1 &_fTrgtDist_dn_dptG_inv,
                     TAxis &_axisPtG)
    : id(_id),
      hMCGen_dn_dptG_inv(_hMCGen_dn_dptG_inv),
      fTrgtDist_dn_dptG_inv(_fTrgtDist_dn_dptG_inv),
      axisPtG(_axisPtG),
      fTF1(*new TF1((id + "_PtWeightTF1").data(),
                    this, &PtWeights::GetPtWeight,
                    axisPtG.GetXmin(), axisPtG.GetXmax(), 0, "PtWeights", "GetPtWeight"))
{
    fTrgtDist_dn_dptG_inv.SetRange(axisPtG.GetXmin(), axisPtG.GetXmax());

    printf("PtWeights::PtWeights(): created instance %s\n", id.data());
}

double PtWeights::GetPtWeight(double *ptG, double *) const
{
    double denom = hMCGen_dn_dptG_inv.Interpolate(*ptG);
    if (!denom)
    {
        printf("GetPtWeight(): weight is 0 for ptG = %f. Name instance: %s\n", *ptG, id.data());
    }
    return denom ? fTrgtDist_dn_dptG_inv.Eval(*ptG) / denom : 0.;
}

TF1 &PtWeights::GetTF1()
{
    return fTF1;
}
