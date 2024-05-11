#include "PtWeights.h"
#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2024.h"

#include "TF1.h"
#include "TF2.h"
#include "TH1.h"

PtWeights::PtWeights(std::string const &_id,
                     TH1 const &_hMCGen_dn_dptG,
                     TF1 const &_fTrgtDist_dn_dptG,
                     TAxis const &_axisPtG)
    : id(_id),
      hMCGen_dn_dptG(*cloneTH1(_hMCGen_dn_dptG,
                               "",
                               (id + "_hMCGen_dn_dptG").data())),
      fTrgtDist_dn_dptG(
          dynamic_cast<TF1 &>(*cloneTNamedObject(_fTrgtDist_dn_dptG,
                                                 "",
                                                 (id + "_fTrgtDist_dn_dptG").data()))),
      axisPtG(_axisPtG),
      fTF1(*new TF1((id + "_PtWeightTF1").data(),
                    this, &PtWeights::GetPtWeight,
                    axisPtG.GetXmin(), axisPtG.GetXmax(), 0, "PtWeights", "GetPtWeight"))
{
    fTrgtDist_dn_dptG.SetRange(axisPtG.GetXmin(), axisPtG.GetXmax());

    printf("PtWeights::PtWeights(): created instance %s\n", id.data());
}

double PtWeights::GetPtWeight(double *ptG, double *) const
{
    double denom = hMCGen_dn_dptG.Interpolate(*ptG);
    if (!denom)
    {
        printf("GetPtWeight(): weight is 0 for ptG = %f. Name instance: %s\n", *ptG, id.data());
    }
    return denom ? fTrgtDist_dn_dptG.Eval(*ptG) / denom : 0.;
}

TF1 &PtWeights::GetTF1()
{
    return fTF1;
}
