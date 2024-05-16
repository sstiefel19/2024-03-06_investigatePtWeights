#include "../include/PtWeights.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"

#include "TF1.h"
#include "TF2.h"
#include "TH1.h"

PtWeights::PtWeights(std::string const &_id,
                     bool _bComputeInInvariantForm,
                     TH1 const &_hMCGen_dn_dptG,
                     TF1 const &_fTrgtDist_dn_dptG,
                     TAxis const &_axisPtG,
                     std::vector<TObject *> *theVAllDrawableObjects /*= nullptr*/)
    : id(_id),
      vAllDrawableObjects(
          theVAllDrawableObjects ? theVAllDrawableObjects : new std::vector<TObject *>()),
      bComputeInInvariantForm(_bComputeInInvariantForm),
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

    vAllDrawableObjects->insert(vAllDrawableObjects->end(),
                                {&hMCGen_dn_dptG, &fTrgtDist_dn_dptG, &fTF1});
    printf("PtWeights::PtWeights(): created instance:\n"
           "\tid: %s\n"
           "\tbComputeInInvariantForm: %d\n"
           "\thMCGen_dn_dptG: %s\n"
           "\tfTrgtDist_dn_dptG: %s\n"
           "\taxisPtG: nBins %d, xmin %f, xmax %f\n",
           id.data(), bComputeInInvariantForm,
           hMCGen_dn_dptG.GetName(), fTrgtDist_dn_dptG.GetName(),
           axisPtG.GetNbins(), axisPtG.GetXmin(), axisPtG.GetXmax());
}

PtWeights::PtWeights(std::string const &_fID)
    : id(_fID + "_dummyConstructor"),
      bComputeInInvariantForm(false),
      hMCGen_dn_dptG(*new TH1D((id + "_hMCGen_dn_dptG").data(), "dummy", 1, 0., 1.)),
      fTrgtDist_dn_dptG(*new TF1((id + "_fTrgtDist_dn_dptG").data(), "1.", 0., 1.)),
      axisPtG(1, 0., 1.),
      fTF1((*new TF1((id + "fTF1").data(), "1.", 0., 1.)))
{
    printf("PtWeights::PtWeights(): created dummy instance:\n"
           "\tid: %s\n"
           "\tbComputeInInvariantForm: %d\n"
           "\thMCGen_dn_dptG: %s\n"
           "\tfTrgtDist_dn_dptG: %s\n"
           "\taxisPtG: nBins %d, xmin %f, xmax %f\n",
           id.data(), bComputeInInvariantForm,
           hMCGen_dn_dptG.GetName(), fTrgtDist_dn_dptG.GetName(),
           axisPtG.GetNbins(), axisPtG.GetXmin(), axisPtG.GetXmax());
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
