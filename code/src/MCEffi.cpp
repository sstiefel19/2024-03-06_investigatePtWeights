#include "../include/MCEffi.h"

#include "/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_plotting.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_computational.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TH1.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"

#include <list>

// MCEffi
// public:
/*  1. std::string sID;
    2. TF1 &fMCGenDist_dn_dptG; // the natural distribution of generated particles in this MC.
    TAxis axisPtR;
    // defining since it holds the pt-weights instance
    3. dN_dptR *tdN_dptR_WW_opt; // to simulate the numbers of reconstructed particles with pt-weights

// expressive members
    4. dN_dptR tdN_dptR_NW; // to simulate the numbers of reconstructed particles without pt-weights
    5. TH1 *hMCGen_NW;
    6. TH1 *hMCRec_NW;
    8. TH1 *hMCGen_WW;             // optional
    9. TH1 *hMCRec_WW;             // optional
    10. TH1 *hMeasuredEffi_NW_2;
    11. TH1 *hMeasuredEffi_WW_2_opt; // optional
    12. TList fListAllInit;*/
MCEffi::MCEffi(std::string const &_id,
               TF1 &_fGenDist_dn_dptG,
               TF1 &_fTrueEffiOverall_dp_dptG,
               utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
               TAxis &_axisPtR,
               PtWeights *_tPtWeights /*= nullptr*/,
               std::vector<TObject *> *theVAllDrawableObjects /*= nullptr*/)
    : sID(_id),
      vAllDrawableObjects(
          theVAllDrawableObjects ? theVAllDrawableObjects : new std::vector<TObject *>()),
      axisPtR(_axisPtR),
      tdN_dptR_NW(_id + "_fdN_dptR",
                  _fGenDist_dn_dptG,
                  _fTrueEffiOverall_dp_dptG,
                  _vFits_ptG_i_dp_dr_wAxis,
                  nullptr /*tPtWeights*/,
                  vAllDrawableObjects),
      tdN_dptR_WW_opt(!_tPtWeights ? nullptr
                                   : new dN_dptR(_id + "_fdN_dptR",
                                                 _fGenDist_dn_dptG,
                                                 _fTrueEffiOverall_dp_dptG,
                                                 _vFits_ptG_i_dp_dr_wAxis,
                                                 _tPtWeights,
                                                 vAllDrawableObjects)),
      hMeasuredEffi_NW_2(nullptr),
      hMeasuredEffi_WW_2_opt(nullptr)
{
    printf("MCEffi::MCEffi(): created instance:\n"
           "\tid: %s\n"
           "\taxisPtR: nBins: %d, xmin: %f, xmax: %f\n"
           "\ttdN_dptR_NW: %s\n"
           "\ttdN_dptR_WW_opt: %s\n",
           sID.data(),
           axisPtR.GetNbins(),
           axisPtR.GetXmin(),
           axisPtR.GetXmax(),
           tdN_dptR_NW.GetID().data(),
           tdN_dptR_WW_opt ? tdN_dptR_WW_opt->GetID().data() : "nullptr");
}

TH1 const &MCEffi::GetMeasuredEffi_NW()
{
    if (hMeasuredEffi_NW_2)
    {
        return *hMeasuredEffi_NW_2;
    }
    SampleMeasuredEffi_NW_2();
    return *hMeasuredEffi_NW_2;
}

TH1 const *MCEffi::GetMeasuredEffi_WW()
{
    if (!tdN_dptR_WW_opt)
    {
        printf("MCEffi::GetMeasuredEffi_WW(): WARNING: Called for instance '%s' without pt-weights. Returning nullptr.\n",
               sID.data());
        return nullptr;
    }
    if (hMeasuredEffi_WW_2_opt)
    {
        return hMeasuredEffi_WW_2_opt;
    }
    SampleMeasuredEffi_WW_2();
    return hMeasuredEffi_WW_2_opt;
}

TH1 &MCEffi::GetMeasuredEffi_NW_clone()
{
    if (hMeasuredEffi_NW_2)
    {
        return *utils_files_strings::CloneTH1(*hMeasuredEffi_NW_2);
    }
    return *utils_files_strings::CloneTH1(SampleMeasuredEffi_NW_2());
}

TH1 *MCEffi::GetMeasuredEffi_WW_clone()
{
    if (hMeasuredEffi_WW_2_opt)
    {
        return utils_files_strings::CloneTH1(*hMeasuredEffi_WW_2_opt);
    }
    TH1 *h = SampleMeasuredEffi_WW_2();
    return h ? utils_files_strings::CloneTH1(*h) : nullptr;
}

TObject *MCEffi::GetObservableObject(std::string const &theObservableName)
{
    {
        if (theObservableName == "MeasuredEffi_NW")
        {
            return &GetMeasuredEffi_NW_clone();
        }
        else if (theObservableName == "MeasuredEffi_WW")
        {
            return GetMeasuredEffi_WW_clone();
        }
        else
        {
            printf("InvPtW::CompareObservables_generic():\n\t"
                   "ERROR: unknown observable: %s\n"
                   "\tReturning dummy TH1.\n",
                   theObservableName.data());
            return new TH1D();
        }
    };
}

// ================== private member functions =======================
TH1 &MCEffi::SampleMeasuredEffi_generic_1D(std::string const &theName,
                                           TF1 const &theNumF,
                                           TF1 const &theDenF,
                                           TAxis const &theAxis) const

{
    TH1 &hNum = *utils_computational::GetSampledH(theNumF, theAxis);
    TH1 &hDen = *utils_computational::GetSampledH(theDenF, theAxis);
    TH1 &hRatio = *utils_TH1::DivideTH1ByTH1(hNum, hDen, "", theName.data(), theName.data());
    delete &hNum;
    delete &hDen;
    return hRatio;
}

TH1 &MCEffi::SampleMeasuredEffi_generic_2D(std::string const &theName,
                                           TF2 const &theNumTF2_dN_dptG_dptR,
                                           TF1 const &theDenTF1_dN_dptR) const
{
    TH1 &hNum = utils_computational::SampleTH1FromTF2_projectionY(theName + "_num",
                                                                  theNumTF2_dN_dptG_dptR,
                                                                  axisPtR);

    TH1 &hDen = *utils_computational::GetSampledH(theDenTF1_dN_dptR, axisPtR);
    TH1 &hRatio = *utils_TH1::DivideTH1ByTH1(hNum, hDen, "", theName.data(), theName.data());
    delete &hNum;
    delete &hDen;
    return hRatio;
}

TH1 &MCEffi::SampleMeasuredEffi_NW_2()
{
    hMeasuredEffi_NW_2 = &SampleMeasuredEffi_generic_2D(sID + "_hMeasuredEffi_NW_2",
                                                        tdN_dptR_NW.GetTF2_dN_dptG_dptR(),
                                                        tdN_dptR_NW.GetGenDist_dn_dptG());
    hMeasuredEffi_NW_2->SetLineColor(kPink);
    hMeasuredEffi_NW_2->SetMarkerColor(kPink);
    vAllDrawableObjects->push_back(hMeasuredEffi_NW_2);
    return *hMeasuredEffi_NW_2;
}
TH1 *MCEffi::SampleMeasuredEffi_WW_2() /*vAllDrawableObjects->push_back()*/
{
    if (!tdN_dptR_WW_opt)
    {
        printf("MCEffi::SampleMeasuredEffi_MC_WW_2(): Called for instance '%s' without pt-weights. returning nullptr.\n",
               sID.data());
        return nullptr;
    }

    hMeasuredEffi_WW_2_opt = &SampleMeasuredEffi_generic_2D(sID + "_hMeasuredEffi_WW_2",
                                                            tdN_dptR_WW_opt->GetTF2_dN_dptG_dptR(),
                                                            tdN_dptR_WW_opt->GetGenDist_dn_dptG());
    hMeasuredEffi_WW_2_opt->SetLineColor(kRed);
    hMeasuredEffi_WW_2_opt->SetMarkerColor(kRed);
    vAllDrawableObjects->push_back(hMeasuredEffi_WW_2_opt);
    return hMeasuredEffi_WW_2_opt;
}

void MCEffi::PlotAll(TLegend *theLeg /*= new TLegend(.73, .64, .90, .90, "")*/)
{
    TCanvas *c0 = new TCanvas((sID + "_PlotAll").data(), (sID + "_PlotAll").data(), 2000, 1000);
    TH2 *hd20 = new TH2F((sID + "_PlotAll_hd").data(), ";pT (GeV);various", 1, 0., 10., 1, 1e-5, 1e6);
    hd20->Draw();
    gPad->SetLogy();

    // for (auto const &i : data)

    for (TObject const *obj : *vAllDrawableObjects)
    {
        // make copy for drawing (which is not const)
        TObject &lObj = *obj->Clone(Form("%s_clone", obj->GetName()));
        printf("MCEffi::PlotAll(): instance %s: drawing object: %s\n",
               sID.data(),
               lObj.GetName());
        lObj.Draw("same");
        if (theLeg)
        {
            theLeg->AddEntry(obj, obj->GetName());
        }
    }
    if (theLeg)
    {
        theLeg->Draw("same");
    }
    utils_plotting::SaveCanvasAs(*c0);
}
