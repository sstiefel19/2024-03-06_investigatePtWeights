#include "../include/MCEffi.h"

#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_sstiefel_2024.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"

// MCEffi
// public:
MCEffi::MCEffi(std::string const &_id,
                         TF1 &_fGenDist_dn_dptG,
                         TF1 &_fEffi_dp_dptG,
                         utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                         TAxis &_axisPtR,
                         PtWeights *_tPtWeights /*= nullptr*/)
    : id(_id),
      fListAllInit(),
      fGenDist_dn_dptG(_fGenDist_dn_dptG),
      tdN_dptR_NW(_id + "_fdN_dptR",
                  _fGenDist_dn_dptG,
                  _fEffi_dp_dptG,
                  _vFits_ptG_i_dp_dr_wAxis),
      axisPtR(_axisPtR),
      tPtWeights(_tPtWeights),
      tdN_dptR_WW(!_tPtWeights ? nullptr
                               : new dN_dptR(_id + "_fdN_dptR",
                                             _fGenDist_dn_dptG,
                                             _fEffi_dp_dptG,
                                             _vFits_ptG_i_dp_dr_wAxis,
                                             _tPtWeights)),
      fGenDist_dn_dptG_WW(
          _tPtWeights ? &utils_TF1::TF1Product(id + "_fGenDist_dn_dptG_WW",
                                               tPtWeights->GetTF1(),
                                               fGenDist_dn_dptG)
                      : nullptr),
      hMCGen_NW(nullptr),
      hMCRec_NW(nullptr),
      hEffiRec_NW(nullptr),

      hMCGen_WW(nullptr),
      hMCRec_WW(nullptr),
      hEffiRec_WW(nullptr)
{
    id + "_fGenDist_dn_dptG_WW",
        fGenDist_dn_dptG.SetLineColor(kBlack);
    fListAllInit.Add(&fGenDist_dn_dptG);

    printf("MCEffi::MCEffi(): created instance %s %s ptWeights.\n",
           id.data(),
           tPtWeights ? "with" : "without");
    printf("MCEffi::MCEffi(): id: %s, axisPtR: xmin: %f, xmax: %f, nBins: %d\ntdN_dptR_NW: %s\n",
           id.data(),
           axisPtR.GetXmin(),
           axisPtR.GetXmax(),
           axisPtR.GetNbins(),
           tdN_dptR_NW.GetID().data());

    if (tPtWeights)
    {
        fListAllInit.Add(fGenDist_dn_dptG_WW);
        tPtWeights->GetTF1().SetLineColor(kBlue);
        fListAllInit.Add(&tPtWeights->GetTF1());
        printf("MCEffi::MCEffi(): id: %s, tdN_dptR_WW: %s, tPtWeights: %s\n",
               id.data(),
               tdN_dptR_WW->GetID().data(),
               tPtWeights->GetID().data());
    }
}

void MCEffi::PlotAll(TLegend *theLeg /*= new TLegend(.73, .64, .90, .90, "")*/)
{
    TCanvas *c0 = new TCanvas((id + "_PlotAll").data(), (id + "_PlotAll").data(), 2000, 1000);
    TH2 *hd20 = new TH2F((id + "_PlotAll_hd").data(), ";pT (GeV);various", 1, 0., 10., 1, 1e-5, 1e6);
    hd20->Draw();
    gPad->SetLogy();

    for (TObject *obj : fListAllInit)
    {
        cout << obj->GetName() << endl;
        obj->Draw("same");
        if (theLeg)
        {
            theLeg->AddEntry(obj, obj->GetName());
        }
    }
    if (theLeg)
    {
        theLeg->Draw("same");
    }
    saveCanvasAs(*c0);
}

// TF1 &MCEffi::MeasuredEffiTF1_NW(Color_t theLineColor /*= kBlue*/)
// {
//     auto &f = *getTF1Division(id + "_GetMeasuredEffiTF1_NW",
//                               &tdN_dptR_NW.GetTF1(),
//                               &fGenDist_dn_dptG);
//     f.SetLineColor(theLineColor);
//     fListAllInit.Add(&f);
//     return f;
// }

// TF1 *MCEffi::MeasuredEffiTF1_WW(Color_t theLineColor /*= kRed*/)
// {
//     if (!(tdN_dptR_WW && fGenDist_dn_dptG_WW))
//     {
//         printf("MCEffi::GetMeasuredEffiTF1_WW() called for instance %s without weights!. Returning nullptr.\n", id.data());
//         return nullptr;
//     }
//     auto *f = getTF1Division(id + "_GetMeasuredEffiTF1_WW",
//                              &tdN_dptR_WW->GetTF1(),
//                              fGenDist_dn_dptG_WW);

//     f->SetLineColor(theLineColor);
//     fListAllInit.Add(f);
//     return f;
// }

TH1 &MCEffi::SampleMeasuredEffi_generic(std::string const &theName,
                                             TF1 &theNumF,
                                             TF1 &theDenF,
                                             TAxis const &theAxis) const
{
    TH1 &hNum = *getSampledH(theNumF, theAxis);
    TH1 &hDen = *getSampledH(theDenF, theAxis);
    TH1 &hRatio = *divideTH1ByTH1(hNum, hDen, "", theName.data());
    delete &hNum;
    delete &hDen;
    return hRatio;
}

// TH1 &MCEffi::SampleMeasuredEffi_NW(Color_t theLineColor /*= kBlue*/)
// {
//     hEffiRec_NW = &SampleMeasuredEffi_generic(id + "_hEffiRec_NW",
//                                               tdN_dptR_NW.GetTF1(),
//                                               fGenDist_dn_dptG,
//                                               axisPtR);
//     hEffiRec_NW->SetLineColor(theLineColor);
//     fListAllInit.Add(hEffiRec_NW);
//     return *hEffiRec_NW;
// }

// TH1 *MCEffi::SampleMeasuredEffi_WW(Color_t theLineColor /*= kRed*/)
// {
//     if (!(tdN_dptR_WW && fGenDist_dn_dptG_WW))
//     {
//         printf("MCEffi::SampleMeasuredEffi_WW(): Called but no tPtWeights initialized. Instance %s. returning nullptr.\n", id.data());
//         return nullptr;
//     }
//     hMCGen_WW = getSampledH(*fGenDist_dn_dptG_WW, axisPtR);
//     hMCRec_WW = getSampledH(tdN_dptR_WW->GetTF1(), axisPtR);
//     hEffiRec_WW = divideTH1ByTH1(*hMCRec_WW, *hMCGen_WW, "", (id + "_hEffiRec_WW").data());
//     hEffiRec_WW->SetLineColor(theLineColor);
//     fListAllInit.Add(hEffiRec_WW);
//     return hEffiRec_WW;
// }

// integrate TF2
TH1 &MCEffi::SampleMeasuredEffi_NW_2(Color_t theLineColor /*= kBlue*/)
{

    if (!hMCGen_NW)
    {
        hMCGen_NW = getSampledH(fGenDist_dn_dptG, axisPtR);
    }
    TH1 &hMCRec_NW_2 = SampleTH1FromTF2(tdN_dptR_NW.GetIntegrand().GetTF2());
    // fListAllInit.Add(&hMCRec_NW_2);

    TH1 &hEffiRec_NW2 = *divideTH1ByTH1(hMCRec_NW_2, *hMCGen_NW, "", (id + "_SampleMeasuredEffi_NW_2").data());
    hEffiRec_NW2.SetLineColor(theLineColor);
    fListAllInit.Add(&hEffiRec_NW2);
    return hEffiRec_NW2;
}

TH1 *MCEffi::SampleMeasuredEffi_WW_2(Color_t theLineColor /*= kBlue*/)
{
    if (!(tdN_dptR_WW && fGenDist_dn_dptG_WW))
    {
        printf("MCEffi::SampleMeasuredEffi_WW_2(): Called but no tPtWeights initialized. Instance %s. returning nullptr.\n", id.data());
        return nullptr;
    }
    if (!hMCGen_WW)
    {
        hMCGen_WW = getSampledH(*fGenDist_dn_dptG_WW, axisPtR);
    }
    TH1 &hMCRec_WW_2 = SampleTH1FromTF2(tdN_dptR_WW->GetIntegrand().GetTF2());

    TH1 &hEffiRec_WW2 = *divideTH1ByTH1(hMCRec_WW_2, *hMCGen_WW, "", (id + "SampleMeasuredEffi_WW_2").data());
    hEffiRec_WW2.SetLineColor(theLineColor);
    fListAllInit.Add(&hEffiRec_WW2);
    return &hEffiRec_WW2;
}

// private:
TH1 &MCEffi::SampleTH1FromTF2(TF2 &theTF2_d_dptG_dptR, std::string name /*=""*/)
{
    if (!name.size())
    {
        name = id + "_SampleTH1FromTF2_" + theTF2_d_dptG_dptR.GetName();
    }
    TH1 &h1_d_dptR = *new TH1D(name.data(),
                               name.data(),
                               axisPtR.GetNbins(),
                               axisPtR.GetXmin(),
                               axisPtR.GetXmax());

    for (int i = 1; i <= h1_d_dptR.GetNbinsX(); ++i)
    {
        double integral = theTF2_d_dptG_dptR.Integral(
            theTF2_d_dptG_dptR.GetXmin(), // ptG axis min
            theTF2_d_dptG_dptR.GetXmax(), // ptG axis max
            axisPtR.GetBinLowEdge(i),
            axisPtR.GetBinUpEdge(i));
        double width = axisPtR.GetBinWidth(i);
        h1_d_dptR.SetBinContent(i, width ? integral / width : 0.);
    }
    // fListAllInit.Add(&h1_d_dptR);
    return h1_d_dptR;
}
