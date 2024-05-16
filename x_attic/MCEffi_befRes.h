#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2023.h"

class MyWeightedDistFunctionObject
{
public:
    MyWeightedDistFunctionObject(TF1 *_dataF, TH1 *_mcH, TF1 *theDistToBeWeighted) : fDataDist(_dataF),
                                                                                     fSampledMCHist(_mcH),
                                                                                     fDist(theDistToBeWeighted) {}
    double operator()(double *x, double *)
    {
        return fDist->Eval(*x) * fDataDist->Eval(*x) / fSampledMCHist->Interpolate(*x);
        // try bin center instead:
        //~ return fDist->Eval(*x) * fDataDist->Eval(*x) / fSampledMCHist->GetBinContent(fSampledMCHist->FindBin(*x));
    }
    TF1 *fDataDist;
    TH1 *fSampledMCHist;
    TF1 *fDist;
};
TF1 *buildWeightedDist(std::string id, TF1 *fD, TH1 *hSampledMC, TF1 *theDistToBeWeighted)
{
    MyWeightedDistFunctionObject *fobj = new MyWeightedDistFunctionObject(fD, hSampledMC, theDistToBeWeighted);
    return new TF1(id.data(), *fobj, theDistToBeWeighted->GetXmin(), theDistToBeWeighted->GetXmax(), 0); // create TF1 class.
}

/* purpose: toy sim for pt weights mechanism in neutral meson reconstruction
 * constructor: MCEffi(std::string _id, TH1* _hGen_i, TF1* _fD_i, TF1* _trueEffi) [_i means pt-invariant]
 *  _hGen_i: a histo of the MC generated spectrum in pt_inv form.
 *           it will only be used to obtain  fGenDist_i via fit
 *  _fD_i:   the true distribution of produced pions in data // plays the role of the target distribution shape
 * _trueEffi: the true reconstruction efficiency

 * */
class MCEffi
{
public:
    MCEffi(std::string _id, TH1 *_hGen_i, TF1 *_fD_i, TF1 *_trueEffi) : fID(_id),
                                                                        fListAllInit(new TList),
                                                                        hGenHistInit_i(_hGen_i),
                                                                        fDataDist_i(_fD_i),
                                                                        fTrueEffi(_trueEffi)

    {
        fDataDist_i->SetLineColor(kGreen);
        fDataDist_i->SetName((fID + "fDataDist_i").data());

        fGenDist_i = (TF1 *)fDataDist_i->Clone(Form("%s_clone", fDataDist_i->GetName())); // initialize fGenDist_i
        fGenDist_i->SetLineColor(kBlue);
        fGenDist_i->SetName((fID + "fGenDist_i").data());
        hGenHistInit_i->Fit(fGenDist_i, "N"); // fit it

        fGenDist = multiplyTF1ByX(*fGenDist_i, (fID + "fGenDist").data());
        fDataDist = multiplyTF1ByX(*fDataDist_i, (fID + "fDataDist").data());
        fDataDist->SetLineColor(kGreen);

        hMCGenSampled = getSampledH(*fGenDist, *hGenHistInit_i /*for binning*/);
        hMCGenSampled_i = divideTH1ByBinCenters(*hMCGenSampled);

        fGenDistWW = buildWeightedDist(fID + "WeightedGenDist", fDataDist_i, hMCGenSampled_i, fGenDist);
        fGenDistWW->SetLineColor(kBlack);

        fMCRecDist = getTF1Product(fID + "fMCRecDist", fTrueEffi, fGenDist);

        fMCRecDistWW = getTF1Product(fID + "fMCRecDistWW", fTrueEffi, fGenDistWW);

        hMCRec = getSampledH(*fMCRecDist, *hGenHistInit_i /*for binning*/);
        hMCGenSampledWW = getSampledH(*fGenDistWW, *hGenHistInit_i /*for binning*/);
        hMCRecWW = getSampledH(*fMCRecDistWW, *hGenHistInit_i /*for binning*/);

        AddToPlotList();
    }
    void AddToPlotList()
    {
        //~ fListAllInit->Add(fDataDist_i);
        fListAllInit->Add(fDataDist);
        //~ fListAllInit->Add(hGenHistInit_i);
        fListAllInit->Add(fGenDist);
        //~ fListAllInit->Add(fGenDist_i);
        fListAllInit->Add(hMCGenSampled);
        fListAllInit->Add(fGenDistWW);

        fListAllInit->Add(fTrueEffi);
        fListAllInit->Add(fMCRecDist);
        fListAllInit->Add(fMCRecDistWW);
    }

    void PlotAll(TLegend *theLeg = nullptr)
    {
        for (TObject *obj : *fListAllInit)
        {
            std::cout << obj->GetName() << std::endl;
            obj->Draw("same");
            if (theLeg)
            {
                theLeg->AddEntry(obj, obj->GetName());
            }
        }
    }

    TH1 *GetMeasuredEffiWoW() { return divideTH1ByTH1(*hMCRec, *hMCGenSampled, "", (fID + "hMeasuredEffi").data()); }
    TH1 *GetMeasuredEffiWW() { return divideTH1ByTH1(*hMCRecWW, *hMCGenSampledWW, "", (fID + "hMeasuredEffiWW").data()); }

    std::string fID;
    TH1 *hGenHistInit_i; // used to fit fGenDist
    TF1 *fGenDist;
    TF1 *fGenDist_i;

    TF1 *fDataDist;
    TF1 *fDataDist_i; // set to our latest fit to data

    TF1 *fGenDistWW;

    // reconstructed stuff starts below

    TF1 *fTrueEffi;

    TF1 *fMCRecDist;
    TF1 *fMCRecDistWW;

    TH1 *hMCGenSampled;   // sampled
    TH1 *hMCGenSampled_i; // calc
    TH1 *hMCRec;          // sampled

    TH1 *hMCGenSampledWW; // sampled
    TH1 *hMCRecWW;        // sampled
    TList *fListAllInit;
};
