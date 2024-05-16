#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2023.h"

class Integrand_dN_dptR__FO
{
public:
    Integrand_dN_dptR__FO(std::string const &_id,
                          TF1 &_fGen_dn_dptG,
                          TF1 &_fEffi_dp_dptG,
                          std::vector<TF1 *> &_vFits_ptG_i_dp_dr,
                          TAxis &_axisPtG) : id(_id),
                                             fGen_dn_dptG(_fGen_dn_dptG),
                                             fEffi_dp_dptG(_fEffi_dp_dptG),
                                             vFits_ptG_i_dp_dr(_vFits_ptG_i_dp_dr),
                                             axisPtG(_axisPtG)
    {
        printf("Integrand_dN_dptR__FO::Integrand_dN_dptR__FO(): created instance %s. vFits_ptG_i_dp_dr.size() = %zu\n", id.data(), vFits_ptG_i_dp_dr.size());
    }

    ~Integrand_dN_dptR__FO() { printf("Integrand_dN_dptR__FO::~Integrand_dN_dptR__FO(): destructor for %s called.\n", id.data()); }

    // x = ptG, p = ptR
    double operator()(double *x, double *p) const
    {

        double result = 0.;

        double ptG = *x;
        double ptR = *p;
        double r = ptR - ptG;
        if (r > -1. && r < 1.)
        {

            size_t iPtGBin = axisPtG.FindBin(ptG);
            bool inRange = iPtGBin < vFits_ptG_i_dp_dr.size();

            auto f_dp_dr = inRange ? vFits_ptG_i_dp_dr.at(iPtGBin) : static_cast<TF1 *>(nullptr);

            if (f_dp_dr)
            {
                printf(" p ");
            }
            double p_GtoR = f_dp_dr ? f_dp_dr->Eval(r) : 0.;
            double result = fGen_dn_dptG.Eval(ptG) * fEffi_dp_dptG.Eval(ptG) * p_GtoR;

            if (result)
            {
                printf("Integrand_dN_dptR__FO::operator(): ptG = %.1f ptR = %.1f r = %.1f i = %zu inRange %d p_GtoR = %f res = %f\n", ptG, ptR, r, iPtGBin, inRange, p_GtoR, result);
            }
        }

        return result;
    }
    TF1 *getTF1() const
    {
        return new TF1(id.data(), this, axisPtG.GetXmin(), axisPtG.GetXmax(), 1);
    }

private:
    std::string id;
    TF1 &fGen_dn_dptG;
    TF1 &fEffi_dp_dptG;
    std::vector<TF1 *> &vFits_ptG_i_dp_dr;
    TAxis &axisPtG;
};

class dN_dptR__FO
{
public:
    dN_dptR__FO(std::string const &_id,
                TF1 &_fGen_dn_dptG,
                TF1 &_fEffi_dp_dptG,
                std::vector<TF1 *> &_vFits_ptG_i_dp_dr,
                TAxis &_axisPtG) : id(_id),
                                   fIntegrand__FO(_id + "_integrand",
                                                  _fGen_dn_dptG,
                                                  _fEffi_dp_dptG,
                                                  _vFits_ptG_i_dp_dr,
                                                  _axisPtG),
                                   fIntegrandTF1(*fIntegrand__FO.getTF1())
    {
        printf("dN_dptR__FO::dN_dptR__FO():created instance %s. _vFits_ptG_i_dp_dr.size() = %zu\n", id.data(), _vFits_ptG_i_dp_dr.size());
    }
    ~dN_dptR__FO() { printf("dN_dptR__FO::~dN_dptR__FO(): destructor for %s called.\n", id.data()); }

    // x = ptR
    double operator()(double *x, double *)
    {
        fIntegrandTF1.SetParameter(0, *x);
        return fIntegrandTF1.Integral(fIntegrandTF1.GetXmin(), fIntegrandTF1.GetXmax());
    }
    TF1 *getTF1()
    {
        return new TF1(id.data(), this, fIntegrandTF1.GetXmin(), fIntegrandTF1.GetXmax(), 0);
    }
    Integrand_dN_dptR__FO const &
    getIntegrand() const { return fIntegrand__FO; }

private:
    std::string id;
    Integrand_dN_dptR__FO const fIntegrand__FO;
    TF1 &fIntegrandTF1;
};

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

        fGenDistWW = buildWeightedDist(fID + "WeightedGenDist",
                                       fDataDist_i,
                                       hMCGenSampled_i,
                                       fGenDist);
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
