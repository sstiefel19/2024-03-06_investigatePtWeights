#include "../include/dN_dptR_integrand.h"
#include "../include/PtWeights.h"

dN_dptR_integrand::dN_dptR_integrand(std::string const &_id,
                                     TF1 &_fGen_dn_dptG,
                                     TF1 &_fEffi_dp_dptG,
                                     utils_fits::TPairFitsWAxis &_vFits_ptG_i_dp_dr_wAxis,
                                     PtWeights *_tPtWeights /*= nullptr*/,
                                     std::vector<TObject *> *theVAllDrawableObjects /*= nullptr*/)
    : id(_id),
      vAllDrawableObjects(
          theVAllDrawableObjects ? theVAllDrawableObjects
                                 : new std::vector<TObject *>()),
      fGen_dn_dptG_NW(_fGen_dn_dptG),
      fEffi_dp_dptG(_fEffi_dp_dptG),
      vFits_ptG_i_dp_dr(_vFits_ptG_i_dp_dr_wAxis.first),
      axisPtG(_vFits_ptG_i_dp_dr_wAxis.second),
      tPtWeights_opt(_tPtWeights),
      fTF1(GetNewTF1(0. /*ptR*/)),
      fTF2(GetNewTF2_dN_dptG_dptR())
{
    // correct names
    fGen_dn_dptG_NW.SetName((id + "_fGen_dn_dptG").data());
    fEffi_dp_dptG.SetName((id + "_fEffi_dp_dptG").data());

    vAllDrawableObjects->insert(
        vAllDrawableObjects->end(),
        std::initializer_list<TObject *>({&fGen_dn_dptG_NW,
                                          &fEffi_dp_dptG,
                                          &fTF1,
                                          &fTF2}));
    printf("dN_dptR_integrand::dN_dptR_integrand(): created instance %s %s ptWeights.\n",
           id.data(), tPtWeights_opt ? "with" : "without");
    printf("dN_dptR_integrand::dN_dptR_integrand(): id: %s, axisPtG: xmin: %f, xmax: %f, nBins: %d\n",
           id.data(), axisPtG.GetXmin(), axisPtG.GetXmax(), axisPtG.GetNbins());
}

dN_dptR_integrand::~dN_dptR_integrand()
{
    printf("dN_dptR_integrand::~dN_dptR_integrand(): destructor for %s called.\n", id.data());
}

// ======================== private member functions ========================
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
        double ptWeight = !p_GtoR ? 0. : tPtWeights_opt ? tPtWeights_opt->GetPtWeight(ptG, nullptr)
                                                        : 1.;
        result = fGen_dn_dptG_NW.Eval(*ptG) * fEffi_dp_dptG.Eval(*ptG) * p_GtoR * ptWeight;

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
        double ptWeight = !p_GtoR ? 0. : tPtWeights_opt ? tPtWeights_opt->GetPtWeight(&ptG, nullptr)
                                                        : 1.;
        result = fGen_dn_dptG_NW.Eval(ptG) * fEffi_dp_dptG.Eval(ptG) * p_GtoR * ptWeight;

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
    vAllDrawableObjects->push_back(&f);
    return f;
}

TF2 &dN_dptR_integrand::GetNewTF2_dN_dptG_dptR() const
{
    auto &f = *new TF2((id + "_TF2").data(), this,
                       &dN_dptR_integrand::Evaluate_2D,
                       axisPtG.GetXmin(), axisPtG.GetXmax(),
                       axisPtG.GetXmin(), axisPtG.GetXmax(),
                       0,
                       "dN_dptR_integrand", "Evaluate_2D");
    vAllDrawableObjects->push_back(&f);
    return f;
}

TF1 const &dN_dptR_integrand::GetGenDist_dn_dptG_NW() const { return fGen_dn_dptG_NW; }
    TF1 &dN_dptR_integrand::GetGenDist_dn_dptG_NW_clone() const
    {
        return *(TF1 *)fGen_dn_dptG_NW.Clone(Form("%s_clone", fGen_dn_dptG_NW.GetName()));
    }
    std::string const &dN_dptR_integrand::GetID() const { return id; }

    TF1 &dN_dptR_integrand::GetTF1Reference() { return fTF1; }
    TF1 const &dN_dptR_integrand::GetTF1WithLastSetPtR() const { return fTF1; }
    TF2 const &dN_dptR_integrand::GetTF2_dN_dptG_dptR() const { return fTF2; }

    void dN_dptR_integrand::SetTF1ParameterPtR(double ptR) { fTF1.SetParameter(0, ptR); }
