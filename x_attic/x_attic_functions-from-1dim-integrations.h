TCanvas &
compareMeasuredEffis_TF1(std::string const &fID,
                         MCEffi_wRes &lMCEffi_AS,
                         MCEffi_wRes &lMCEffi_D)
{

    auto &fEffi_AS_NW = lMCEffi_AS.MeasuredEffiTF1_NW(kBlue);
    auto &fEffi_AS_WW = *lMCEffi_AS.MeasuredEffiTF1_WW(kGreen);
    auto &fEffi_D_NW = lMCEffi_D.MeasuredEffiTF1_NW(kRed);

    auto &fEffi_AS_NW_over_D = *getTF1Division(fID + "_fEffi_AS_NW_over_D", &fEffi_AS_NW, &fEffi_D_NW);
    auto &fEffi_AS_WW_over_D = *getTF1Division(fID + "_fEffi_AS_WW_over_D", &fEffi_AS_WW, &fEffi_D_NW);

    auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
    {
        o.Draw("same");
        if (leg)
            leg->AddEntry(&o, o.GetName(), "l");
    };
    // setup canvas
    TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TF1").data(),
                               (fID + "compareMeasuredEffis_TF1").data(), 2000, 1000);

    c1.Divide(1, 2);
    c1.cd(1);
    gPad->SetLogy();
    TH2 *hd21 = new TH2F("hd121", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 4e-3);
    hd21->Draw();

    auto leg = new TLegend(.73, .64, .90, .90, "");
    drawAndAdd(fEffi_AS_NW, leg);
    drawAndAdd(fEffi_AS_WW, leg);
    drawAndAdd(fEffi_D_NW, leg);
    leg->Draw("same");

    c1.cd(2);
    fEffi_AS_WW_over_D.Draw();
    // fEffi_AS_WW_over_D.Dump();

    //~ gPad->SetLogy();
    TH2 *hd22 = new TH2F("hd122", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
    hd22->Draw();

    auto leg2 = new TLegend();
    drawAndAdd(fEffi_AS_NW_over_D);
    drawAndAdd(fEffi_AS_WW_over_D);
    leg2->Draw("same");

    saveCanvasAs(c1);
    return c1;
}

TCanvas &
compareMeasuredEffis_TH1(std::string const &fID,
                         MCEffi_wRes &lMCEffi_AS,
                         MCEffi_wRes &lMCEffi_D)
{

    auto &hEffi_AS_NW = lMCEffi_AS.SampleMeasuredEffi_NW(kBlue);
    auto &hEffi_AS_WW = *lMCEffi_AS.SampleMeasuredEffi_WW(kGreen);
    auto &hEffi_D_NW = lMCEffi_D.SampleMeasuredEffi_NW(kRed);

    auto &hEffi_AS_NW_over_D = *divideTH1ByTH1(hEffi_AS_NW, hEffi_D_NW, "", "hEffi_AS_NW_over_D");
    auto &hEffi_AS_WW_over_D = *divideTH1ByTH1(hEffi_AS_WW, hEffi_D_NW, "", "hEffi_AS_WW_over_D");

    auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
    {
        o.Draw("same");
        if (leg)
            leg->AddEntry(&o, o.GetName(), "l");
    };
    // setup canvas
    TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TH1").data(),
                               (fID + "compareMeasuredEffis_TH1").data(), 2000, 1000);

    c1.Divide(1, 2);
    c1.cd(1);
    gPad->SetLogy();
    TH2 *hd21 = new TH2F("hd21", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 8e-3);
    hd21->Draw();

    auto leg = new TLegend(.73, .64, .90, .90, "");
    drawAndAdd(hEffi_AS_NW, leg);
    drawAndAdd(hEffi_AS_WW, leg);
    drawAndAdd(hEffi_D_NW, leg);
    leg->Draw("same");

    c1.cd(2);
    //~ gPad->SetLogy();
    TH2 *hd22 = new TH2F("hd22", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
    hd22->Draw();

    auto leg2 = new TLegend();
    drawAndAdd(hEffi_AS_NW_over_D);
    drawAndAdd(hEffi_AS_WW_over_D);
    //~ leg2->Draw("same");

    saveCanvasAs(c1);
    return c1;
}

TCanvas &
compareMeasuredEffis_TH1_New(std::string const &fID,
                             MCEffi_wRes &lMCEffi_AS,
                             MCEffi_wRes &lMCEffi_D)
{

    auto &hEffi_AS_NW = lMCEffi_AS.SampleMeasuredEffi_NW_2(kPink);
    auto &hEffi_AS_WW = *lMCEffi_AS.SampleMeasuredEffi_WW_2(kCyan);
    auto &hEffi_D_NW = lMCEffi_D.SampleMeasuredEffi_NW_2(kMagenta);

    auto &hEffi_AS_NW_over_D = *divideTH1ByTH1(hEffi_AS_NW, hEffi_D_NW, "", "hEffi_AS_NW_over_D_2");
    auto &hEffi_AS_WW_over_D = *divideTH1ByTH1(hEffi_AS_WW, hEffi_D_NW, "", "hEffi_AS_WW_over_D_2");

    auto drawAndAdd = [](TObject &o, TLegend *leg = nullptr)
    {
        o.Draw("same");
        if (leg)
            leg->AddEntry(&o, o.GetName(), "l");
    };
    // setup canvas
    TCanvas &c1 = *new TCanvas((fID + "compareMeasuredEffis_TH1_New").data(),
                               (fID + "compareMeasuredEffis_TH1_New").data(), 2000, 1000);

    c1.Divide(1, 2);
    c1.cd(1);
    gPad->SetLogy();
    TH2 *hd21 = new TH2F("hd21", ";pT (GeV);effis", 1, 0., 10., 1, 1e-5, 8e-3);
    hd21->Draw();

    auto leg = new TLegend(.73, .64, .90, .90, "");
    drawAndAdd(hEffi_AS_NW, leg);
    drawAndAdd(hEffi_AS_WW, leg);
    drawAndAdd(hEffi_D_NW, leg);
    leg->Draw("same");

    c1.cd(2);
    //~ gPad->SetLogy();
    TH2 *hd22 = new TH2F("hd22", ";pT (GeV);ratios", 1, 0., 10., 1, 0.5, 2.);
    hd22->Draw();

    auto leg2 = new TLegend();
    drawAndAdd(hEffi_AS_NW_over_D);
    drawAndAdd(hEffi_AS_WW_over_D);
    //~ leg2->Draw("same");

    saveCanvasAs(c1);
    return c1;
}


