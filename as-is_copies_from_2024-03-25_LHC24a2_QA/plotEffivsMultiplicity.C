#include "/analysisSoftware/SupportingMacros/utils_sstiefel_2023.h"

#include <iostream>
#include <string.h>
#include <vector>


#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"

struct infoPack {
    TH1* hEffi;
    Double_t multMeas;
    std::string label;
};



infoPack getTrueEfficiencyHisto(std::string meson, std::string cent, int nR = 1, bool theMB=true){
    
    
    std::string config(theMB ? "994" : meson=="Pi0" ? "995" : "996");
    std::string fname(theMB ? 
        Form("/trains/2023-11-30_mc_ptw_0b/%s/GCo_994.root", (cent=="101") ? "ab" :  
                                                             (cent=="135") ? "ac" : "2004a") : 
        Form("/trains/2024-01-29_LHC24a1_QA_ptw_0b/GCo_%s_%s.root", config.data(), meson.data()));
    std::string evtcut(cent + ((cent=="101" || cent=="135") ? "3" : "1") + "00" + (theMB ? "53" : "23"));
    
    std::string wholecut(evtcut+"_0d200009ab770c00amd0404000_0152101500000000");
    std::string histoname(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    std::string fullpathinfile(Form("GammaConvV1_%s/Cut Number %s/%s True histograms/%s", config.data(), wholecut.data(), wholecut.data(), histoname.data()));
    std::string fullpathinfileDenom(Form("GammaConvV1_%s/Cut Number %s/%s MC histograms/MC_%sInAcc_Pt", config.data(), wholecut.data(), wholecut.data(), meson.data()));
    std::string fullpathinfileGenWOW(Form("GammaConvV1_%s/Cut Number %s/%s MC histograms/MC_%sWOWeightInAcc_Pt", config.data(), wholecut.data(), wholecut.data(), meson.data()));
    TH1* hGoodESDTracks = (TH1*)getObjectFromPathInFile(fname, 
        Form("GammaConvV1_%s/Cut Number %s/%s ESD histograms/GoodESDTracks",  config.data(), wholecut.data(), wholecut.data()));
    Double_t nEvents = hGoodESDTracks->GetEntries();
    cout << nEvents << endl;
    
    Double_t lMeanESDTracks = hGoodESDTracks->GetMean();

    TH1* hMesonGenInAcc = (TH1*)getObjectFromPathInFile(fname,fullpathinfileDenom,"MB");
    TH1* hMesonWOWGenInAcc_notNorm = (TH1*)getObjectFromPathInFile(fname,fullpathinfileGenWOW,"MB");
    hMesonGenInAcc->Scale(1./nEvents);
    TH2* h2 = (TH2*)getObjectFromPathInFile(fname, fullpathinfile,"cloneMB");
    h2->Scale(1./nEvents);
    TH1* h1TrueMeson = h2->ProjectionX();
    
    h1TrueMeson->Rebin(nR);
    hMesonGenInAcc->Rebin(nR);
    hMesonWOWGenInAcc_notNorm->Rebin(nR);
    TH1* hMyTrueEffi = divideTH1ByTH1(*h1TrueMeson, *hMesonGenInAcc,"", (meson + "_" + cent + "_hMyTrueEffi").data());
    
    infoPack lRes;
    lRes.hEffi = hMyTrueEffi;
    lRes.multMeas = hGoodESDTracks->GetMean();
    lRes.label = buildCentString(cent);

    return lRes;
}

TGraphAsymmErrors* makeGraphAndPlot(std::vector<infoPack> &vInfoPerCent, int thePtBin){
    
    int nPoints = vInfoPerCent.size();
    // plot Text next to markers
    auto labelPoints = [&](TGraphAsymmErrors *g){
        Double_t x,y;
        TLatex *l;
        for (int i=0; i<nPoints; i++) {
            g->GetPoint(i,x,y);
            l = new TLatex();
            l->SetTextSize(0.035);
            l->SetTextFont(42);
            l->SetTextAlign(11);
            l->DrawLatex(x, y, Form(" %s", vInfoPerCent[i].label.data()));
        }
    };
    
    TH1 *hI = vInfoPerCent[0].hEffi;
    std::string id(Form("gPtBin_%d_%.1f-%.1f_GeV", thePtBin, hI->GetBinLowEdge(thePtBin), hI->GetBinLowEdge(thePtBin+1)));
    std::vector<Double_t> x;
    std::vector<Double_t> y;
    std::vector<Double_t> ey;
    for (auto p : vInfoPerCent){
        x.push_back(p.multMeas);
        y.push_back(p.hEffi->GetBinContent(thePtBin));
        ey.push_back(p.hEffi->GetBinError(thePtBin));
    }
    auto gr = new TGraphAsymmErrors(nPoints, x.data(), y.data(), nullptr, nullptr, ey.data(), ey.data());
    gr->SetName(id.data());
    gr->SetTitle(id.data());
    gr->GetXaxis()->SetTitle("Mean hGoodESDTracks");
    gr->GetYaxis()->SetTitle("True Efficiency");
    
    //~ TCanvas *c1 = new TCanvas(id.data(), id.data(), 2000, 1000);
    //~ gr->Draw("ALP");
    //~ labelPoints(gr);
    //~ c1->SaveAs(Form("png/%s.png", id.data()));
    return gr;
}

std::vector<std::vector<TF1*>>&
    getLinearExtrapolationsBetweenCentralities(std::string           meson,
                                               std::vector<infoPack> vInfoPerCent,
                                               int                   iMax) 
{
    std::vector<TGraphAsymmErrors*> vGraphs{};
    std::vector<Double_t> vMMs{}; // mins and maxes
    for (int iPt=1; iPt< iMax; ++iPt){
        vGraphs.push_back(makeGraphAndPlot(vInfoPerCent, iPt));
        vMMs.insert(vMMs.end(), {vGraphs.back()->GetYaxis()->GetXmin(), 
                                 vGraphs.back()->GetYaxis()->GetXmax()});
    }
    const auto [lGMin, lGMax] = std::minmax_element(begin(vMMs), end(vMMs));
    *lGMin = 0.;
    *lGMax = .0041;
    printf("GMin Max: %f and %f\n", *lGMin, *lGMax);
    
    // plot all fits in original dimensions
    TCanvas *c1 = new TCanvas(Form("plotEffiVsMultiplicity_%s", meson.data()),
                              Form("plotEffiVsMultiplicity_%s", meson.data()), 2000, 1000);
    auto g = (TGraphAsymmErrors*)vGraphs.front()->Clone("cloneForDrawingAxis");
    //~ g->SetMarkerColor(kWhite);
    g->SetTitle(Form("LHC20e3: %s true effis in diff. centralities (x-axis) for diff. pT bins (colors)", meson.data()));
    g->Draw("ALP");
    g->GetYaxis()->SetRangeUser(0., .002);
    
    auto leg = new TLegend(0.6, .63, .9,.92);
    leg->Draw("same");
    
    // outer over pt bins, inner over cent classes, starting with 010 - 3050. at last a fit over whole range 
    std::vector<std::vector<TF1*>> &vFits = *new std::vector<std::vector<TF1*>>{};
    
    int iPt =0;
    for (auto &ig : vGraphs){
        //~ if (ig!=vGraphs.front()){leg->AddEntry(ig, Form(), "l");}
        Color_t color = iPt;
        ig->SetLineColor(color);
        ig->SetMarkerColor(color);
        ig->Draw("SAME P");
        vFits.emplace_back(std::vector<TF1*>{});
        
        for (int iP=0; iP<ig->GetN(); ++iP){
            bool last = iP==ig->GetN()-1;
            Double_t* xs = ig->GetX();
            Double_t xmin = last ? xs[iP] : xs[iP+1];
            Double_t xmax = last ? xs[0]  : xs[iP];
            TF1* f = new TF1(Form("fit_%s_%d", ig->GetName(), iP), "[0] + [1]*x", xmin, xmax);
            ig->Fit(f,"NR");
            f->SetLineColor(ig->GetMarkerColor());
            if (!last && ig!=vGraphs.front() /*&& f->GetParameter(1)<=0*/){f->Draw("same");}
            vFits.back().push_back(f);
        }
        if (ig!=vGraphs.front()){
            TF1* f = vFits.back().back();
            leg->AddEntry(ig, f->GetName(), "lp");
        }

        ++iPt;
    }
    // plot Text next to markers
    auto labelPoints = [&](TGraphAsymmErrors *g){
        Double_t x,y;
        TLatex *l;
        for (int i=0; i<g->GetN(); i++) {
            g->GetPoint(i,x,y);
            l = new TLatex();
            l->SetTextSize(0.035);
            l->SetTextFont(42);
            l->SetTextAlign(12);
            l->DrawLatex(x-75, .0018, Form(" %s", vInfoPerCent[i].label.data()));
        }
    };
    labelPoints(vGraphs.back());
    
    
    c1->SaveAs(Form("%s.png", c1->GetName()));

    if (false)
    {
        // move all fits upwards such they intersect in the middle between 
        TCanvas *c2 = new TCanvas((meson+"_compareSlopes").data(),(meson+"_compareSlopes").data(), 2000, 1000);
        auto g2 = (TGraphAsymmErrors*)g->Clone("g2");
        g2->Draw("AP");
        g2->GetYaxis()->SetRangeUser(.0012, .0030);
        
        Double_t yFixPoint = 0.5 * (*lGMin + *lGMax);    
        for (int iM=0;iM < vFits.front().size()-1; ++iM){
            TF1* fFirst = vFits[0][iM];
            Double_t xFixPoint = 0.5 * (fFirst->GetXmin() + fFirst->GetXmax());
            for (int iPt=0; iPt < vFits.size(); ++iPt){
                TF1 *fOrig = vFits[iPt][iM];
                TF1 *fUp = (TF1*)fOrig->Clone(Form("%s_upped", fOrig->GetName()));
                Double_t m = fUp->GetParameter(1);
                fUp->SetParameter(0, yFixPoint - m*xFixPoint);
                
                fUp->SetLineColor(kBlue+iPt);
                fUp->Draw("same");
            }
        }
        c2->SaveAs(Form("%s.png", c2->GetName()));
    }
    return vFits;
}

void drawAndAdd(TH1* h, TLegend* leg=nullptr, const char* str=nullptr){
    h->Draw("same");
    if (leg) {leg->AddEntry(h, str ? str : h->GetName());}
}

void plotEffivsMultiplicity_(std::string meson, 
                             std::string centPlot,
                             std::vector<infoPack> &vInfosMB,
                             int nRebin){
    int lUseWhichFit = centPlot=="101" ? 0 : 1;
    int iCent =        centPlot=="101" ? 0 : 2;
    
    int iMax = vInfosMB[0].hEffi->FindBin(8.99)+1;        
    std::vector<std::vector<TF1*>> &vFits = getLinearExtrapolationsBetweenCentralities(meson, vInfosMB, iMax);
    
    std::string nameCanvas(Form("plotEffivsMultiplicity_extrapolated-effis_%s_%s", meson.data(), centPlot.data()));
    auto c = new TCanvas(nameCanvas.data(), nameCanvas.data(), 2000, 1000);
    c->Divide(1,2);
    c->cd(1);
    gPad->SetLogy();
    TH2* hd4 = new TH2F("hd4",Form("%s %s true efficiencies;MC pT(GeV);True MesonEffi", 
                                   buildCentString(centPlot).data(), meson.data()),1,0.,10.,1.,1e-5,2e-3);
    hd4->Draw();
    auto leg = new TLegend();
    leg->Draw("same");
    
    infoPack infoAS = getTrueEfficiencyHisto(meson, centPlot, nRebin, false);
    TH1* hEffiMB = vInfosMB[iCent].hEffi;
    TH1* hEffiAS = infoAS.hEffi;
    hEffiAS->SetLineColor(kRed);
    
    drawAndAdd(hEffiMB, leg, "MB LHC20e3");
    drawAndAdd(hEffiAS, leg, "AS LHC24a1");
    
    auto createEffiHistoFromFits = [&vFits](int iFit, Double_t atMult, TH1* theHForB){
        printf("createEffiHistoFromFits with iFit %d at multi = %f\n", iFit, atMult);
        TH1* hEffiFromFit = cloneTH1(*theHForB, "_extr");
        for (int iPt=0; iPt < vFits.size(); ++iPt){
            TF1* f = vFits[iPt][iFit];
            printf("evaluating fit %s at multi = %f\n", f->GetName(), atMult);
            hEffiFromFit->SetBinContent(iPt+1, f->Eval(atMult)); 
            hEffiFromFit->SetBinError(iPt+1, 0.);
        }
        return hEffiFromFit;
    };
    
    Double_t correctedMulti = vInfosMB[iCent].multMeas + infoAS.multMeas;
    TH1* hEffiMB_ext = createEffiHistoFromFits(lUseWhichFit, correctedMulti, hEffiMB);
    hEffiMB_ext->SetLineColor(kBlack);
    drawAndAdd(hEffiMB_ext, leg, Form("MB extrapolated to 'Mean hGoodESDTracks' = %.1f", correctedMulti));

    
    c->cd(2);
    gPad->SetGridy();
    TH2* hd5 = new TH2F("hd5",Form("normalized by MB ;MC pT(GeV); / MB"),1,0.,10.,1.,0.6,1.4   );
    hd5->Draw();
    auto leg2 = new TLegend();
    leg2->Draw("same");
    
    TH1* hASoverMB = divideTH1ByTH1(*hEffiAS, *hEffiMB, "", "AS");
    drawAndAdd(hASoverMB, leg2);
    
    TH1* hExtoverMB = divideTH1ByTH1(*hEffiMB_ext, *hEffiMB, "", "MB extrapolated");
    reduceErrorsByFactor(*hExtoverMB);
    hExtoverMB->SetLineColor(kBlack);
    drawAndAdd(hExtoverMB, leg2);
    
    c->SaveAs(Form("%s.png", c->GetName()));
}



void plotEffivsMultiplicity_(std::string meson){
    
    int nRebin =2;
    
    std::vector<std::string> vCents{"101", "113", "135", "159"};
    std::vector<infoPack> vInfosMB{};
    for (auto cent : vCents){
        vInfosMB.push_back(getTrueEfficiencyHisto(meson, cent, nRebin));
    }

    plotEffivsMultiplicity_(meson, "101", vInfosMB, nRebin);
    plotEffivsMultiplicity_(meson, "135", vInfosMB, nRebin);

}

void plotEffivsMultiplicity(){
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    plotEffivsMultiplicity_("Pi0");    
}
