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



struct GCo {
    std::string fname;
    std::string mDir;
    std::string evCut;
    std::string cutNo;
    
    TObject* GetFromESD(std::string oName){
        return getObjectFromPathInFile(fname, mDir+"Cut Number " + cutNo + "/" + cutNo + " ESD histograms/" + oName);}
    
    TObject* GetFromEvt(std::string oName){
        return getObjectFromPathInFile(fname, mDir+"Cut Number " + cutNo + "/ConvEventCuts_" + evCut + "/" + oName);}
};

template <typename T>
double median(const T* h1) { 

   int n = h1->GetXaxis()->GetNbins();  
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray(); 
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]); 
}

// TFormula* neutralsF  = new TFormula("neutrals",  "max(1.,470.*(x<5.)+62.*(x>7.5)*(x<12.5))");
void plotTPCoccupancy(){
    
    gROOT->Reset();
    //~ gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    std::vector<std::string> vCents{"1013", "1353"};
    
    
    TCanvas* c = new TCanvas("plotTPCoccupancy", "plotTPCoccupancy", 2000, 1000);
    c->Divide(3,2);
    
    std::string dirnameMB("/trains/2024-02-26_allMCs_ptw_0b/LHC20e3/");
    

    int iCent = 1;
    for (auto cent : vCents){
        
        std::string fnameMB(dirnameMB + (iCent==1 ? "ab/" : "ac/" ) + "GCo_994.root");
        GCo MB({fnameMB, "GammaConvV1_994/", cent+"0053", cent+"0053_0d200009ab770c00amd0404000_0152101500000000"});
    
        //~ std::string fname24(Form("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/2042_child%d/GCo_995.root", iCent));
        std::string fname24(Form("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_995.root"));
        
        GCo AS({fname24, "GammaConvV1_995/", cent+"0023", cent+"0023_0d200009ab770c00amd0404000_0152101500000000"});
        
        TH1* hOut = (TH1*)MB.GetFromEvt("hTPCoutTracks");
        TH1* hClusters = (TH1*)MB.GetFromEvt("hTPCclusters");
        TH1* hGoodESDtracks = (TH1*)MB.GetFromESD("GoodESDTracks");
        
        
        TH1* hOut24 = (TH1*)AS.GetFromEvt("hTPCoutTracks");
        TH1* hClusters24 = (TH1*)AS.GetFromEvt("hTPCclusters");
        TH1* hGoodESDtracks24 = (TH1*)AS.GetFromESD("GoodESDTracks");
        
        hOut24->SetTitle(Form("%s %s", buildCentString(cent).data(), hOut->GetTitle()));
        hClusters->SetTitle(Form("%s %s", buildCentString(cent).data(), hClusters->GetTitle()));
        hGoodESDtracks24->SetTitle(Form("%s %s", buildCentString(cent).data(), hGoodESDtracks->GetTitle()));
        
        
        Double_t scaleFactor = hClusters24->GetEntries() / (Double_t)hClusters->GetEntries();
        
        auto drawMean = [](TH1* h){
            Double_t x = h->GetMean(); 
            Double_t y = h->GetMaximum();
            auto l = new TLine(x,0.,x,y);
            l->SetLineColor(h->GetLineColor());
            printf("trying to to paint line with x,y: %f, %f\n", x, y);
            l->Draw("same");};
        
        hOut->Scale(scaleFactor);
        hClusters->Scale(scaleFactor);
        hGoodESDtracks->Scale(scaleFactor);
        
        auto makeLegend = [&iCent](TH1* hMB, TH1* hAS){
            auto l = new TLegend();
            l->AddEntry(hMB, "LHC20e3 scaled by n.o.e. LHC24a1");
            l->AddEntry(hAS, Form("LHC24a1 (number of injected pi0s, etas : %d)", iCent==1 ? 470 : 62));
            l->SetTextSize(0.03);
            l->Draw("same");
            
            hAS->SetName(Form("LHC24a1"));
        };
        
        c->cd(3*iCent-2);
        gPad->SetLogy();
        hOut24->Draw(iCent>1 ? "same" : "");
        hOut24->SetLineColor(kRed);
        hOut->Draw("same");
        makeLegend(hOut, hOut24);
        
        
        c->cd(3*iCent-1);
        //~ gPad->SetLogy();
        hClusters->Draw(iCent>1 ? "same" : "");
        hClusters24->SetLineColor(kRed);
        hClusters24->Draw("same");
        makeLegend(hOut, hOut24);

        
        //~ drawMean(hClusters);
        //~ drawMean(hClusters24);
        
        
        c->cd(3*iCent);
        gPad->SetLogy();
        hGoodESDtracks24->Draw(iCent>1 ? "same" : "");
        hGoodESDtracks24->SetLineColor(kRed);
        hGoodESDtracks->Draw("same");
        makeLegend(hOut, hOut24);

        
        ++iCent;
    }


    c->SaveAs(Form("%s.pdf", c->GetName()));
}











