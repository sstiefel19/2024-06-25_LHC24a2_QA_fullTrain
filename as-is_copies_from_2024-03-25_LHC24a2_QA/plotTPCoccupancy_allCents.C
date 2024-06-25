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

void plotTPCoccupancy(){
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    std::vector<std::string> vCents{"1013", "1131", "1353", "1591"};
    std::vector<std::string> vCents{"1013", "1131", "1353", "1591"};
    std::vector<infoPack> vInfosMB{};
    
    TCanvas* c = new TCanvas("plotTPCoccupancy", "plotTPCoccupancy", 2000, 1000);
    c->Divide(3,4);
    
    std::string dirnameMB("/trains/2024-02-26_allMCs_ptw_0b/LHC20e3/");
    std::string fname24("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/GCo_995.root");

    int iCent = 1;
    for (auto cent : vCents){
        bool my = !((iCent-1)%2);
        
        std::string fnameMB(dirnameMB + (!my ? "2050a/" : iCent==1 ? "ab/" : "ac/" ) + "GCo_994.root");
        GCo MB({fnameMB, "GammaConvV1_994/", cent+"0053", cent+"0053_0d200009ab770c00amd0404000_0152101500000000"});
        GCo AS({fname24, "GammaConvV1_995/", cent+"0023", cent+"0023_0d200009ab770c00amd0404000_0152101500000000"});
        
        TH1* hOut = (TH1*)MB.GetFromEvt("hTPCoutTracks");
        TH1* hClusters = (TH1*)MB.GetFromEvt("hTPCclusters");
        TH1* hOut24 = (TH1*)AS.GetFromEvt("hTPCoutTracks");
        TH1* hClusters24 = (TH1*)AS.GetFromEvt("hTPCclusters");
        
        TH1* hGoodESDtracks = (TH1*)MB.GetFromESD("GoodESDTracks");
        TH1* hGoodESDtracks24 = (TH1*)AS.GetFromESD("GoodESDTracks");
        
        
        c->cd(3*iCent-2);
        hOut->Draw(iCent>1 ? "same" : "");
        hOut24->SetLineColor(kRed);
        if (my)hOut24->Draw("same");
        
        c->cd(3*iCent-1);
        hClusters->Draw(iCent>1 ? "same" : "");
        hClusters24->SetLineColor(kRed);
        if (my) hClusters24->Draw("same");
        
        c->cd(3*iCent);
        hGoodESDtracks->Draw(iCent>1 ? "same" : "");
        hGoodESDtracks24->SetLineColor(kRed);
        if (my)hGoodESDtracks24->Draw("same");
        
        ++iCent;
    }

}











