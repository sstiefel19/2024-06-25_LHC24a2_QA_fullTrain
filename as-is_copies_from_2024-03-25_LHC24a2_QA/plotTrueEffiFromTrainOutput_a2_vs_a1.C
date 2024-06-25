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

void plotTrueEffiFromTrainOutput_(std::string meson="Pi0", std::string cent="101"){
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    std::string configA( meson=="Pi0" ? "995" : "996");
    //std::string configA( "997");
    
    int nR = meson=="Pi0" ? 3 : 2;
    
    
    std::string evtcutA(cent + "30023");
    
    std::string fnameA1(Form("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/%s/GCo_%s.root", 
        cent=="101" ? "child_1_010" : "child_2_3050", configA.data()));
    std::string fnameA2(Form("/trains/2024-03-26_LHC24a2_ptw_0b/%s/GCo_%s.root", 
        cent=="101" ? "child_1_010" : "child_2_3050", configA.data()));
    
    std::string restCutNo("_0d200009ab770c00amd0404000_0152101500000000");
    
    GCo AS1(fnameA1, Form("GammaConvV1_%s/", configA.data()), evtcutA, restCutNo);
    GCo AS2(fnameA2, Form("GammaConvV1_%s/", configA.data()), evtcutA, restCutNo);
    
    auto scaleRebinDiv = [nR](TH1 &h, double nEvts, const char *tag){
        h.Scale(1./nEvts);
        h.Rebin(nR);
        return divideTH1ByBinWidths(h, tag);
    };
    
    //AS1
    Double_t nEventsA1       = ((TH1*)AS1.GetFromESD("VertexZ"))->GetEntries(); cout << nEventsA1 << endl;
    TH1* hMesonGenInAccA1    = (TH1*)AS1.GetFromMC(Form("MC_%sInAcc_Pt", meson.data()));
    TH1* hMesonWOWGenInAccA1 = (TH1*)AS1.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data()));    
    TH2* h2A1                = (TH2*)AS1.GetFromTrue(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    TH1* h1TrueMesonA1 = h2A1->ProjectionX();
    
    hMesonGenInAccA1 = scaleRebinDiv(*hMesonGenInAccA1, nEventsA1, "A1");
    h1TrueMesonA1 = scaleRebinDiv(*h1TrueMesonA1, nEventsA1, "A1");

    //AS2
    Double_t nEventsA2       = ((TH1*)AS2.GetFromESD("VertexZ"))->GetEntries(); cout << nEventsA2 << endl;
    TH1* hMesonGenInAccA2    = (TH1*)AS2.GetFromMC(Form("MC_%sInAcc_Pt", meson.data()));
    TH1* hMesonWOWGenInAccA2 = (TH1*)AS2.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data()));    
    TH2* h2A2                = (TH2*)AS2.GetFromTrue(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    TH1* h1TrueMesonA2 = h2A2->ProjectionX();    
    
    hMesonWOWGenInAccA2 = scaleRebinDiv(*hMesonWOWGenInAccA2, nEventsA2, "A2");
    hMesonGenInAccA2 = scaleRebinDiv(*hMesonGenInAccA2, nEventsA2, "A2");
    h1TrueMesonA2 = scaleRebinDiv(*h1TrueMesonA2, nEventsA2, "A2");
    
    
    std::string nameC(Form("effis_from_trainoutput_a2_vs_a1_%s_%s_reb%d_%s", meson.data(), cent.data(), nR, configA.data()));
    TCanvas *c1 = new TCanvas(nameC.data(), nameC.data(), 2000, 1000);
    c1->Divide(2, 2, 0.01, 0.002);
    
    c1->cd(1);
    gPad->SetLogy();
    TH2* hd1 = new TH2F("hd1",Form("%s: True%sYields (x-proj. of ESD_TruePrimary%s_MCPt_ResolPt);MC pT(GeV);dN/dpT (1/Nevt.)", buildCentString(evtcutA).data(), meson.data(), meson.data()),1,0.,10.,1.,1e-6,4e-1);
    hd1->Draw();
    std::string tagASMC(cent=="101" ? "LHC24a1a2 AS MC WW" : "LHC24a1b2 AS MC WW");
    
    {
        auto *leg = new TLegend();
        drawAndAdd(*h1TrueMesonA1, "same", kBlue, leg, "A1 WW");
        drawAndAdd(*h1TrueMesonA2, "same", kRed, leg, "A2 WW");
        leg->Draw("same");
    }   
    
    c1->cd(2);
    gPad->SetLogy();
    TH2* hd2 = new TH2F("hd2",Form("True%sYields / MC_%sInAcc_Pt;MC pT(GeV);trueEffi", meson.data(), meson.data()),1,0.,10.,1.,1e-5,2e-3);
    hd2->Draw();
    TH1* hMyTrueEffiA1 = divideTH1ByTH1(*h1TrueMesonA1, *hMesonGenInAccA1,"", "hMyTrueEffiA1","hMyTrueEffiA1");
    TH1* hMyTrueEffiA2 = divideTH1ByTH1(*h1TrueMesonA2, *hMesonGenInAccA2,"", "hMyTrueEffiA2","hMyTrueEffiA2");
    
    {
        auto *leg = new TLegend();
        drawAndAdd(*hMyTrueEffiA1, "same", kBlue, leg, "A1 WW");
        drawAndAdd(*hMyTrueEffiA2, "same", kRed, leg,  "A2 WW");
        leg->Draw("same");
    }

    c1->cd(3);
    gPad->SetLogy();
    gPad->SetGridy();
    TH2* hd3 = new TH2F("hd3",Form("MC_%sInAcc_Pt (generated %ss in accep.);MC pT(GeV);(pt-weighted) counts (1/Nevt.)", meson.data(), meson.data()),1,0.,10.,1.,1e-2,1e3);
    hd3->Draw();    
    {
        auto *leg = new TLegend();
        drawAndAdd(*hMesonGenInAccA1, "same", kBlue, leg, "A1 WW");
        drawAndAdd(*hMesonGenInAccA2, "same", kRed, leg, "A2 WW");
        leg->Draw("same");
    }
    
    c1->cd(4);
    //~ gPad->SetLogy();
    TH2* hd4 = new TH2F("hd4",";MC pT(GeV);AS2 / AS1",1,0.,10.,1.,0.5,2.);
    hd4->Draw();
    TH1* hEffiA2overA1 = divideTH1ByTH1(*hMyTrueEffiA2, *hMyTrueEffiA1, "", "hEffiA2overA1");
    hEffiA2overA1->Draw("same");
    TH1* hGenA2overGenA1 = divideTH1ByTH1(*hMesonGenInAccA2, *hMesonGenInAccA1, "","hGenA2overGenA1");
    hGenA2overGenA1->SetLineColor(kBlue);
    hGenA2overGenA1->Draw("same");
    
    {
        auto *leg = new TLegend();
        leg->AddEntry(hEffiA2overA1, "true meson efficiencies");
        leg->AddEntry(hGenA2overGenA1, Form("MC_%sInAcc_Pt", meson.data()));
        leg->Draw("same");
    }
    
    gPad->SetGridy();
    
    // gSystem->Exec("mkdir root");
    // TFile outfile(Form("root/%s.root", nameC.data()), "RECREATE");
    // hMesonGenInAcc->Write();
    // h1TrueMeson->Write();
    // hMyTrueEffi->Write();
    // hMesonGenInAccA->Write();
    // h1TrueMesonA->Write();
    // hMyTrueEffiA->Write();
    // hEffiAoverMB->Write();
    // hGenAoverGenMB->Write();
    // c1->Write();
    // outfile.Close();
    
    saveCanvasAs(*c1, "png");
}

void plotTrueEffiFromTrainOutput_a2_vs_a1(std::string meson="Pi0", std::string cent="101"){
        
    // plotTrueEffiFromTrainOutput_("Pi0", "101");
    // return;
    
    
    for (auto m : std::vector<std::string>({"Pi0", "Eta"})){
        
        for (auto c : std::vector<std::string>({"101", "135"})){
                
            plotTrueEffiFromTrainOutput_(m, c);
        }
    }
    
}
