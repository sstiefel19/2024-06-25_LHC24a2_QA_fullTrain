#include "/analysisSoftware/utils_sstiefel_2024/include/GCo.h"           
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_computational.h" 
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_files_strings.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_fits.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_plotting.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_utils.h"
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TF1.h"           
#include "/analysisSoftware/utils_sstiefel_2024/include/utils_TH1.h"

#include <iostream>
#include <string.h>
#include <vector>

#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"

void plotTrueEffiFromTrainOutput_(std::string meson, 
                                  std::string cent,
                                  std::string configA, 
                                  int nR=2,
                                  bool keepFilesOpen=true){
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    std::string restCutNo("_0d200009ab770c00amd0404000_0152101500000000");
    std::string evtcut(cent + "30053");
    std::string evtcutA(cent + "30023");
    
    // set filenames to be used here
    std::string fname(Form("/trains/2023-11-30_mc_ptw_0b/%s/GCo_994.root", 
        (cent=="101") ? "ab" : "ac"));

    std::string fnameA(
        Form("/trains/2024-06-25_invPtW_all-AS-MCs/2101_20240624-1718/merge_runlist_1/root_archive/GCo_%s.root", 
            configA.data()));
   
    GCo MB(fname, "GammaConvV1_994/", evtcut, restCutNo, keepFilesOpen);
    GCo AS(fnameA, Form("GammaConvV1_%s/", configA.data()), evtcutA, restCutNo, keepFilesOpen);
        
    // MB
    Double_t nEvents    = ((TH1*)MB.GetFromESD("VertexZ"))->GetEntries(); cout << nEvents << endl;
    TH1* hMesonWOWGenInAcc = (TH1*)MB.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data()));    
    hMesonWOWGenInAcc->Rebin(nR);
    
    TH2* h2_wow             = (TH2*)MB.GetFromTrue(Form("ESD_TruePrimaryMotherW0Weights_InvMass_Pt"));
    TH1* h1TrueMeson_wow    = h2_wow->ProjectionY();    
    h1TrueMeson_wow->Rebin(nR);
    
    //AS
    Double_t nEventsA       = ((TH1*)AS.GetFromESD("VertexZ"))->GetEntries(); cout << nEventsA << endl;
    TH1* hMesonWOWGenInAccA = (TH1*)AS.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data())); 
    hMesonWOWGenInAccA->Rebin(nR);   
    
    TH2* h2A                = (TH2*)AS.GetFromTrue(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    TH1* h1TrueMeson_wowA = h2A->ProjectionX();    
    h1TrueMeson_wowA->Rebin(nR);
    
        
    std::string nameC(Form("unweighted-counts_from_trainoutput_%s_%s_reb%d_%s", meson.data(), cent.data(), nR, configA.data()));
    TCanvas *c1 = new TCanvas(nameC.data(), nameC.data(), 2000, 1000);
    c1->Divide(1, 2, 0.01, 0.002);
    
    c1->cd(1);
    gPad->SetLogy();
    TH2* hd1 = new TH2F("hd1",
                        Form("%s: True%sYields (y-proj. of ESD_TruePrimaryMotherW0Weights_InvMass_Pt;MC pT(GeV);N", 
                             utils_files_strings::BuildCentString(evtcut).data(), meson.data()),
                        1,0.,10.,1.,1,4e5);
    hd1->Draw();
        std::string tagASMC(cent=="101" ? "LHC24a1a2 AS MC WW" : "LHC24a1b2 AS MC WW");
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*h1TrueMeson_wow, 
                                   "same", 
                                   kBlue,
                                   1., 
                                   leg, 
                                   "LHC20e3ab MB MC WW",
                                   "lp");
        utils_plotting::DrawAndAdd(*h1TrueMeson_wowA, 
                                    "same", 
                                    kRed,
                                    1., 
                                    leg, 
                                    tagASMC,
                                    "lp");                                   
        leg->Draw("same");

        auto &lPav = utils_plotting::SetupTPaveTextAndAddOneLine(nameC, 0.5, 0.8, 0.9, 0.9, 0.04);
        lPav.Draw("same");
    }   
    
    c1->cd(2);
    gPad->SetLogy();
    
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*hMesonWOWGenInAcc, "same", kBlue, 1., leg, "LHC20e3ab MB MC WW", "lp");
        utils_plotting::DrawAndAdd(*hMesonWOWGenInAccA, "same", kRed, 1., leg, tagASMC, "lp");
        leg->Draw("same");
        auto &lPav = utils_plotting::SetupTPaveTextAndAddOneLine(nameC, 0.5, 0.8, 0.9, 0.9, 0.04);
        lPav.Draw("same");
    }
    
    gPad->SetGridy();
    
    //gSystem->Exec("mkdir root");
    //TFile outfile(Form("root/%s.root", nameC.data()), "RECREATE");
    //hMesonGenInAcc->Write();
    //h1TrueMeson_wow->Write();
    //hMyTrueEffi->Write();
    //hMesonGenInAccA->Write();
    //h1TrueMeson_wowA->Write();
    //hMyTrueEffiA->Write();
    //hEffiAoverMB->Write();
    //hGenAoverGenMB->Write();
    //c1->Write();
    //outfile.Close();
    
    utils_plotting::SaveCanvasAs(*c1, "png");

}

void plotUnweightedCountsForAssesingUncs(){
        
    int nR = 1;

    for (auto m : std::vector<std::string>({"Pi0"})){
        for (auto c : std::vector<std::string>({"135"})){
        
    //for (auto m : std::vector<std::string>({"Pi0", "Eta"})){
        //for (auto c : std::vector<std::string>({"101", "135"})){
        
            
            // select configAs based on cent and meson

            // std::vector<std::string> configAs = (m=="Pi0") 
            //     ? {"997"}
            //     : (cent=="101")
            //         ? {"998"}
                    // : {"999"};
            std::string configA = (m=="Pi0") 
                ? "997"
                : (c=="101")
                    ? "998"
                    : "999";
            plotTrueEffiFromTrainOutput_(m, c, configA, nR);     
        }
    }
}
