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

void plotTrueEffiFromTrainOutput_(std::string meson="Pi0", 
                                  std::string cent="101",
                                  bool use997=false, 
                                  int nR=2,
                                  bool keepFilesOpen=true){
    
    gROOT->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    std::string restCutNo("_0d200009ab770c00amd0404000_0152101500000000");
    std::string evtcut(cent + "30053");
    std::string evtcutA(cent + "30023");
    std::string configA( use997 ? "997" : meson=="Pi0" ? "995" : "996");
    
    // config filenames to use here
    std::string fname(Form("/trains/2023-11-30_mc_ptw_0b/%s/GCo_994.root", 
        (cent=="101") ? "ab" : "ac"));
    std::string fnameA(Form("/trains/2024-03-26_LHC24a2_ptw_0b/%s/GCo_%s.root", 
     cent=="101" ? "child_1_010" : "child_2_3050", configA.data()));
    
    // old a1
    //std::string fnameA(Form("/trains/2024-02-26_allMCs_ptw_0b/LHC24a1/%s/GCo_%s.root", 
    //    cent=="101" ? "child_1_010" : "child_2_3050", configA.data()));
   

    GCo MB(fname, "GammaConvV1_994/", evtcut, restCutNo, keepFilesOpen);
    GCo AS(fnameA, Form("GammaConvV1_%s/", configA.data()), evtcutA, restCutNo, keepFilesOpen);
    
    auto scaleRebinDiv = [nR](TH1 &h, double nEvts, const char *tag){
        h.Scale(1./nEvts);
        h.Rebin(nR);
        return utils_TH1::DivideTH1ByBinWidths(h, tag, "", "");
    };
    
    // MB
    TH1* hMesonGenInAcc            = (TH1*)MB.GetFromMC(Form("MC_%sInAcc_Pt", meson.data()));
    TH1* hMesonWOWGenInAcc_notNorm = (TH1*)MB.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data()));    
    hMesonWOWGenInAcc_notNorm->Rebin(nR);
    TH2* h2             = (TH2*)MB.GetFromTrue(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    TH1* h1TrueMeson    = h2->ProjectionX();    
    Double_t nEvents    = ((TH1*)MB.GetFromESD("VertexZ"))->GetEntries(); cout << nEvents << endl;

    hMesonGenInAcc = scaleRebinDiv(*hMesonGenInAcc, nEvents, "MB");
    h1TrueMeson = scaleRebinDiv(*h1TrueMeson, nEvents, "MB");

    //AS
    Double_t nEventsA       = ((TH1*)AS.GetFromESD("VertexZ"))->GetEntries(); cout << nEventsA << endl;
    TH1* hMesonGenInAccA    = (TH1*)AS.GetFromMC(Form("MC_%sInAcc_Pt", meson.data()));
    TH1* hMesonWOWGenInAccA = (TH1*)AS.GetFromMC(Form("MC_%sWOWeightInAcc_Pt", meson.data()));    
    TH2* h2A                = (TH2*)AS.GetFromTrue(Form("ESD_TruePrimary%s_MCPt_ResolPt", meson.data()));
    TH1* h1TrueMesonA = h2A->ProjectionX();    
    
    hMesonWOWGenInAccA = scaleRebinDiv(*hMesonWOWGenInAccA, nEventsA, "AS");
    hMesonGenInAccA = scaleRebinDiv(*hMesonGenInAccA, nEventsA, "AS");
    h1TrueMesonA = scaleRebinDiv(*h1TrueMesonA, nEventsA, "AS");
        
    std::string nameC(Form("effis_from_trainoutput_%s_%s_reb%d_%s", meson.data(), cent.data(), nR, configA.data()));
    TCanvas *c1 = new TCanvas(nameC.data(), nameC.data(), 2000, 1000);
    c1->Divide(2, 2, 0.01, 0.002);
    
    c1->cd(1);
    gPad->SetLogy();
    TH2* hd1 = new TH2F("hd1",
                        Form("%s: True%sYields (x-proj. of ESD_TruePrimary%s_MCPt_ResolPt);MC pT(GeV);dN/dpT (1/Nevt.)", 
                             utils_files_strings::BuildCentString(evtcut).data(), meson.data(), meson.data()),
                        1,0.,10.,1.,1e-6,4e-1);
    hd1->Draw();
        std::string tagASMC(cent=="101" ? "LHC24a1a2 AS MC WW" : "LHC24a1b2 AS MC WW");
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*h1TrueMeson, 
                                   "same", 
                                   kBlue,
                                   1., 
                                   leg, 
                                   "LHC20e3ab MB MC WW",
                                   "lp");
        utils_plotting::DrawAndAdd(*h1TrueMeson, 
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
    // gPad->SetLogy();
    TH2* hd2 = new TH2F("hd2",Form("True%sYields / MC_%sInAcc_Pt;MC pT(GeV);trueEffi", meson.data(), meson.data()),1,0.,10.,1.,1e-5,2e-3);
    hd2->Draw();
    TH1* hMyTrueEffi = utils_TH1::utils_TH1::DivideTH1ByTH1(*h1TrueMeson, *hMesonGenInAcc,"", "hMyTrueEffi","hMyTrueEffi");
    TH1* hMyTrueEffiA = utils_TH1::utils_TH1::DivideTH1ByTH1(*h1TrueMesonA, *hMesonGenInAccA,"", "hMyTrueEffiA","hMyTrueEffiA");
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*hMyTrueEffi, "same", kBlue, 1., leg, "LHC20e3ab MB MC WW", "lp");
        utils_plotting::DrawAndAdd(*hMyTrueEffiA, "same", kRed, 1., leg, tagASMC, "lp");
        leg->Draw("same");
        auto &lPav = utils_plotting::SetupTPaveTextAndAddOneLine(nameC, 0.5, 0.8, 0.9, 0.9, 0.04);
        lPav.Draw("same");
    }

    c1->cd(3);
    gPad->SetLogy();
    gPad->SetGridy();
    TH2* hd3 = new TH2F("hd3",Form("MC_%sInAcc_Pt (generated %ss in accep.);MC pT(GeV);(pt-weighted) counts (1/Nevt.)", meson.data(), meson.data()),1,0.,10.,1.,1e-2,1e3);
    hd3->Draw();
    TH1*  hTrueMeson_estimateWOW = utils_utils::CloneTH1(*hMesonWOWGenInAcc_notNorm, nullptr, "hTrueMeson_estimateWOW_notNorm");
    if (!hTrueMeson_estimateWOW->Multiply(hMyTrueEffiA)) {cout << "\n\n\n MULTIPLICATION FAILED\n\n\n\n";}
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*hMesonGenInAcc, "same", kBlue, 1., 
                                   leg, "LHC20e3ab MB MC WW");
        utils_plotting::DrawAndAdd(*hMesonGenInAccA, "same", kRed, 1., 
                                   leg, tagASMC);
        utils_plotting::DrawAndAdd(*hTrueMeson_estimateWOW, "same", kBlack, 1., 
                                   leg, Form("MB MC: est. total number of true %ss wo/ weights.", meson.data()));
        leg->Draw("same");
        auto &lPav = utils_plotting::SetupTPaveTextAndAddOneLine(nameC, 0.5, 0.8, 0.9, 0.9, 0.04);
        lPav.Draw("same");
    }

    // c2
    {
        // TCanvas *c2 = new TCanvas("c2", "c2", 2000, 1000);
        // //~ c2->Divide(1,2);
        // //~ c2->cd(1)
        // gPad->SetLogy();
        // TH2* hd21 = new TH2F("hd21",";MC pT(GeV);dN/dpT",1,0.,10.,1.,1e-6,1e4);
        // hd21->Draw();
    
        // TH1* h1TrueMeson_diff = utils_TH1::DivideTH1ByBinWidths(*h1TrueMeson);
        // TH1* hMesonGenInAcc_diff = utils_TH1::DivideTH1ByBinWidths(*hMesonGenInAcc);
        // TH1* h1TrueMesonA_diff = utils_TH1::DivideTH1ByBinWidths(*h1TrueMesonA);
        // TH1* hMesonGenInAccA_diff = utils_TH1::DivideTH1ByBinWidths(*hMesonGenInAccA);
        // TH1* hMesonWOWGenInAccA_diff = utils_TH1::DivideTH1ByBinWidths(*hMesonWOWGenInAccA);
        
        // auto leg = new TLegend();
        // utils_plotting::DrawAndAdd(*hMesonGenInAcc_diff, "same", kBlue, 1., leg);
        // utils_plotting::DrawAndAdd(*hMesonWOWGenInAccA_diff, "same", kBlack, 1., leg);
        // utils_plotting::DrawAndAdd(*hMyTrueEffiA, "same", kRed, 1., leg);
        // leg->Draw("same");
        
        //~ hMesonGenInAcc_diff->Draw("same");
        //~ hMesonWOWGenInAccA_diff->Draw("same");
        //~ h1TrueMeson_diff->Draw("same");
        //~ hMesonGenInAccA_diff->Draw("same");
        //~ h1TrueMesonA_diff->Draw("same");
        //~ hMyTrueEffiA->Draw("same");
        
    }
    
    c1->cd(4);
    //~ gPad->SetLogy();
    TH2* hd4 = new TH2F("hd4",";MC pT(GeV);AS / MB",1,0.,10.,1.,0.5,2.);
    hd4->Draw();
    TH1* hEffiAoverMB = utils_TH1::utils_TH1::DivideTH1ByTH1(*hMyTrueEffiA, *hMyTrueEffi, "", "hEffiAoverMB", "hEffiAoverMB");
    hEffiAoverMB->Draw("same");
    TH1* hGenAoverGenMB = utils_TH1::utils_TH1::DivideTH1ByTH1(*hMesonGenInAccA, *hMesonGenInAcc, "","hGenAoverGenMB", "hGenAoverGenMB");
    
    {
        auto *leg = new TLegend();
        utils_plotting::DrawAndAdd(*hEffiAoverMB, "same", kRed, 1., leg, "true meson efficiencies");
        utils_plotting::DrawAndAdd(*hGenAoverGenMB, "same", kBlue, 1., leg, Form("MC_%sInAcc_Pt", meson.data()));
        leg->Draw("same");

        auto &lPav = utils_plotting::SetupTPaveTextAndAddOneLine(nameC, 0.5, 0.8, 0.9, 0.9, 0.04);
        lPav.Draw("same");
    }
    
    gPad->SetGridy();
    
    gSystem->Exec("mkdir root");
    TFile outfile(Form("root/%s.root", nameC.data()), "RECREATE");
    hMesonGenInAcc->Write();
    h1TrueMeson->Write();
    hMyTrueEffi->Write();
    hMesonGenInAccA->Write();
    h1TrueMesonA->Write();
    hMyTrueEffiA->Write();
    hEffiAoverMB->Write();
    hGenAoverGenMB->Write();
    c1->Write();
    outfile.Close();
    
    utils_plotting::SaveCanvasAs(*c1, "png");

}

void plotTrueEffiFromTrainOutput(){
        
    bool use997 = false;
    int nR = 2;
    //plotTrueEffiFromTrainOutput_("Pi0", "101", use997, nR);
    //return;
    
    for (auto m : std::vector<std::string>({"Pi0", "Eta"})){
    //for (auto m : std::vector<std::string>({"Pi0"})){
        for (auto c : std::vector<std::string>({"101", "135"})){
        //for (auto c : std::vector<std::string>({"101"})){
            if (m=="Pi0"){
                //for (bool use997 : std::vector<bool>{false, true}){
                for (bool use997 : std::vector<bool>{false}){
                     //for (int nR : std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8}){
                    
                       plotTrueEffiFromTrainOutput_(m, c, use997, nR); 
                    //}
                }
            }
            else {
                plotTrueEffiFromTrainOutput_(m, c, false, nR);
            }    
        }
    }
    
}
