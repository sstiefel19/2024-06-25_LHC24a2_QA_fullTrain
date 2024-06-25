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


TCanvas* plotRatioDirectC(TH1 &hNom, TH1 &hDen, std::string theTitle="", std::string theLegStr="", float yMin=0.5, float yMax=2.){

    TH1* hRatioDirect = (TH1*)hNom.Clone("hRatioDirect");
    hRatioDirect->Divide(&hNom, &hDen,1.,1.);
    TCanvas* c3 = new TCanvas("c3","c3",800,600);
    
    TH2 *hDraw = new TH2F("hdr", theTitle.data(), 1, 0., 12., 1, yMin, yMax);
    hDraw->SetStats(0);
    hDraw->Draw();
    hRatioDirect->Draw("same");
    gPad->SetGridy();
    TLegend* leg = new TLegend(0.2,0.8,0.9,0.9);
    //leg->SetHeader("legend title","C"); // option "C" allows to center the header
    leg->AddEntry(hRatioDirect, theLegStr.data());
    leg->Draw();
    
    return c3;
}


void compare_mcs_after_afterburner(){
    
    //~ std::string evtcut("10130e03");
    //~ std::string meson("Eta");
    std::string histoName("TrueMesonEffiPt");
    
    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    c->Divide(2,2);
    
    int iPad=0;
    for (auto meson : std::vector<std::string>({"Pi0", "Eta"})){
    
        for (auto evtcut : std::vector<std::string>({"10130e03", "13530e03"})){
    
            TH1* h24 = (TH1*)getObjectFromPathInFile(Form("/afterburner/2024-01-30_MBptw0b_a24ptw0b/standard/%s_0d200009ab770c00amd0404000_0152101500000000/PbPb_5.02TeV/%s_MC_GammaConvV1CorrectionHistosAddSig_%s_0d200009ab770c00amd0404000_0152101500000000.root",evtcut.data(), meson.data(), evtcut.data()), histoName.data());
            
            TH1* hMB = (TH1*)getObjectFromPathInFile(Form("/afterburner/2024-01-30_MBptw0b_a24ptw0b/standard/%s_0d200009ab770c00amd0404000_0152101500000000/PbPb_5.02TeV/%s_MC_GammaConvV1CorrectionHistosMinBias_%s_0d200009ab770c00amd0404000_0152101500000000.root",evtcut.data(), meson.data(), evtcut.data()), histoName.data());
    
            auto p1 = plotRatioDirectC(*h24, *hMB, Form(";pT (GeV);LHC24a1 / LHCe3a"), Form("%s TrueEfficiency %s", meson.data(), buildCentString(evtcut).data()));
    
            ++iPad;
            c->cd(iPad);
            p1->DrawClonePad();
            
        }
    }
    
    
    
    
    
}
