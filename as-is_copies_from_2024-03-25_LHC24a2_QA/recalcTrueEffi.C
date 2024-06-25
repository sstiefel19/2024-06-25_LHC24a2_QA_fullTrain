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


void recalcTrueEffi(){
    
    //~ std::string evtcut("10130e03");
    //~ std::string meson("Eta");
    std::string histoName("TrueMesonEffiPt");
    
    TCanvas *c = new TCanvas("c", "c", 1000, 1000);
    //~ gPad->SetLogy();
    //~ c->Divide(1,2);
    
    int iPad=0;
    for (auto meson : std::vector<std::string>({"Pi0"})){
    
        for (auto evtcut : std::vector<std::string>({"10130e03", "13530e03"})){
            
            std::string fname("/afterburner/2023_11_24_LHC20e3ab_wPtw_0b/10130e03_0d200009ab770c00amd0404000_0152101500000000/PbPb_5.02TeV/Pi0_MC_GammaConvV1WithoutCorrection_10130e03_0d200009ab770c00amd0404000_0152101500000000.root");
    
            TH1* hTrueYieldDiff = (TH1*)getObjectFromPathInFile(fname, "histoYieldTrueMeson");
            TH1* hTrueYieldCounts = multiplyTH1ByBinWidths(*hTrueYieldDiff, "", "Pi0_TrueYieldCounts");
    
            //~ TH1* hGen = (TH1*)getObjectFromPathInFile(fname, "MC_Pi0InAcc_Pt");
            //~ TH1* hGenCounts = multiplyTH1ByBinWidths(*hGen);

            //~ TH1* hTrueEffiB = cloneTH1(*hGen, nullptr, "hTrueEffiB");
            //~ hTrueEffiB->Divide(hTrueYieldCounts, hGenCounts, 1., 1., "B");
            
            //~ TH1* hTrueEffi = (TH1*)getObjectFromPathInFile("/afterburner/2023_11_24_LHC20e3ab_wPtw_0b/10130e03_0d200009ab770c00amd0404000_0152101500000000/PbPb_5.02TeV/Pi0_MC_GammaConv_OnlyCorrectionFactor_10130e03_0d200009ab770c00amd0404000_0152101500000000.root", "TrueMesonEffiPt");
            //~ hTrueEffi->SetLineColor(kRed);
            
            

            //~ hTrueEffiB->DrawClone();
            //~ hTrueEffi->DrawClone("same");
            ++iPad;
            //~ c->cd(iPad);
            
            hTrueYieldCounts->DrawClone();
            gPad->SetLogy();
            c->SaveAs("Pi0_TrueYieldCounts_LHC20e3ab.pdf");
            break;
        }
        break;      
    }
    
    
    
    
    
}
