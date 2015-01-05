#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLorentzVector.h"

#include "classes/DelphesClasses.h"


/////////////////////////////////////////////////////                                                                                                                        

int main (int argc, char** argv){

  if(argc < 3 ) {
    std::cerr<<" to be used as: ./<exe file> <cfg file> <outputPlot> "<<std::endl;
    return -1;
  }

  // Setting for style                                                                                                                                                         
  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C(1)").c_str());
  }

  gStyle->SetOptStat(111110);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);

  // output folders                                                                                                                                                             
  std::string outputFileDirectory = argv[2];

  system(("mkdir -p "+outputFileDirectory).c_str());
  system(("rm -r "   +outputFileDirectory+"/*").c_str());

  // input files and chain to be analyzed 
  std::string inputFileDirectory = argv[1]; 

  TChain *inputChain = new TChain("Delphes");
  inputChain->Add((inputFileDirectory+"/*1.root").c_str());
  inputChain->Add((inputFileDirectory+"/*2.root").c_str());
  inputChain->Add((inputFileDirectory+"/*3.root").c_str());
  
  std::cout<<"number of events to analyze : "<<inputChain->GetEntries()<<std::endl;

  // threshold for jet and leptons
  float leptonPtThreshold = 0;
  if(argc > 2) leptonPtThreshold = atof(argv[3]);

  float jetPtThreshold = 0;
  if(argc > 3) jetPtThreshold = atof(argv[4]);


  // set the branch address to be used
  TClonesArray* GenMET = new TClonesArray("MissingET");
  inputChain->SetBranchAddress("GenMissingET",&GenMET);

  TClonesArray* RecoMET = new TClonesArray("MissingET");
  inputChain->SetBranchAddress("MissingET",&RecoMET);

  TClonesArray* PuppiMET = new TClonesArray("MissingET");
  inputChain->SetBranchAddress("PuppiMissingET",&PuppiMET);

  TClonesArray* GenHT = new TClonesArray("ScalarHT");
  inputChain->SetBranchAddress("GenHT",&GenHT);

  TClonesArray* RecoHT = new TClonesArray("ScalarHT");
  inputChain->SetBranchAddress("HT",&RecoHT);

  TClonesArray* PuppiHT = new TClonesArray("ScalarHT");
  inputChain->SetBranchAddress("PuppiHT",&PuppiHT);

  TClonesArray* Muons = new TClonesArray("Muon");
  inputChain->SetBranchAddress("Muon",&Muons);

  TClonesArray* Electrons = new TClonesArray("Electron");
  inputChain->SetBranchAddress("Electron",&Electrons);

  TClonesArray* Jets = new TClonesArray("Jet");
  inputChain->SetBranchAddress("JetPUID",&Jets);

  inputChain->SetBranchStatus("LHE*",0);
  inputChain->SetBranchStatus("PuppiJet*",0);
  inputChain->SetBranchStatus("LHE*",0);
  inputChain->SetBranchStatus("Track*",0);
  inputChain->SetBranchStatus("*Rho*",0);
  inputChain->SetBranchStatus("GenJet*",0);
  inputChain->SetBranchStatus("GenMissingET",1);
  inputChain->SetBranchStatus("MissingET",1);
  inputChain->SetBranchStatus("PuppiMissingET",1);
  inputChain->SetBranchStatus("GenHT",1);
  inputChain->SetBranchStatus("HT",1);
  inputChain->SetBranchStatus("PuppiHT",1);
  inputChain->SetBranchStatus("Muon",1);
  inputChain->SetBranchStatus("Electron",1);
  inputChain->SetBranchStatus("JetPUID",1);

 
  // histogram to fill on the met
  TH1F* ScalarHT_Gen     = new TH1F ("ScalarHT_Gen","",300,0,10000);
  TH1F* MissingET_Gen    = new TH1F ("MissingET_Gen","",50, 0,600);
  TH1F* MissingETPhi_Gen = new TH1F ("MissingETPhi_Gen","",50,-3.14,3.14);

  TH1F* ScalarHT_Reco     = new TH1F ("ScalarHT_Reco","",300,0,10000);
  TH1F* MissingET_Reco    = new TH1F ("MissingET_Reco","",50, 0,600);
  TH1F* MissingETPhi_Reco = new TH1F ("MissingETPhi_Reco","",50,-3.14,3.14);

  TH1F* ScalarHT_Puppi     = new TH1F ("ScalarHT_Puppi","",300,0,10000);
  TH1F* MissingET_Puppi    = new TH1F ("MissingET_Puppi","",50, 0,600);
  TH1F* MissingETPhi_Puppi = new TH1F ("MissingETPhi_Puppi","",50,-3.14,3.14);

  TH1F* ScalarHT_Resp     = new TH1F ("ScalarHT_Resp","",100,-1000,1000);
  TH1F* MissingET_Resp    = new TH1F ("MissingET_Resp","",50,-250,250);
  TH1F* MissingETPhi_Resp = new TH1F ("MissingETPhi_Resp","",50,-0.7,0.7);

  TH1F* ScalarHT_Resp_Puppi     = new TH1F ("ScalarHT_Resp_Puppi","",125,-1500,1500);
  TH1F* MissingET_Resp_Puppi    = new TH1F ("MissingET_Resp_Puppi","",50,-250,250);
  TH1F* MissingETPhi_Resp_Puppi = new TH1F ("MissingETPhi_Resp_Puppi","",50,-0.7,0.7);

  // loop on the events and fill the histos
  for(int iEntry = 0; iEntry < inputChain->GetEntries() ; iEntry++){

    if(iEntry%10000 == 0) std::cout<<"reading entry "<<iEntry<<std::endl;

    inputChain->GetEntry(iEntry);

    int nLep = 0;
    int nJet = 0;

    // require at least two leptons over a fixed PT cut, no isolation required
    for(int iMuon = 0; iMuon < Muons->GetEntries(); iMuon++){
      Muon* muon = (Muon*) Muons->At(iMuon);
      if(muon->PT >= leptonPtThreshold) nLep++; 
    }

    for(int iElectron = 0; iElectron < Electrons->GetEntries(); iElectron++){
      Electron* elec = (Electron*) Electrons->At(iElectron);
      if(elec->PT >= leptonPtThreshold) nLep++; 
    }

    if(nLep < 2) continue ;

    // require at least two jets (reco) over a fixed PT cut
    for(int iJet = 0; iJet < Jets->GetEntries(); iJet++){
      Jet* jet = (Jet*) Jets->At(iJet);
      if(jet->PT >= jetPtThreshold) nJet++; 
    }

    if(nJet < 2) continue ;

    // fill histo
    ScalarHT_Gen->Fill(dynamic_cast<ScalarHT*>(GenHT->At(0))->HT);
    ScalarHT_Reco->Fill(dynamic_cast<ScalarHT*>(RecoHT->At(0))->HT);
    ScalarHT_Puppi->Fill(dynamic_cast<ScalarHT*>(PuppiHT->At(0))->HT);

    ScalarHT_Resp->Fill(dynamic_cast<ScalarHT*>(RecoHT->At(0))->HT-dynamic_cast<ScalarHT*>(GenHT->At(0))->HT);
    ScalarHT_Resp_Puppi->Fill(dynamic_cast<ScalarHT*>(PuppiHT->At(0))->HT-dynamic_cast<ScalarHT*>(GenHT->At(0))->HT);

    MissingET_Gen->Fill(dynamic_cast<MissingET*>(GenMET->At(0))->MET);
    MissingET_Reco->Fill(dynamic_cast<MissingET*>(RecoMET->At(0))->MET);
    MissingET_Puppi->Fill(dynamic_cast<MissingET*>(PuppiMET->At(0))->MET);

    MissingET_Resp->Fill(dynamic_cast<MissingET*>(RecoMET->At(0))->MET-dynamic_cast<MissingET*>(GenMET->At(0))->MET);
    MissingET_Resp_Puppi->Fill(dynamic_cast<MissingET*>(PuppiMET->At(0))->MET-dynamic_cast<MissingET*>(GenMET->At(0))->MET);

    MissingETPhi_Gen->Fill(dynamic_cast<MissingET*>(GenMET->At(0))->Phi);
    MissingETPhi_Reco->Fill(dynamic_cast<MissingET*>(RecoMET->At(0))->Phi);
    MissingETPhi_Puppi->Fill(dynamic_cast<MissingET*>(PuppiMET->At(0))->Phi);

    MissingETPhi_Resp->Fill(dynamic_cast<MissingET*>(RecoMET->At(0))->Phi-dynamic_cast<MissingET*>(GenMET->At(0))->Phi);
    MissingETPhi_Resp_Puppi->Fill(dynamic_cast<MissingET*>(PuppiMET->At(0))->Phi-dynamic_cast<MissingET*>(GenMET->At(0))->Phi);

  }

  
  ///////////////////////////////////////////////                                                                                                                       
  // Plot                                                                                                                                                                       
  ///////////////////////////////////////////////                                                                                                                            

  // make the canvas and basic banners                                                                                                                                       
  TCanvas *cCanvas = new TCanvas("cCanvas","",180,52,550,550);
  cCanvas->SetTicks();
  cCanvas->SetFillColor(0);
  cCanvas->SetBorderMode(0);
  cCanvas->SetBorderSize(2);
  cCanvas->SetTickx(1);
  cCanvas->SetTicky(1);
  cCanvas->SetRightMargin(0.05);
  cCanvas->SetBottomMargin(0.12);
  cCanvas->SetFrameBorderMode(0);

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"Delphes");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.286,0.92,"Simulation Preliminary");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  TLegend* legend = new TLegend(0.16,0.78,0.36,0.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.031);
  legend->SetTextFont(42);

  // scalar HT plots
  ScalarHT_Reco->SetLineWidth(2);
  ScalarHT_Reco->SetLineColor(kBlue);
  ScalarHT_Reco->GetXaxis()->SetTitle("HT (GeV)");
  ScalarHT_Reco->GetYaxis()->SetTitle("Entries");

  ScalarHT_Reco->Draw("hist");

  gPad->Update();
  TPaveStats *tps1 = (TPaveStats*) ScalarHT_Reco->FindObject("stats");
  tps1->SetTextColor(kBlue);
  tps1->SetLineColor(kBlue);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  double X1 = tps1->GetX1NDC();
  double Y1 = tps1->GetY1NDC();
  double X2 = tps1->GetX2NDC();
  double Y2 = tps1->GetY2NDC();
  
  ScalarHT_Puppi->SetLineWidth(2);
  ScalarHT_Puppi->SetLineColor(kRed);
  ScalarHT_Puppi->Draw("hist");
  
  gPad->Update();
  TPaveStats *tps2 = (TPaveStats*) ScalarHT_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  ScalarHT_Gen->SetLineWidth(2);
  ScalarHT_Gen->SetLineColor(kBlack);
  ScalarHT_Gen->Draw("hist");
  
  gPad->Update();
  TPaveStats *tps3 = (TPaveStats*) ScalarHT_Gen->FindObject("stats");
  tps3->SetTextColor(kBlack);
  tps3->SetLineColor(kBlack);
  tps3->SetX1NDC(X1);
  tps3->SetX2NDC(X2);
  tps3->SetY1NDC(Y1-(Y2-Y1));
  tps3->SetY2NDC(Y1);

  ScalarHT_Gen->GetYaxis()->SetRangeUser(0.001,std::max(ScalarHT_Gen->GetMaximum(),std::max(ScalarHT_Reco->GetMaximum(),ScalarHT_Puppi->GetMaximum()))*1.25);

  ScalarHT_Gen->Draw("hist");
  ScalarHT_Reco->Draw("hist same");
  ScalarHT_Puppi->Draw("hist same");


  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);
  tps3->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");
  tps3->Draw("same");

  legend->AddEntry(ScalarHT_Reco,"Reco HT","l");
  legend->AddEntry(ScalarHT_Puppi,"Puppi HT","l");
  legend->AddEntry(ScalarHT_Gen,"Gen HT","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();
  
  // missing ET plots
  MissingET_Reco->SetLineWidth(2);
  MissingET_Reco->SetLineColor(kBlue);
  MissingET_Reco->GetXaxis()->SetTitle("E_{T}^{Miss} (GeV)");
  MissingET_Reco->GetYaxis()->SetTitle("Entries");

  MissingET_Reco->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) MissingET_Reco->FindObject("stats");
  tps1->SetTextColor(kBlue);
  tps1->SetLineColor(kBlue);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();

  MissingET_Puppi->SetLineWidth(2);
  MissingET_Puppi->SetLineColor(kRed);
  MissingET_Puppi->Draw("hist");

  gPad->Update();
  tps2 = (TPaveStats*) MissingET_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  MissingET_Gen->SetLineWidth(2);
  MissingET_Gen->SetLineColor(kBlack);
  MissingET_Gen->Draw("hist");

  gPad->Update();
  tps3 = (TPaveStats*) MissingET_Gen->FindObject("stats");
  tps3->SetTextColor(kBlack);
  tps3->SetLineColor(kBlack);
  tps3->SetX1NDC(X1);
  tps3->SetX2NDC(X2);
  tps3->SetY1NDC(Y1-(Y2-Y1));
  tps3->SetY2NDC(Y1);

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);
  tps3->SetFillStyle(0);

  MissingET_Gen->GetYaxis()->SetRangeUser(0.001,std::max(MissingET_Gen->GetMaximum(),std::max(MissingET_Reco->GetMaximum(),MissingET_Puppi->GetMaximum()))*1.25);

  MissingET_Gen->Draw("hist");
  MissingET_Reco->Draw("hist same");
  MissingET_Puppi->Draw("hist same");

  tps1->Draw("same");
  tps2->Draw("same");
  tps3->Draw("same");

  legend->AddEntry(MissingET_Reco,"Reco MET","l");
  legend->AddEntry(MissingET_Puppi,"Puppi MET","l");
  legend->AddEntry(MissingET_Gen,"Gen MET","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET.root").c_str(),"root");
  cCanvas->SetLogy();
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  // missing ET Phi plots
  MissingETPhi_Reco->SetLineWidth(2);
  MissingETPhi_Reco->SetLineColor(kBlue);
  MissingETPhi_Reco->GetXaxis()->SetTitle("#phi");
  MissingETPhi_Reco->GetYaxis()->SetTitle("Entries");

  MissingETPhi_Reco->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) MissingETPhi_Reco->FindObject("stats");
  tps1->SetTextColor(kBlue);
  tps1->SetLineColor(kBlue);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();

  MissingETPhi_Puppi->SetLineWidth(2);
  MissingETPhi_Puppi->SetLineColor(kRed);
  MissingETPhi_Puppi->Draw("hist");

  gPad->Update();
  tps2 = (TPaveStats*) MissingETPhi_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  MissingETPhi_Gen->SetLineWidth(2);
  MissingETPhi_Gen->SetLineColor(kBlack);
  MissingETPhi_Gen->Draw("hist");

  gPad->Update();
  tps3 = (TPaveStats*) MissingETPhi_Gen->FindObject("stats");
  tps3->SetTextColor(kBlack);
  tps3->SetLineColor(kBlack);
  tps3->SetX1NDC(X1);
  tps3->SetX2NDC(X2);
  tps3->SetY1NDC(Y1-(Y2-Y1));
  tps3->SetY2NDC(Y1);

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);
  tps3->SetFillStyle(0);

  MissingETPhi_Gen->GetYaxis()->SetRangeUser(0.001,std::max(MissingETPhi_Gen->GetMaximum(),std::max(MissingETPhi_Reco->GetMaximum(),MissingETPhi_Puppi->GetMaximum()))*1.25);


  MissingETPhi_Gen->Draw("hist");
  MissingETPhi_Reco->Draw("hist same");
  MissingETPhi_Puppi->Draw("hist same");

  tps1->Draw("same");
  tps2->Draw("same");
  tps3->Draw("same");

  legend->AddEntry(MissingETPhi_Reco,"Reco MET","l");
  legend->AddEntry(MissingETPhi_Puppi,"Puppi MET","l");
  legend->AddEntry(MissingETPhi_Gen,"Gen MET","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi.root").c_str(),"root");
  cCanvas->SetLogy();
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  // Response plot
   
  ScalarHT_Resp->SetLineWidth(2);
  ScalarHT_Resp->SetLineColor(kBlack);
  ScalarHT_Resp->GetXaxis()->SetTitle("HT-HT^{Gen} (GeV)");
  ScalarHT_Resp->GetYaxis()->SetTitle("Entries");

  ScalarHT_Resp->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ScalarHT_Resp->FindObject("stats");
  tps1->SetTextColor(kBlack);
  tps1->SetLineColor(kBlack);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();

  ScalarHT_Resp_Puppi->SetLineWidth(2);
  ScalarHT_Resp_Puppi->SetLineColor(kRed);
  ScalarHT_Resp_Puppi->Draw("hist");

  gPad->Update();
  tps2 = (TPaveStats*) ScalarHT_Resp_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  ScalarHT_Resp->GetYaxis()->SetRangeUser(0.001,std::max(ScalarHT_Resp->GetMaximum(),ScalarHT_Resp_Puppi->GetMaximum())*1.25);

  ScalarHT_Resp->Draw("hist");
  ScalarHT_Resp_Puppi->Draw("hist same");

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ScalarHT_Resp,"Reco HT","l");
  legend->AddEntry(ScalarHT_Resp_Puppi,"Puppi HT","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp.root").c_str(),"root");
  cCanvas->SetLogy();
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/ScalarHT_Resp_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  /////
  MissingET_Resp->SetLineWidth(2);
  MissingET_Resp->SetLineColor(kBlack);
  MissingET_Resp->GetXaxis()->SetTitle("MET-MET^{Gen} (GeV)");
  MissingET_Resp->GetYaxis()->SetTitle("Entries");

  MissingET_Resp->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) MissingET_Resp->FindObject("stats");
  tps1->SetTextColor(kBlack);
  tps1->SetLineColor(kBlack);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();

  MissingET_Resp_Puppi->SetLineWidth(2);
  MissingET_Resp_Puppi->SetLineColor(kRed);
  MissingET_Resp_Puppi->Draw("hist");

  gPad->Update();
  tps2 = (TPaveStats*) MissingET_Resp_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  MissingET_Resp->GetYaxis()->SetRangeUser(0.001,std::max(MissingET_Resp->GetMaximum(),MissingET_Resp_Puppi->GetMaximum())*1.25);

  MissingET_Resp->Draw("hist");
  MissingET_Resp_Puppi->Draw("hist same");

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(MissingET_Resp,"Reco MET","l");
  legend->AddEntry(MissingET_Resp_Puppi,"Puppi MET","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp.root").c_str(),"root");
  cCanvas->SetLogy();
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingET_Resp_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update(); 

  legend->Clear();

  ///
  MissingETPhi_Resp->SetLineWidth(2);
  MissingETPhi_Resp->SetLineColor(kBlack);
  MissingETPhi_Resp->GetXaxis()->SetTitle("#phi(MET)-#phi(MET^{Gen})");
  MissingETPhi_Resp->GetYaxis()->SetTitle("Entries");

  MissingETPhi_Resp->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) MissingETPhi_Resp->FindObject("stats");
  tps1->SetTextColor(kBlack);
  tps1->SetLineColor(kBlack);
  tps1->SetX1NDC(0.72);
  tps1->SetY1NDC(0.68);
  tps1->SetX2NDC(0.93);
  tps1->SetY2NDC(0.89);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();

  MissingETPhi_Resp_Puppi->SetLineWidth(2);
  MissingETPhi_Resp_Puppi->SetLineColor(kRed);
  MissingETPhi_Resp_Puppi->Draw("hist");

  gPad->Update();
  tps2 = (TPaveStats*) MissingETPhi_Resp_Puppi->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  MissingETPhi_Resp->GetYaxis()->SetRangeUser(0.001,std::max(MissingETPhi_Resp->GetMaximum(),MissingETPhi_Resp_Puppi->GetMaximum())*1.25);

  MissingETPhi_Resp->Draw("hist");
  MissingETPhi_Resp_Puppi->Draw("hist same");

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(MissingETPhi_Resp,"Reco MET","l");
  legend->AddEntry(MissingETPhi_Resp_Puppi,"Puppi MET","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp.root").c_str(),"root");
  cCanvas->SetLogy();
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp_log.png").c_str(),"png");
  cCanvas->SaveAs(std::string(outputFileDirectory+"/MissingETPhi_Resp_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  return 0;

}
