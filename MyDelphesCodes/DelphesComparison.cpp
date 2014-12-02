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

#include "classes/DelphesClasses.h"

// g++ -Wall -o DelphesComparison `root-config --glibs --libs --cflags` -lTreePlayer DelphesComparison.cpp
typedef std::pair<TH1F*,TH1F*> histoPair ;
typedef std::pair<std::string, histoPair > mapElement;

float deltaPhi (float a, float b){
  if(fabs(a-b) > TMath::Pi()) return 2*TMath::Pi()-fabs(a-b);
  else return fabs(a-b);

}

void getHistogram(std::map<std::string,histoPair> & map){

  map.insert(mapElement("NPU",histoPair(new TH1F("NPU_1","",100,0,100),new TH1F("NPU_2","",100,0,100))));  
  /*
  map.insert(mapElement("Rho",histoPair(new TH1F("Rho_1","",30,0,100),new TH1F("Rho_2","",30,0,100))));
  map.insert(mapElement("PuppiRho",histoPair(new TH1F("PuppiRho_1","",30,0,100),new TH1F("PuppiRho_2","",30,0,100))));

  map.insert(mapElement("puppiParticlePT",histoPair(new TH1F("puppiParticlePT_1","",30,0,100),new TH1F("puppiParticlePT_2","",30,0,100))));
  map.insert(mapElement("puppiParticleEta",histoPair(new TH1F("puppiParticleEta_1","",25,-5,5),new TH1F("puppiParticleEta_2","",25,-5,5))));

  map.insert(mapElement("JetPtLead",  histoPair(new TH1F("JetPtLead_1","",50,10,1000),new TH1F("JetPtLead_2","",50,0,100))));
  map.insert(mapElement("JetEtaLead", histoPair(new TH1F("JetEtaLead_1","",50,-5,5),  new TH1F("JetEtaLead_2","",50,-5,5))));
  map.insert(mapElement("JetMassLead",histoPair(new TH1F("JetMassLead_1","",30,0,150),new TH1F("JetMassLead_2","",30,0,150))));

  map.insert(mapElement("JetPtSecond",  histoPair(new TH1F("JetPtSecond_1","",50,10,1000),new TH1F("JetPtSecond_2","",50,0,100))));
  map.insert(mapElement("JetEtaSecond", histoPair(new TH1F("JetEtaSecond_1","",50,-5,5),  new TH1F("JetEtaSecond_2","",50,-5,5))));
  map.insert(mapElement("JetMassSecond",histoPair(new TH1F("JetMassSecond_1","",30,0,150),new TH1F("JetMassSecond_2","",30,0,150))));

  map.insert(mapElement("JetPtThird",  histoPair(new TH1F("JetPtThird_1","",50,10,1000),new TH1F("JetPtThird_2","",50,0,100))));
  map.insert(mapElement("JetEtaThird", histoPair(new TH1F("JetEtaThird_1","",50,-5,5),  new TH1F("JetEtaThird_2","",50,-5,5))));
  map.insert(mapElement("JetMassThird",histoPair(new TH1F("JetMassThird_1","",30,0,150),new TH1F("JetMassThird_2","",30,0,150))));

  map.insert(mapElement("PuppiJetPtLead",  histoPair(new TH1F("PuppiJetPtLead_1","",50,10,1000),new TH1F("PuppiJetPtLead_2","",50,0,100))));
  map.insert(mapElement("PuppiJetEtaLead", histoPair(new TH1F("PuppiJetEtaLead_1","",50,-5,5),  new TH1F("PuppiJetEtaLead_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassLead",histoPair(new TH1F("PuppiJetMassLead_1","",30,0,150),new TH1F("PuppiJetMassLead_2","",30,0,150))));

  map.insert(mapElement("PuppiJetPtSecond",  histoPair(new TH1F("PuppiJetPtSecond_1","",50,10,1000),new TH1F("PuppiJetPtSecond_2","",50,0,100))));
  map.insert(mapElement("PuppiJetEtaSecond", histoPair(new TH1F("PuppiJetEtaSecond_1","",50,-5,5),  new TH1F("PuppiJetEtaSecond_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassSecond",histoPair(new TH1F("PuppiJetMassSecond_1","",30,0,150),new TH1F("PuppiJetMassSecond_2","",30,0,150))));

  map.insert(mapElement("PuppiJetPtThird",  histoPair(new TH1F("PuppiJetPtThird_1","",50,10,1000),new TH1F("PuppiJetPtThird_2","",50,0,100))));
  map.insert(mapElement("PuppiJetEtaThird", histoPair(new TH1F("PuppiJetEtaThird_1","",50,-5,5),  new TH1F("PuppiJetEtaThird_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassThird",histoPair(new TH1F("PuppiJetMassThird_1","",30,0,150),new TH1F("PuppiJetMassThird_2","",30,0,150))));

  map.insert(mapElement("GenMissingET",histoPair(new TH1F("GenMissingET_1","",30,0,600),new TH1F("GenMissingET_2","",30,0,600))));
  map.insert(mapElement("RecoMissingET",histoPair(new TH1F("RecoMissingET_1","",30,0,600),new TH1F("RecoMissingET_2","",30,0,600))));
  map.insert(mapElement("PuppiMissingET",histoPair(new TH1F("PuppiMissingET_1","",30,0,600),new TH1F("PuppiMissingET_2","",30,0,600))));

  map.insert(mapElement("MuonPt",histoPair(new TH1F("MuonPt_1","",50,10,1000),new TH1F("MuonPt_2","",50,0,100))));
  map.insert(mapElement("MuonEta",histoPair(new TH1F("MuonEta_1","",50,-5,5),new TH1F("MuonPt_2","",50,-5,5))));

  map.insert(mapElement("ElectronPt",histoPair(new TH1F("ElectronPt_1","",50,10,1000),new TH1F("ElectronPt_2","",50,0,100))));
  map.insert(mapElement("ElectronEta",histoPair(new TH1F("ElectronEta_1","",50,-5,5),new TH1F("ElectronPt_2","",50,-5,5))));  
  */
}

void getResponseHistogram( std::map<std::string,histoPair> & map) {
  /*
  map.insert(mapElement("MetResp",histoPair(new TH1F("MetResp_1","",50,-100,100),new TH1F("MetPtResp_2","",50,-100,100))));
  map.insert(mapElement("puppiMetResp",histoPair(new TH1F("puppiMetResp_1","",50,-100,100),new TH1F("puppiMetPtResp_2","",50,-100,100))));

  map.insert(mapElement("JetPtRespLead",histoPair(new TH1F("JetPtRespLead_1","",50,-100,100),new TH1F("JetPtRespLead_2","",50,-100,100))));
  map.insert(mapElement("JetEtaRespLead",histoPair(new TH1F("JetEtaRespLead_1","",25,-5,5),new TH1F("JetEtaRespLead_2","",25,-5,5))));
  map.insert(mapElement("JetMassRespLead",histoPair(new TH1F("JetMassRespLead_1","",25,-50,50),new TH1F("JetMassRespLead_2","",25,-50,50))));

  map.insert(mapElement("JetPtRespSecond",histoPair(new TH1F("JetPtRespSecond_1","",50,-100,100),new TH1F("JetPtRespSecond_2","",50,-100,100))));
  map.insert(mapElement("JetEtaRespSecond",histoPair(new TH1F("JetEtaRespSecond_1","",25,-5,5),new TH1F("JetEtaRespSecond_2","",25,-5,5))));
  map.insert(mapElement("JetMassRespSecond",histoPair(new TH1F("JetMassRespSecond_1","",25,-50,50),new TH1F("JetMassRespSecond_2","",25,-50,50))));

  map.insert(mapElement("JetPtRespThird",histoPair(new TH1F("JetPtRespThird_1","",50,-100,100),new TH1F("JetPtRespThird_2","",50,-100,100))));
  map.insert(mapElement("JetEtaRespThird",histoPair(new TH1F("JetEtaRespThird_1","",25,-5,5),new TH1F("JetEtaRespThird_2","",25,-5,5))));
  map.insert(mapElement("JetMassRespThird",histoPair(new TH1F("JetMassRespThird_1","",25,-50,50),new TH1F("JetMassRespThird_2","",25,-50,50))));

  map.insert(mapElement("PuppiJetPtRespLead",histoPair(new TH1F("PuppiJetPtRespLead_1","",50,-100,100),new TH1F("PuppiJetPtRespLead_2","",50,-100,100))));
  map.insert(mapElement("PuppiJetEtaRespLead",histoPair(new TH1F("PuppiJetEtaRespLead_1","",25,-5,5),new TH1F("PuppiJetEtaRespLead_2","",25,-5,5))));
  map.insert(mapElement("PuppiJetMassRespLead",histoPair(new TH1F("PuppiJetMassRespLead_1","",25,-50,50),new TH1F("PuppiJetMassRespLead_2","",25,-50,50))));

  map.insert(mapElement("PuppiJetPtRespSecond",histoPair(new TH1F("PuppiJetPtRespSecond_1","",50,-100,100),new TH1F("PuppiJetPtRespSecond_2","",50,-100,100))));
  map.insert(mapElement("PuppiJetEtaRespSecond",histoPair(new TH1F("PuppiJetEtaRespSecond_1","",25,-5,5),new TH1F("PuppiJetEtaRespSecond_2","",25,-5,5))));
  map.insert(mapElement("PuppiJetMassRespSecond",histoPair(new TH1F("PuppiJetMassRespSecond_1","",25,-50,50),new TH1F("PuppiJetMassRespSecond_2","",25,-50,50))));

  map.insert(mapElement("PuppiJetPtRespThird",histoPair(new TH1F("PuppiJetPtRespThird_1","",50,-100,100),new TH1F("PuppiJetPtRespThird_2","",50,-100,100))));
  map.insert(mapElement("PuppiJetEtaRespThird",histoPair(new TH1F("PuppiJetEtaRespThird_1","",25,-5,5),new TH1F("PuppiJetEtaRespThird_2","",25,-5,5))));
  map.insert(mapElement("PuppiJetMassRespThird",histoPair(new TH1F("PuppiJetMassRespThird_1","",25,-50,50),new TH1F("PuppiJetMassRespThird_2","",25,-50,50))));
  */
}


/////////////////////////////////////////////////////

int main (int argc, char** argv){

  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  }
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);

  std::string outputFileDirectory = "MyDelphesCodes/outputPlots";
  system(("mkdir -p "+outputFileDirectory).c_str());
  system(("rm -r "+outputFileDirectory+"/*").c_str());

  std::vector<std::string> inputFileList_1 ;
  std::vector<std::string> inputFileList_2 ;

  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_0.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_2.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_5.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_6.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_7.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_10.root");
  inputFileList_1.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev/outputtree_14.root");

  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_0.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_2.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_5.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_6.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_7.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_10.root");
  inputFileList_2.push_back("eos/cms/store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_SethZenz/outputtree_14.root");

  TChain* chain_1 = new TChain("Delphes");
  TChain* chain_2 = new TChain("Delphes");


  for( size_t iVec = 0; iVec < inputFileList_1.size() ; iVec++){
    chain_1->Add(inputFileList_1.at(iVec).c_str());
  }
  
  for( size_t iVec = 0; iVec < inputFileList_2.size() ; iVec++){
    chain_2->Add(inputFileList_2.at(iVec).c_str());
  }
  
  // Take input informations
  TClonesArray* NPU_1 = new TClonesArray("ScalarHT");
  TClonesArray* NPU_2 = new TClonesArray("ScalarHT");
  chain_1->SetBranchAddress("NPU",&NPU_1);    
  chain_2->SetBranchAddress("NPU",&NPU_2);    

  TClonesArray* Rho_1 = new TClonesArray("Rho");
  TClonesArray* Rho_2 = new TClonesArray("Rho");
  chain_1->SetBranchAddress("RhoKt4",&Rho_1);    
  chain_2->SetBranchAddress("RhoKt4",&Rho_2);    

  TClonesArray* GenJets_1 = new TClonesArray("Jet");
  TClonesArray* GenJets_2 = new TClonesArray("Jet");
  chain_1->SetBranchAddress("GenJet",&GenJets_1);    
  chain_2->SetBranchAddress("GenJet",&GenJets_2);    

  TClonesArray* Jets_1 = new TClonesArray("Jet");
  TClonesArray* Jets_2 = new TClonesArray("Jet");
  chain_1->SetBranchAddress("JetPUID",&Jets_1);    
  chain_2->SetBranchAddress("Jet",&Jets_2);    

  TClonesArray* puppiParticle_1 = new TClonesArray("GenParticle");
  TClonesArray* puppiParticle_2 = new TClonesArray("GenParticle");
  chain_1->SetBranchAddress("puppiParticles",&puppiParticle_1);    
  chain_2->SetBranchAddress("puppiParticles",&puppiParticle_2);    

  TClonesArray* puppiRho_1 = new TClonesArray("Rho");
  TClonesArray* puppiRho_2 = new TClonesArray("Rho");
  chain_1->SetBranchAddress("PuppiRhoKt4",&puppiRho_1);    
  chain_2->SetBranchAddress("PuppiRhoKt4",&puppiRho_2);    

  TClonesArray* puppiJets_1 = new TClonesArray("Jet");
  TClonesArray* puppiJets_2 = new TClonesArray("Jet");
  chain_1->SetBranchAddress("PuppiJetPUID",&puppiJets_1);    
  chain_2->SetBranchAddress("PuppiJet",&puppiJets_2);    

  TClonesArray* GenMET_1 = new TClonesArray("MissingET");
  TClonesArray* GenMET_2 = new TClonesArray("MissingET");
  chain_1->SetBranchAddress("GenMissingET",&GenMET_1);    
  chain_2->SetBranchAddress("GenMissingET",&GenMET_2);    

  TClonesArray* RecoMET_1 = new TClonesArray("MissingET");
  TClonesArray* RecoMET_2 = new TClonesArray("MissingET");
  chain_1->SetBranchAddress("MissingET",&RecoMET_1);    
  chain_2->SetBranchAddress("MissingET",&RecoMET_2);    

  TClonesArray* PuppiMET_1 = new TClonesArray("MissingET");
  TClonesArray* PuppiMET_2 = new TClonesArray("MissingET");
  chain_1->SetBranchAddress("PuppiMissingET",&PuppiMET_1);    
  chain_2->SetBranchAddress("PuppiMissingET",&PuppiMET_2);    

  TClonesArray* Muon_1 = new TClonesArray("Muon");
  TClonesArray* Muon_2 = new TClonesArray("Muon");
  chain_1->SetBranchAddress("Muon",&Muon_1);    
  chain_2->SetBranchAddress("Muon",&Muon_2);    

  TClonesArray* Electron_1 = new TClonesArray("Electron");
  TClonesArray* Electron_2 = new TClonesArray("Electron");
  chain_1->SetBranchAddress("Electron",&Electron_1);    
  chain_2->SetBranchAddress("Electron",&Electron_2);    

  // Loop on the events
  int maxEvents = 0;
  // take the max number of events
  (chain_1->GetEntries() < chain_2->GetEntries()) ? maxEvents = chain_1->GetEntries() : maxEvents = chain_2->GetEntries();
  // List of histograms
  std::map<std::string,histoPair> histogramSingleVariables ;
  getHistogram(histogramSingleVariables) ;
  std::map<std::string,histoPair> histogramResponse ;
  getResponseHistogram(histogramResponse);

  // set errors
  for(std::map<std::string,histoPair>::iterator itMap = histogramSingleVariables.begin(); itMap != histogramSingleVariables.end(); itMap++){
    itMap->second.first->Sumw2();
    itMap->second.second->Sumw2();
  }

  for(std::map<std::string,histoPair>::iterator itMap = histogramResponse.begin(); itMap != histogramResponse.end(); itMap++){
    itMap->second.first->Sumw2();
    itMap->second.second->Sumw2();
  }
  
  std::cout<<" Max Entries "<<maxEvents<<std::endl;  
  for(int iEntry = 0; iEntry < maxEvents/15 ; iEntry++){

    if(iEntry%1000 == 0) std::cout<<" reading entry "<<iEntry<<std::endl;

    // Read the entry
    chain_1->GetEntry(iEntry);
    chain_2->GetEntry(iEntry);
    /*   
    // fill PU
    for( int i = 0; i < NPU_1->GetEntries() and i < NPU_2->GetEntries(); i++){
      histogramSingleVariables["NPU"].first->Fill(dynamic_cast<ScalarHT*>(NPU_1->At(i))->HT);
      histogramSingleVariables["NPU"].second->Fill(dynamic_cast<ScalarHT*>(NPU_2->At(i))->HT);
    }
    /*
    // Rho Kt4
    for( int i = 0; i < Rho_1->GetEntries() and i < Rho_2->GetEntries() ; i++){
      histogramSingleVariables["Rho"].first->Fill(dynamic_cast<Rho*>(Rho_1->At(i))->Rho);
      histogramSingleVariables["Rho"].second->Fill(dynamic_cast<Rho*>(Rho_2->At(i))->Rho);
    }

    // Puppi Rho Kt4
    for( int i = 0; i < puppiRho_1->GetEntries() and i < puppiRho_2->GetEntries(); i++){
      histogramSingleVariables["PuppiRho"].first->Fill(dynamic_cast<Rho*>(puppiRho_1->At(i))->Rho);
      histogramSingleVariables["PuppiRho"].second->Fill(dynamic_cast<Rho*>(puppiRho_2->At(i))->Rho);
    }

    //Muon 
    for( int i = 0; i < Muon_1->GetEntries() and i < Muon_2->GetEntries() ; i++){
      histogramSingleVariables["MuonPt"].first->Fill(dynamic_cast<Muon*>(Muon_1->At(i))->PT);
      histogramSingleVariables["MuonPt"].second->Fill(dynamic_cast<Muon*>(Muon_2->At(i))->PT);

      histogramSingleVariables["MuonEta"].first->Fill(dynamic_cast<Muon*>(Muon_1->At(i))->Eta);
      histogramSingleVariables["MuonEta"].second->Fill(dynamic_cast<Muon*>(Muon_2->At(i))->Eta);
    }

    //Electron 
    for( int i = 0; i< Electron_1->GetEntries() and i < Electron_2->GetEntries(); i++){
      histogramSingleVariables["ElectronPt"].first->Fill(dynamic_cast<Electron*>(Electron_1->At(i))->PT);
      histogramSingleVariables["ElectronPt"].second->Fill(dynamic_cast<Electron*>(Electron_2->At(i))->PT);

      histogramSingleVariables["ElectronEta"].first->Fill(dynamic_cast<Electron*>(Electron_1->At(i))->Eta);
      histogramSingleVariables["ElectronEta"].second->Fill(dynamic_cast<Electron*>(Electron_2->At(i))->Eta);
    }

    // GEN MET
    for(int i = 0; i < GenMET_1->GetEntries() and i < GenMET_2->GetEntries(); i++){
      histogramSingleVariables["GenMissingET"].first->Fill(dynamic_cast<MissingET*>(GenMET_1->At(i))->MET);
      histogramSingleVariables["GenMissingET"].second->Fill(dynamic_cast<MissingET*>(GenMET_2->At(i))->MET);
    }

    // RECO MET
    for(int i = 0; i < RecoMET_1->GetEntries() and i < RecoMET_2->GetEntries(); i++){
      histogramSingleVariables["RecoMissingET"].first->Fill(dynamic_cast<MissingET*>(RecoMET_1->At(i))->MET);
      histogramSingleVariables["RecoMissingET"].second->Fill(dynamic_cast<MissingET*>(RecoMET_2->At(i))->MET);
    }

    // Puppi MET
    for(int i = 0; i < PuppiMET_1->GetEntries() and i < PuppiMET_2->GetEntries(); i++){
      histogramSingleVariables["PuppiMissingET"].first->Fill(dynamic_cast<MissingET*>(PuppiMET_1->At(i))->MET);
      histogramSingleVariables["PuppiMissingET"].second->Fill(dynamic_cast<MissingET*>(PuppiMET_2->At(i))->MET);
    }

    // puppi particle
    for(int i = 0; i < puppiParticle_1->GetEntries() ; i++){
      histogramSingleVariables["puppiParticlesPT"].first->Fill(dynamic_cast<GenParticle*>(puppiParticle_1->At(i))->PT);
      histogramSingleVariables["puppiParticlesEta"].first->Fill(dynamic_cast<GenParticle*>(puppiParticle_1->At(i))->Eta);
    }

    for(int i = 0; i < puppiParticle_2->GetEntries() ; i++){
      histogramSingleVariables["puppiParticlesPT"].second->Fill(dynamic_cast<GenParticle*>(puppiParticle_2->At(i))->PT);
      histogramSingleVariables["puppiParticlesEta"].second->Fill(dynamic_cast<GenParticle*>(puppiParticle_2->At(i))->Eta);
    }

    // JET 1 
    for(int i = 0; i < Jets_1->GetEntries() and i < 2; i++){
      if(i==0){
       histogramSingleVariables["JetPtLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT);
       histogramSingleVariables["JetEtaLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta);
       histogramSingleVariables["JetMassLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["JetPtSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT);
       histogramSingleVariables["JetEtaSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta);
       histogramSingleVariables["JetMassSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["JetPtThrid"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass);
      }
      else break;
    }

    // JET 2
    for(int i = 0; i < Jets_2->GetEntries() and i < 2; i++){
      if(i==0){
       histogramSingleVariables["JetPtLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT);
       histogramSingleVariables["JetEtaLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta);
       histogramSingleVariables["JetMassLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["JetPtSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT);
       histogramSingleVariables["JetEtaSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta);
       histogramSingleVariables["JetMassSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["JetPtThrid"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass);
      }
      else break;
    }


    // PUPPI JET 1 
    for(int i = 0; i < puppiJets_1->GetEntries() and i < 2; i++){
      if(i==0){
       histogramSingleVariables["JetPtLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["JetEtaLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["JetMassLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["JetPtSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["JetEtaSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["JetMassSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["JetPtThrid"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
      else break;
    }

    // PUPPI JET 2
    for(int i = 0; i < puppiJets_2->GetEntries() and i < 2; i++){
      if(i==0){
       histogramSingleVariables["JetPtLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["JetEtaLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["JetMassLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["JetPtSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["JetEtaSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["JetMassSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["JetPtThrid"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
      else break;
    }
     
    ///////////////////////////////////
    // RESPONSE
    /////////////////////////////////////

    //  MET
    for(int i = 0; i < RecoMET_1->GetEntries() and i < GenMET_1->GetEntries(); i++){
      histogramResponse["MetResp"].first->Fill(dynamic_cast<MissingET*>(RecoMET_1->At(i))->MET-dynamic_cast<MissingET*>(GenMET_1->At(i))->MET);
    }
    for(int i = 0; i < RecoMET_2->GetEntries() and i < GenMET_2->GetEntries(); i++){
      histogramResponse["MetResp"].second->Fill(dynamic_cast<MissingET*>(RecoMET_2->At(i))->MET-dynamic_cast<MissingET*>(GenMET_2->At(i))->MET);
    }

    for(int i = 0; i < PuppiMET_1->GetEntries() and i < GenMET_1->GetEntries(); i++){
      histogramResponse["puppiMetResp"].first->Fill(dynamic_cast<MissingET*>(PuppiMET_1->At(i))->MET-dynamic_cast<MissingET*>(GenMET_1->At(i))->MET);
    }
    for(int i = 0; i < PuppiMET_2->GetEntries() and i < GenMET_2->GetEntries(); i++){
      histogramResponse["puppiMetResp"].second->Fill(dynamic_cast<MissingET*>(PuppiMET_2->At(i))->MET-dynamic_cast<MissingET*>(GenMET_2->At(i))->MET);
    }

    for(int i = 0; i < Jets_1->GetEntries() and i < 2 ; i++){
      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_1->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(dynamic_cast<Jet*>(Jets_1->At(i))->Eta+dynamic_cast<Jet*>(GenJets_1->At(j))->Eta),2)+pow(deltaPhi(dynamic_cast<Jet*>(Jets_1->At(i))->Phi,dynamic_cast<Jet*>(GenJets_1->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["JetPtRespLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespLead"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["JetPtRespSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespSecond"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["JetPtRespThrid"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
       else break;
      }
    }

    for(int i = 0; i < Jets_2->GetEntries() and i < 2 ; i++){
      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_2->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(dynamic_cast<Jet*>(Jets_2->At(i))->Eta+dynamic_cast<Jet*>(GenJets_2->At(j))->Eta),2)+pow(deltaPhi(dynamic_cast<Jet*>(Jets_2->At(i))->Phi,dynamic_cast<Jet*>(GenJets_2->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["JetPtRespLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespLead"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["JetPtRespSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespSecond"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["JetPtRespThrid"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
       else break;
      }
    }

    /////////// PUPPI JET 

    for(int i = 0; i < puppiJets_1->GetEntries() and i < 2 ; i++){
      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_1->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta+dynamic_cast<Jet*>(GenJets_1->At(j))->Eta),2)+pow(deltaPhi(dynamic_cast<Jet*>(puppiJets_1->At(i))->Phi,dynamic_cast<Jet*>(GenJets_1->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["JetPtRespLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["JetPtRespSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["JetPtRespThrid"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
       else break;
      }
    }

    for(int i = 0; i < puppiJets_2->GetEntries() and i < 2 ; i++){
      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_2->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta+dynamic_cast<Jet*>(GenJets_2->At(j))->Eta),2)+pow(deltaPhi(dynamic_cast<Jet*>(puppiJets_2->At(i))->Phi,dynamic_cast<Jet*>(GenJets_2->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["JetPtRespLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["JetPtRespSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["JetPtRespThrid"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
       else break;
      }
    }      
    */
  
  }
  ///////////////////////////////////////////////
  // Plot
  ///////////////////////////////////////////////
  // make the plot vs PT
  /*
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
  
  TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);   

  for(std::map<std::string,histoPair>::const_iterator itMap = histogramSingleVariables.begin(); itMap != histogramSingleVariables.end(); itMap++){
    cCanvas->cd();
    legend->Clear();
    itMap->second.first->Scale(1/itMap->second.first->Integral());   
    itMap->second.second->Scale(1/itMap->second.second->Integral());   

    itMap->second.first->SetLineColor(kBlack);
    itMap->second.first->SetMarkerColor(kBlack);
    itMap->second.first->SetLineWidth(2);
    itMap->second.first->GetXaxis()->SetTitle(itMap->second.first->GetName());
    itMap->second.first->GetYaxis()->SetTitle("a.u.");
    itMap->second.first->SetMarkerStyle(20);
    itMap->second.first->Draw("p");
    itMap->second.second->SetLineColor(kRed);
    itMap->second.second->SetMarkerColor(kRed);
    itMap->second.second->SetLineWidth(2);
    itMap->second.second->SetMarkerStyle(22);
    itMap->second.second->Draw("psame");

    tex->Draw("same");
    tex2->Draw("same");
    tex3->Draw("same");

    legend->AddEntry(itMap->second.first,"New Delphes","pl");
    legend->AddEntry(itMap->second.second,"Old Delphes","pl");
    legend->Draw("same");

    std::string VariableName = itMap->second.first->GetName();

    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".png").c_str(),"png");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".pdf").c_str(),"pdf");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".root").c_str(),"root");
    
  }   

  for(std::map<std::string,histoPair>::const_iterator itMap = histogramResponse.begin(); itMap != histogramResponse.end(); itMap++){

    cCanvas->cd();
    legend->Clear();
    itMap->second.first->Scale(1/itMap->second.first->Integral());   
    itMap->second.second->Scale(1/itMap->second.second->Integral());   

    itMap->second.first->SetLineColor(kBlack);
    itMap->second.first->SetMarkerColor(kBlack);
    itMap->second.first->SetLineWidth(2);
    itMap->second.first->GetXaxis()->SetTitle(itMap->second.first->GetName());
    itMap->second.first->GetYaxis()->SetTitle("a.u.");
    itMap->second.first->SetMarkerStyle(20);
    itMap->second.first->Draw("p");
    itMap->second.second->SetLineColor(kRed);
    itMap->second.second->SetMarkerColor(kRed);
    itMap->second.second->SetLineWidth(2);
    itMap->second.second->SetMarkerStyle(22);
    itMap->second.second->Draw("psame");

    tex->Draw("same");
    tex2->Draw("same");
    tex3->Draw("same");

    legend->AddEntry(itMap->second.first,"New Delphes","pl");
    legend->AddEntry(itMap->second.second,"Old Delphes","pl");
    legend->Draw("same");
    
    std::string VariableName = itMap->second.first->GetName();
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".png").c_str(),"png");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".pdf").c_str(),"pdf");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".root").c_str(),"root");
  */  
  
  return 0 ;

}
