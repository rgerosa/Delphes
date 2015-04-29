#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TClonesArray.h"

#include "classes/DelphesClasses.h"

#define etaBinCut_1   2.5
#define etaBinCut_2   3.0
#define matchingCone  0.3
#define ptBinCut      120
#define leptonIsoThreshold 0.5

using namespace std;

/////////////////////////////////////////////////////                                                                                                                        

int main (int argc, char** argv){

  if(argc < 3 ) {
    cerr<<" to be used as: ./<exe file> <directory with delphes trees> <outputPlot> <letponPT> <jetPT> <PU>"<<endl;
    return -1;
  }

  // Setting for style                                                                                                                                                         
  string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C(1)").c_str());
  }

  gStyle->SetOptStat(111110);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);

  // output folders                                                                                                                                                             
  string outputFileDirectory = argv[2];

  system(("mkdir -p "+outputFileDirectory).c_str());
  system(("rm -r "   +outputFileDirectory+"/*").c_str());

  // input files and chain to be analyzed 
  string inputFileDirectory = argv[1]; 

  TChain *inputChain = new TChain("Delphes");
  inputChain->Add((inputFileDirectory+"/*.root").c_str());

  cout<<"number of events to analyze : "<<inputChain->GetEntries()<<endl;
  
  // threshold for jet and leptons
  float leptonPtThreshold = 0;
  if(argc > 3) 
    leptonPtThreshold = atof(argv[3]);

  float jetPtThreshold = 0;
  if(argc > 4) 
    jetPtThreshold = atof(argv[4]);

  float PU = 0;
  if(argc > 5)
    PU = atof(argv[5]);


  // histogram to fill : basic jet
  TH1F* leadingJetPt  = new TH1F ("leadingJetPt","",35,0,500);
  TH1F* leadingJetEta = new TH1F ("leadingJetEta","",35,-5,5);

  leadingJetPt->Sumw2();
  leadingJetEta->Sumw2();

  TH1F* trailingJetPt  = new TH1F ("trailingJetPt","",35,0,500);
  TH1F* trailingJetEta = new TH1F ("trailingJetEta","",35,-5,5);

  trailingJetPt->Sumw2();
  trailingJetEta->Sumw2();

  TH1F* leadingPuppiJetPt  = new TH1F ("leadingPuppiJetPt","",35,0,500);
  TH1F* leadingPuppiJetEta = new TH1F ("leadingPuppiJetEta","",35,-5,5);

  leadingPuppiJetPt->Sumw2();
  leadingPuppiJetEta->Sumw2();

  TH1F* trailingPuppiJetPt  = new TH1F ("trailingPuppiJetPt","",35,0,500);
  TH1F* trailingPuppiJetEta = new TH1F ("trailingPuppiJetEta","",35,-5,5);

  trailingPuppiJetPt->Sumw2();
  trailingPuppiJetEta->Sumw2();

  // histogram to fill : numerator efficiency
  TH1F* NumeratorReco_vsEta  = new TH1F ("NumeratorReco_vsEta","",35,-5,5);
  TH1F* NumeratorPuppi_vsEta  = new TH1F ("NumeratorPuppi_vsEta","",35,-5,5);

  NumeratorReco_vsEta->Sumw2();
  NumeratorPuppi_vsEta->Sumw2();

  TH1F* NumeratorReco_vsEta_purity   = new TH1F ("NumeratorReco_vsEta_purity","",35,-5,5);
  TH1F* NumeratorPuppi_vsEta_purity  = new TH1F ("NumeratorPuppi_vsEta_purity","",35,-5,5);

  NumeratorReco_vsEta_purity->Sumw2();
  NumeratorPuppi_vsEta_purity->Sumw2();

  TH1F* NumeratorReco_vsPt_central   = new TH1F ("NumeratorReco_vsPt_central","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_central  = new TH1F ("NumeratorPuppi_vsPt_central","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_central->Sumw2();
  NumeratorPuppi_vsPt_central->Sumw2();

  TH1F* NumeratorReco_vsPt_central_purity   = new TH1F ("NumeratorReco_vsPt_central_purity","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_central_purity  = new TH1F ("NumeratorPuppi_vsPt_central_purity","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_central_purity->Sumw2();
  NumeratorPuppi_vsPt_central_purity->Sumw2();

  TH1F* NumeratorReco_vsPt_forward_1   = new TH1F ("NumeratorReco_vsPt_forward_1","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_forward_1  = new TH1F ("NumeratorPuppi_vsPt_forward_1","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_forward_1->Sumw2();
  NumeratorPuppi_vsPt_forward_1->Sumw2();

  TH1F* NumeratorReco_vsPt_forward_2   = new TH1F ("NumeratorReco_vsPt_forward_2","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_forward_2  = new TH1F ("NumeratorPuppi_vsPt_forward_2","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_forward_2->Sumw2();
  NumeratorPuppi_vsPt_forward_2->Sumw2();

  TH1F* NumeratorReco_vsPt_forward_purity_1   = new TH1F ("NumeratorReco_vsPt_forward_purity_1","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_forward_purity_1  = new TH1F ("NumeratorPuppi_vsPt_forward_purity_1","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_forward_purity_1->Sumw2();
  NumeratorPuppi_vsPt_forward_purity_1->Sumw2();

  TH1F* NumeratorReco_vsPt_forward_purity_2   = new TH1F ("NumeratorReco_vsPt_forward_purity_2","",30,jetPtThreshold,500);
  TH1F* NumeratorPuppi_vsPt_forward_purity_2  = new TH1F ("NumeratorPuppi_vsPt_forward_purity_2","",30,jetPtThreshold,500);

  NumeratorReco_vsPt_forward_purity_2->Sumw2();
  NumeratorPuppi_vsPt_forward_purity_2->Sumw2();

  TH1F* NumeratorReco_vsPU_central   = new TH1F ("NumeratorReco_vsPU_central","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_central  = new TH1F ("NumeratorPuppi_vsPU_central","",40,PU-40,PU+40);

  NumeratorReco_vsPU_central->Sumw2();
  NumeratorPuppi_vsPU_central->Sumw2();

  TH1F* NumeratorReco_vsPU_central_purity   = new TH1F ("NumeratorReco_vsPU_central_purity","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_central_purity  = new TH1F ("NumeratorPuppi_vsPU_central_purity","",40,PU-40,PU+40);

  NumeratorReco_vsPU_central_purity->Sumw2();
  NumeratorPuppi_vsPU_central_purity->Sumw2();

  TH1F* NumeratorReco_vsPU_forward_1   = new TH1F ("NumeratorReco_vsPU_forward_1","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_forward_1  = new TH1F ("NumeratorPuppi_vsPU_forward_1","",40,PU-40,PU+40);

  NumeratorReco_vsPU_forward_1->Sumw2();
  NumeratorPuppi_vsPU_forward_1->Sumw2();

  TH1F* NumeratorReco_vsPU_forward_2   = new TH1F ("NumeratorReco_vsPU_forward_2","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_forward_2  = new TH1F ("NumeratorPuppi_vsPU_forward_2","",40,PU-40,PU+40);

  NumeratorReco_vsPU_forward_2->Sumw2();
  NumeratorPuppi_vsPU_forward_2->Sumw2();

  TH1F* NumeratorReco_vsPU_forward_purity_1   = new TH1F ("NumeratorReco_vsPU_forward_purity_1","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_forward_purity_1  = new TH1F ("NumeratorPuppi_vsPU_forward_purity_1","",40,PU-40,PU+40);

  NumeratorReco_vsPU_forward_purity_1->Sumw2();
  NumeratorPuppi_vsPU_forward_purity_1->Sumw2();

  TH1F* NumeratorReco_vsPU_forward_purity_2   = new TH1F ("NumeratorReco_vsPU_forward_purity_2","",40,PU-40,PU+40);
  TH1F* NumeratorPuppi_vsPU_forward_purity_2  = new TH1F ("NumeratorPuppi_vsPU_forward_purity_2","",40,PU-40,PU+40);

  NumeratorReco_vsPU_forward_purity_2->Sumw2();
  NumeratorPuppi_vsPU_forward_purity_2->Sumw2();

  // efficiency denominator

  TH1F* DenominatorReco_vsEta  = new TH1F ("DenominatorReco_vsEta","",35,-5,5);
  TH1F* DenominatorPuppi_vsEta  = new TH1F ("DenominatorPuppi_vsEta","",35,-5,5);

  DenominatorReco_vsEta->Sumw2();
  DenominatorPuppi_vsEta->Sumw2();

  TH1F* DenominatorReco_vsEta_purity  = new TH1F ("DenominatorReco_vsEta_purity","",35,-5,5);
  TH1F* DenominatorPuppi_vsEta_purity  = new TH1F ("DenominatorPuppi_vsEta_purity","",35,-5,5);

  DenominatorReco_vsEta_purity->Sumw2();
  DenominatorPuppi_vsEta_purity->Sumw2();

  TH1F* DenominatorReco_vsPt_central   = new TH1F ("DenominatorReco_vsPt_central","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_central  = new TH1F ("DenominatorPuppi_vsPt_central","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_central->Sumw2();
  DenominatorPuppi_vsPt_central->Sumw2();


  TH1F* DenominatorReco_vsPt_central_purity   = new TH1F ("DenominatorReco_vsPt_central_purity","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_central_purity  = new TH1F ("DenominatorPuppi_vsPt_central_purity","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_central_purity->Sumw2();
  DenominatorPuppi_vsPt_central_purity->Sumw2();

  TH1F* DenominatorReco_vsPt_forward_1   = new TH1F ("DenominatorReco_vsPt_forward_1","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_forward_1  = new TH1F ("DenominatorPuppi_vsPt_forward_1","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_forward_1->Sumw2();
  DenominatorPuppi_vsPt_forward_1->Sumw2();

  TH1F* DenominatorReco_vsPt_forward_2   = new TH1F ("DenominatorReco_vsPt_forward_2","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_forward_2  = new TH1F ("DenominatorPuppi_vsPt_forward_2","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_forward_2->Sumw2();
  DenominatorPuppi_vsPt_forward_2->Sumw2();

  TH1F* DenominatorReco_vsPt_forward_purity_1   = new TH1F ("DenominatorReco_vsPt_forward_purity_1","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_forward_purity_1  = new TH1F ("DenominatorPuppi_vsPt_forward_purity_1","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_forward_purity_1->Sumw2();
  DenominatorPuppi_vsPt_forward_purity_1->Sumw2();

  TH1F* DenominatorReco_vsPt_forward_purity_2   = new TH1F ("DenominatorReco_vsPt_forward_purity_2","",30,jetPtThreshold,500);
  TH1F* DenominatorPuppi_vsPt_forward_purity_2  = new TH1F ("DenominatorPuppi_vsPt_forward_purity_2","",30,jetPtThreshold,500);

  DenominatorReco_vsPt_forward_purity_2->Sumw2();
  DenominatorPuppi_vsPt_forward_purity_2->Sumw2();


  TH1F* DenominatorReco_vsPU_central   = new TH1F ("DenominatorReco_vsPU_central","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_central  = new TH1F ("DenominatorPuppi_vsPU_central","",40,PU-40,PU+40);

  DenominatorReco_vsPU_central->Sumw2();
  DenominatorPuppi_vsPU_central->Sumw2();

  TH1F* DenominatorReco_vsPU_central_purity   = new TH1F ("DenominatorReco_vsPU_central_purity","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_central_purity  = new TH1F ("DenominatorPuppi_vsPU_central_purity","",40,PU-40,PU+40);

  DenominatorReco_vsPU_central_purity->Sumw2();
  DenominatorPuppi_vsPU_central_purity->Sumw2();

  TH1F* DenominatorReco_vsPU_forward_1   = new TH1F ("DenominatorReco_vsPU_forward_1","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_forward_1  = new TH1F ("DenominatorPuppi_vsPU_forward_1","",40,PU-40, PU+40);

  DenominatorReco_vsPU_forward_1->Sumw2();
  DenominatorPuppi_vsPU_forward_1->Sumw2();

  TH1F* DenominatorReco_vsPU_forward_2   = new TH1F ("DenominatorReco_vsPU_forward_2","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_forward_2  = new TH1F ("DenominatorPuppi_vsPU_forward_2","",40,PU-40, PU+40);

  DenominatorReco_vsPU_forward_2->Sumw2();
  DenominatorPuppi_vsPU_forward_2->Sumw2();

  TH1F* DenominatorReco_vsPU_forward_purity_1   = new TH1F ("DenominatorReco_vsPU_forward_purity_1","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_forward_purity_1  = new TH1F ("DenominatorPuppi_vsPU_forward_purity_1","",40,PU-40, PU+40);

  DenominatorReco_vsPU_forward_purity_1->Sumw2();
  DenominatorPuppi_vsPU_forward_purity_1->Sumw2();

  TH1F* DenominatorReco_vsPU_forward_purity_2   = new TH1F ("DenominatorReco_vsPU_forward_purity_2","",40,PU-40,PU+40);
  TH1F* DenominatorPuppi_vsPU_forward_purity_2  = new TH1F ("DenominatorPuppi_vsPU_forward_purity_2","",40,PU-40, PU+40);

  DenominatorReco_vsPU_forward_purity_2->Sumw2();
  DenominatorPuppi_vsPU_forward_purity_2->Sumw2();

  // response

  TH1F* ResponseReco_lowPt_central     = new TH1F ("ResponseReco_lowPt_central","",50,-0.6,0.6);
  TH1F* ResponseReco_lowPt_forward_1   = new TH1F ("ResponseReco_lowPt_forward_1","",50,-0.6,0.6);
  TH1F* ResponseReco_lowPt_forward_2   = new TH1F ("ResponseReco_lowPt_forward_2","",50,-0.6,0.6);
  TH1F* ResponsePuppi_lowPt_central    = new TH1F ("ResponsePuppi_lowPt_central","",50,-0.6,0.6);
  TH1F* ResponsePuppi_lowPt_forward_1  = new TH1F ("ResponsePuppi_lowPt_forward_1","",50,-0.6,0.6);
  TH1F* ResponsePuppi_lowPt_forward_2  = new TH1F ("ResponsePuppi_lowPt_forward_2","",50,-0.6,0.6);

  ResponseReco_lowPt_central->Sumw2();
  ResponseReco_lowPt_forward_1->Sumw2();
  ResponseReco_lowPt_forward_2->Sumw2();
  ResponsePuppi_lowPt_central->Sumw2();
  ResponsePuppi_lowPt_forward_1->Sumw2();
  ResponsePuppi_lowPt_forward_2->Sumw2();

  TH1F* ResponseReco_highPt_central     = new TH1F ("ResponseReco_highPt_central","",50,-0.6,0.6);
  TH1F* ResponseReco_highPt_forward_1   = new TH1F ("ResponseReco_highPt_forward_1","",50,-0.6,0.6);
  TH1F* ResponseReco_highPt_forward_2   = new TH1F ("ResponseReco_highPt_forward_2","",50,-0.6,0.6);
  TH1F* ResponsePuppi_highPt_central    = new TH1F ("ResponsePuppi_highPt_central","",50,-0.6,0.6);
  TH1F* ResponsePuppi_highPt_forward_1  = new TH1F ("ResponsePuppi_highPt_forward_1","",50,-0.6,0.6);
  TH1F* ResponsePuppi_highPt_forward_2  = new TH1F ("ResponsePuppi_highPt_forward_2","",50,-0.6,0.6);

  ResponseReco_highPt_central->Sumw2();
  ResponseReco_highPt_forward_1->Sumw2();
  ResponseReco_highPt_forward_2->Sumw2();
  ResponsePuppi_highPt_central->Sumw2();
  ResponsePuppi_highPt_forward_1->Sumw2();
  ResponsePuppi_highPt_forward_2->Sumw2();

  TH1F* ResponseReco_deta   = new TH1F ("ResponseReco_deta","",50,-0.2,0.2);
  TH1F* ResponsePuppi_deta  = new TH1F ("ResponsePuppi_deta","",50,-0.2,0.2);

  ResponseReco_deta->Sumw2();
  ResponsePuppi_deta->Sumw2();

  TH1F* ResponseReco_mjj   = new TH1F ("ResponseReco_mjj","",50,-0.6,0.6);
  TH1F* ResponsePuppi_mjj  = new TH1F ("ResponsePuppi_mjj","",50,-0.6,0.6);

  ResponseReco_mjj->Sumw2();
  ResponsePuppi_mjj->Sumw2();

  // set the branch address to be used
  TClonesArray* Muons = new TClonesArray("Muon");
  inputChain->SetBranchAddress("Muon",&Muons);

  TClonesArray* Electrons = new TClonesArray("Electron");
  inputChain->SetBranchAddress("Electron",&Electrons);

  TClonesArray* RecoJets = new TClonesArray("Jet");
  inputChain->SetBranchAddress("JetPUID",&RecoJets);

  TClonesArray* PuppiJets = new TClonesArray("Jet");
  inputChain->SetBranchAddress("PuppiJetPUID",&PuppiJets);

  TClonesArray* GenJets = new TClonesArray("Jet");
  inputChain->SetBranchAddress("GenJet",&GenJets);

  TClonesArray* NPU = new TClonesArray("ScalarHT");
  inputChain->SetBranchAddress("NPU",&NPU);

  inputChain->SetBranchStatus("LHE*",0);
  inputChain->SetBranchStatus("Track*",0);
  inputChain->SetBranchStatus("Rho*",0);
  inputChain->SetBranchStatus("*MissingET*",0);
  
 
  // loop on the events and fill the histos
  vector<Muon*> tightMuons;
  vector<Electron*> tightElectrons;

  vector<Jet*> cleanedJets;
  vector<Jet*> cleanedPuppiJets;
  vector<Jet*> cleanedGenJets;

  for(int iEntry = 0; iEntry < inputChain->GetEntries() ; iEntry++){

    tightMuons.clear();
    tightElectrons.clear();
    cleanedJets.clear();
    cleanedPuppiJets.clear();
    cleanedGenJets.clear();

    if(iEntry%10000 == 0) cout<<"reading entry "<<iEntry<<endl;

    inputChain->GetEntry(iEntry);

    int nLep = 0;
    int nJet = 0;
    
    // require at least two leptons over a fixed PT cut, no isolation required    
    for(int iMuon = 0; iMuon < Muons->GetEntries(); iMuon++){
      Muon* muon = (Muon*) Muons->At(iMuon);      
      if(muon->PT >= leptonPtThreshold and muon->IsolationVarRhoCorr <= leptonIsoThreshold){	
	nLep++; 
        tightMuons.push_back(muon);
      }
    }

    for(int iElectron = 0; iElectron < Electrons->GetEntries(); iElectron++){
      Electron* electron = (Electron*) Electrons->At(iElectron);
      if(electron->PT >= leptonPtThreshold and electron->IsolationVarRhoCorr <= leptonIsoThreshold){
	nLep++; 
        tightElectrons.push_back(electron);
      }
    }

    if(nLep < 2 or nLep >=3) continue ; // require two tight leptons (no loose lepton veto here)
    
    // require at least two jets (reco) over a fixed PT cut cleaed
    for(int iJet = 0; iJet < RecoJets->GetEntries(); iJet++){

      Jet* jet = (Jet*) RecoJets->At(iJet);

      if(jet->PT >= jetPtThreshold){
        bool discardJet = false;
        
        for(size_t iMuon = 0; iMuon < tightMuons.size(); iMuon++){
	  if(tightMuons.at(iMuon)->P4().DeltaR(jet->P4()) < matchingCone){
	    discardJet = true;
	  }
	}

        if(discardJet) continue ;

        for(size_t iElectron = 0; iElectron < tightElectrons.size(); iElectron++){
	  if(tightElectrons.at(iElectron)->P4().DeltaR(jet->P4()) < matchingCone){
	    discardJet = true;
	  }
	}
        if(discardJet) continue ;
	nJet++; 
	cleanedJets.push_back(jet);
      }
    }

    if(nJet < 2) continue ;

    // require at least two jets (reco) over a fixed PT cut
    nJet = 0;

    for(int iJet = 0; iJet < PuppiJets->GetEntries(); iJet++){

      Jet* jet = (Jet*) PuppiJets->At(iJet);

      if(jet->PT >= jetPtThreshold){

        bool discardJet = false;

        for(size_t iMuon = 0; iMuon < tightMuons.size(); iMuon++){
	  if(tightMuons.at(iMuon)->P4().DeltaR(jet->P4()) < matchingCone){
	    discardJet = true;
	  }
	}

        if(discardJet) continue ;

        for(size_t iElectron = 0; iElectron < tightElectrons.size(); iElectron++){
	  if(tightElectrons.at(iElectron)->P4().DeltaR(jet->P4()) < matchingCone){
	    discardJet = true;
	  }
	}
        if(discardJet) continue ;
	nJet++; 
	cleanedPuppiJets.push_back(jet);
      }
    }
      

    if(nJet < 2) continue ;

    
    // fill gen jets removing leptons
    for(int iJet = 0; iJet < GenJets->GetEntries(); iJet++){

      Jet* jet = (Jet*) GenJets->At(iJet);

      bool discardJet = false;

      for(size_t iMuon = 0; iMuon < tightMuons.size(); iMuon++){
	if(tightMuons.at(iMuon)->P4().DeltaR(jet->P4()) < matchingCone){
	  discardJet = true;
	}
      }

      if(discardJet) continue ;
      
      for(size_t iElectron = 0; iElectron < tightElectrons.size(); iElectron++){
	if(tightElectrons.at(iElectron)->P4().DeltaR(jet->P4()) < matchingCone){
	  discardJet = true;
	}
      }
      if(discardJet) continue ;

      cleanedGenJets.push_back(jet);
    } 

    if(cleanedGenJets.size() < 2) continue ;

    // fill histo
    ScalarHT* npu = (ScalarHT*) NPU->At(0);
      
    // leading and trailing reco jets
    leadingJetPt->Fill(cleanedJets.at(0)->PT);
    leadingJetEta->Fill(cleanedJets.at(0)->Eta);

    trailingJetPt->Fill(cleanedJets.at(1)->PT);
    trailingJetEta->Fill(cleanedJets.at(1)->Eta);

    // leading and trailing puppi jets
    leadingPuppiJetPt->Fill(cleanedPuppiJets.at(0)->PT);
    leadingPuppiJetEta->Fill(cleanedPuppiJets.at(0)->Eta);

    trailingPuppiJetPt->Fill(cleanedPuppiJets.at(1)->PT);
    trailingPuppiJetEta->Fill(cleanedPuppiJets.at(1)->Eta);

    // numerartors
    DenominatorReco_vsEta->Fill(cleanedJets.at(0)->Eta);
    DenominatorReco_vsEta->Fill(cleanedJets.at(1)->Eta);
    DenominatorPuppi_vsEta->Fill(cleanedPuppiJets.at(0)->Eta);
    DenominatorPuppi_vsEta->Fill(cleanedPuppiJets.at(1)->Eta);

    DenominatorReco_vsEta_purity->Fill(cleanedGenJets.at(0)->Eta);
    DenominatorReco_vsEta_purity->Fill(cleanedGenJets.at(1)->Eta);
    DenominatorPuppi_vsEta_purity->Fill(cleanedGenJets.at(0)->Eta);
    DenominatorPuppi_vsEta_purity->Fill(cleanedGenJets.at(1)->Eta);

    // leading    
    if(fabs(cleanedJets.at(0)->Eta) < etaBinCut_1){       
      DenominatorReco_vsPt_central->Fill(cleanedJets.at(0)->PT);
      DenominatorReco_vsPU_central->Fill(npu->HT);
    }
    else if(fabs(cleanedJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(0)->Eta) < etaBinCut_2){
      DenominatorReco_vsPt_forward_1->Fill(cleanedJets.at(0)->PT);
      DenominatorReco_vsPU_forward_1->Fill(npu->HT);
    }
    else if(fabs(cleanedJets.at(0)->Eta) >= etaBinCut_2){
      DenominatorReco_vsPt_forward_2->Fill(cleanedJets.at(0)->PT);
      DenominatorReco_vsPU_forward_2->Fill(npu->HT);
    }

    ////////  
    if(fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_1){       
      DenominatorPuppi_vsPt_central->Fill(cleanedPuppiJets.at(0)->PT);
      DenominatorPuppi_vsPU_central->Fill(npu->HT);
    }
    else if(fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_2){
      DenominatorPuppi_vsPt_forward_1->Fill(cleanedPuppiJets.at(0)->PT);
      DenominatorPuppi_vsPU_forward_1->Fill(npu->HT);
    }
    else if(fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_2){
      DenominatorPuppi_vsPt_forward_2->Fill(cleanedPuppiJets.at(0)->PT);
      DenominatorPuppi_vsPU_forward_2->Fill(npu->HT);
    }

    /////
    if(fabs(cleanedJets.at(1)->Eta) < etaBinCut_1) {      
      DenominatorReco_vsPt_central->Fill(cleanedJets.at(1)->PT);
      DenominatorReco_vsPU_central->Fill(npu->HT);
    }
    else if(fabs(cleanedJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(1)->Eta) < etaBinCut_2){
      DenominatorReco_vsPt_forward_1->Fill(cleanedJets.at(1)->PT);
      DenominatorReco_vsPU_forward_1->Fill(npu->HT);
    }
    else if(fabs(cleanedJets.at(1)->Eta) >= etaBinCut_2){
      DenominatorReco_vsPt_forward_2->Fill(cleanedJets.at(1)->PT);
      DenominatorReco_vsPU_forward_2->Fill(npu->HT);
    }

    /////
    if(fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_1 ) {
      DenominatorPuppi_vsPt_central->Fill(cleanedPuppiJets.at(1)->PT);
      DenominatorPuppi_vsPU_central->Fill(npu->HT);
    }
    else if(fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_2){
      DenominatorPuppi_vsPt_forward_1->Fill(cleanedPuppiJets.at(1)->PT);
      DenominatorPuppi_vsPU_forward_1->Fill(npu->HT);
    }
    else if(fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_2){
      DenominatorPuppi_vsPt_forward_2->Fill(cleanedPuppiJets.at(1)->PT);
      DenominatorPuppi_vsPU_forward_2->Fill(npu->HT);
    }

    // leading    
    if(fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_1 ){       
      DenominatorReco_vsPt_central_purity->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPt_central_purity->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPU_central_purity->Fill(npu->HT);
      DenominatorReco_vsPU_central_purity->Fill(npu->HT);
    }
    else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_2){
      DenominatorReco_vsPt_forward_purity_1->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPt_forward_purity_1->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPU_forward_purity_1->Fill(npu->HT);
      DenominatorReco_vsPU_forward_purity_1->Fill(npu->HT);
    }
    else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_2){
      DenominatorReco_vsPt_forward_purity_2->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPt_forward_purity_2->Fill(cleanedGenJets.at(0)->PT);
      DenominatorPuppi_vsPU_forward_purity_2->Fill(npu->HT);
      DenominatorReco_vsPU_forward_purity_2->Fill(npu->HT);
    }
    
    // trailing
    if(fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_1 ){       
      DenominatorReco_vsPt_central_purity->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPt_central_purity->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPU_central_purity->Fill(npu->HT);
      DenominatorReco_vsPU_central_purity->Fill(npu->HT);
    }
    else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_2){
      DenominatorReco_vsPt_forward_purity_1->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPt_forward_purity_1->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPU_forward_purity_1->Fill(npu->HT);
      DenominatorReco_vsPU_forward_purity_1->Fill(npu->HT);
    }
    else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_2){
      DenominatorReco_vsPt_forward_purity_2->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPt_forward_purity_2->Fill(cleanedGenJets.at(1)->PT);
      DenominatorPuppi_vsPU_forward_purity_2->Fill(npu->HT);
      DenominatorReco_vsPU_forward_purity_2->Fill(npu->HT);
    }

    // event loop
    bool recoMatched_1  = false;
    bool recoMatched_2  = false;

    bool puppiMatched_1 = false;
    bool puppiMatched_2 = false;

    Jet* matchedGenJet_1 = 0;
    Jet* matchedGenJet_2 = 0;

    Jet* matchedGenJet_puppi_1 = 0;
    Jet* matchedGenJet_puppi_2 = 0;

    for(size_t iJet = 0; iJet < cleanedGenJets.size(); iJet++){
      Jet* jet = cleanedGenJets.at(iJet);

      if(jet->P4().DeltaR(cleanedJets.at(0)->P4()) < matchingCone and !recoMatched_1){ // positive match with leading jet

        recoMatched_1 = true;
        matchedGenJet_1 = jet;

	NumeratorReco_vsEta->Fill(cleanedJets.at(0)->Eta);
 
        if(cleanedJets.at(0)->PT < ptBinCut and fabs(cleanedJets.at(0)->Eta) < etaBinCut_1) {
	  ResponseReco_lowPt_central->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);	
	}
        else if(cleanedJets.at(0)->PT < ptBinCut and fabs(cleanedJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(0)->Eta) < etaBinCut_2){
	  ResponseReco_lowPt_forward_1->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);       
	}
	else if(cleanedJets.at(0)->PT < ptBinCut and fabs(cleanedJets.at(0)->Eta) >= etaBinCut_2){
	  ResponseReco_lowPt_forward_2->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);       
	}
        else if(cleanedJets.at(0)->PT >= ptBinCut and fabs(cleanedJets.at(0)->Eta) < etaBinCut_1){
	  ResponseReco_highPt_central->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedJets.at(0)->PT >= ptBinCut and fabs(cleanedJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(0)->Eta) < etaBinCut_2){
	  ResponseReco_highPt_forward_1->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedJets.at(0)->PT >= ptBinCut and fabs(cleanedJets.at(0)->Eta) >= etaBinCut_2){
	  ResponseReco_highPt_forward_2->Fill((cleanedJets.at(0)->PT-jet->PT)/jet->PT);
	}


	if(fabs(cleanedJets.at(0)->Eta) < etaBinCut_1 ){
	  NumeratorReco_vsPt_central->Fill(cleanedJets.at(0)->PT);
	  NumeratorReco_vsPU_central->Fill(npu->HT);
	}
	else if(fabs(cleanedJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(0)->Eta) < etaBinCut_2){
	  NumeratorReco_vsPt_forward_1->Fill(cleanedJets.at(0)->PT);
	  NumeratorReco_vsPU_forward_1->Fill(npu->HT);
	}
	else if(fabs(cleanedJets.at(0)->Eta) >= etaBinCut_2){
	  NumeratorReco_vsPt_forward_2->Fill(cleanedJets.at(0)->PT);
	  NumeratorReco_vsPU_forward_2->Fill(npu->HT);
	}
      } 
	  
      else if(jet->P4().DeltaR(cleanedJets.at(1)->P4()) < matchingCone  and !recoMatched_2){

        recoMatched_2 = true;
        matchedGenJet_2 = jet ;
	NumeratorReco_vsEta->Fill(cleanedJets.at(1)->Eta);

        if(cleanedJets.at(1)->PT < ptBinCut and fabs(cleanedJets.at(1)->Eta) < etaBinCut_1) {
	  ResponseReco_lowPt_central->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);	
	}
        else if(cleanedJets.at(1)->PT < ptBinCut and fabs(cleanedJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(1)->Eta) < etaBinCut_2){
	  ResponseReco_lowPt_forward_1->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);       
	}
	else if(cleanedJets.at(1)->PT < ptBinCut and fabs(cleanedJets.at(1)->Eta) >= etaBinCut_2){
	  ResponseReco_lowPt_forward_2->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);       
	}
        else if(cleanedJets.at(1)->PT >= ptBinCut and fabs(cleanedJets.at(1)->Eta) < etaBinCut_1){
	  ResponseReco_highPt_central->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedJets.at(1)->PT >= ptBinCut and fabs(cleanedJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(1)->Eta) < etaBinCut_2){
	  ResponseReco_highPt_forward_1->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedJets.at(1)->PT >= ptBinCut and fabs(cleanedJets.at(1)->Eta) >= etaBinCut_2){
	  ResponseReco_highPt_forward_2->Fill((cleanedJets.at(1)->PT-jet->PT)/jet->PT);
	}


	if(fabs(cleanedJets.at(1)->Eta) < etaBinCut_1 ){
	  NumeratorReco_vsPt_central->Fill(cleanedJets.at(1)->PT);
	  NumeratorReco_vsPU_central->Fill(npu->HT);
	}
	else if(fabs(cleanedJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedJets.at(1)->Eta) < etaBinCut_2){
	  NumeratorReco_vsPt_forward_1->Fill(cleanedJets.at(1)->PT);
	  NumeratorReco_vsPU_forward_1->Fill(npu->HT);
	}
	else if(fabs(cleanedJets.at(1)->Eta) >= etaBinCut_2){
	  NumeratorReco_vsPt_forward_2->Fill(cleanedJets.at(1)->PT);
	  NumeratorReco_vsPU_forward_2->Fill(npu->HT);
	}
      }	  

      if(jet->P4().DeltaR(cleanedPuppiJets.at(0)->P4()) < matchingCone and !puppiMatched_1 ){ // positive match with leading jet

        puppiMatched_1 = true;
        matchedGenJet_puppi_1 = jet;

	NumeratorPuppi_vsEta->Fill(cleanedPuppiJets.at(0)->Eta);

        if(cleanedPuppiJets.at(0)->PT < ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_1) {
	  ResponsePuppi_lowPt_central->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);	
	}
        else if(cleanedPuppiJets.at(0)->PT < ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_2){
	  ResponsePuppi_lowPt_forward_1->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);       
	}
	else if(cleanedPuppiJets.at(0)->PT < ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_2){
	  ResponsePuppi_lowPt_forward_2->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);       
	}
        else if(cleanedPuppiJets.at(0)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_1){
	  ResponsePuppi_highPt_central->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedPuppiJets.at(0)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_2){
	  ResponsePuppi_highPt_forward_1->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedPuppiJets.at(0)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_2){
	  ResponsePuppi_highPt_forward_2->Fill((cleanedPuppiJets.at(0)->PT-jet->PT)/jet->PT);
	}


	if(fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_1 ){
	  NumeratorPuppi_vsPt_central->Fill(cleanedPuppiJets.at(0)->PT);
	  NumeratorPuppi_vsPU_central->Fill(npu->HT);
	}
	else if(fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(0)->Eta) < etaBinCut_2){
	  NumeratorPuppi_vsPt_forward_1->Fill(cleanedPuppiJets.at(0)->PT);
	  NumeratorPuppi_vsPU_forward_1->Fill(npu->HT);
	}
	else if(fabs(cleanedPuppiJets.at(0)->Eta) >= etaBinCut_2){
	  NumeratorPuppi_vsPt_forward_2->Fill(cleanedPuppiJets.at(0)->PT);
	  NumeratorPuppi_vsPU_forward_2->Fill(npu->HT);
	}
      }

      else if(jet->P4().DeltaR(cleanedPuppiJets.at(1)->P4()) < matchingCone and !puppiMatched_2){

	puppiMatched_2 = true ;
        matchedGenJet_puppi_2 = jet;

	NumeratorPuppi_vsEta->Fill(cleanedPuppiJets.at(1)->Eta);

        if(cleanedPuppiJets.at(1)->PT < ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_1) {
	  ResponsePuppi_lowPt_central->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);	
	}
        else if(cleanedPuppiJets.at(1)->PT < ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_2){
	  ResponsePuppi_lowPt_forward_1->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);       
	}
	else if(cleanedPuppiJets.at(1)->PT < ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_2){
	  ResponsePuppi_lowPt_forward_2->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);       
	}
        else if(cleanedPuppiJets.at(1)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_1){
	  ResponsePuppi_highPt_central->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedPuppiJets.at(1)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_2){
	  ResponsePuppi_highPt_forward_1->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);
	}
        else if(cleanedPuppiJets.at(1)->PT >= ptBinCut and fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_2){
	  ResponsePuppi_highPt_forward_2->Fill((cleanedPuppiJets.at(1)->PT-jet->PT)/jet->PT);
	}


	if(fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_1 ){
	  NumeratorPuppi_vsPt_central->Fill(cleanedPuppiJets.at(1)->PT);
	  NumeratorPuppi_vsPU_central->Fill(npu->HT);
	}
	else if(fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedPuppiJets.at(1)->Eta) < etaBinCut_2){
	  NumeratorPuppi_vsPt_forward_1->Fill(cleanedPuppiJets.at(1)->PT);
	  NumeratorPuppi_vsPU_forward_1->Fill(npu->HT);
	}
	else if(fabs(cleanedPuppiJets.at(1)->Eta) >= etaBinCut_2){
	  NumeratorPuppi_vsPt_forward_2->Fill(cleanedPuppiJets.at(1)->PT);
	  NumeratorPuppi_vsPU_forward_2->Fill(npu->HT);
	}

      }
    }
	  
    // response for deta and mjj
    if(matchedGenJet_1!=0 and matchedGenJet_2!=0){
      ResponseReco_deta->Fill((fabs(cleanedJets.at(0)->Eta-cleanedJets.at(1)->Eta)-fabs(matchedGenJet_1->Eta-matchedGenJet_2->Eta))/fabs(matchedGenJet_1->Eta-matchedGenJet_2->Eta));
      ResponseReco_mjj->Fill(((cleanedJets.at(0)->P4()+cleanedJets.at(1)->P4()).M()-(matchedGenJet_1->P4()+matchedGenJet_2->P4()).M())/(matchedGenJet_1->P4()+matchedGenJet_2->P4()).M());
    }
 
    if(matchedGenJet_puppi_1!=0 and matchedGenJet_puppi_2!=0){

      ResponsePuppi_deta->Fill((fabs(cleanedPuppiJets.at(0)->Eta-cleanedPuppiJets.at(1)->Eta)-fabs(matchedGenJet_puppi_1->Eta-matchedGenJet_puppi_2->Eta))/fabs(matchedGenJet_puppi_1->Eta-matchedGenJet_puppi_2->Eta));

      ResponsePuppi_mjj->Fill(((cleanedPuppiJets.at(0)->P4()+cleanedPuppiJets.at(1)->P4()).M()-(matchedGenJet_puppi_1->P4()+matchedGenJet_puppi_2->P4()).M())/(matchedGenJet_puppi_1->P4()+matchedGenJet_puppi_2->P4()).M());
    }
    
    // for purity
    bool genMatched_1 = false;
    bool genMatched_2 = false;
    
    for(size_t iJet = 0; iJet < cleanedJets.size(); iJet++){

      Jet* jet = cleanedJets.at(iJet);

      if(jet->P4().DeltaR(cleanedGenJets.at(0)->P4()) < matchingCone and !genMatched_1){ // positive match with leading jet

        genMatched_1 = true;

	NumeratorReco_vsEta_purity->Fill(cleanedGenJets.at(0)->Eta);
 
	if(fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_1 ){
  	  NumeratorReco_vsPt_central_purity->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorReco_vsPU_central_purity->Fill(npu->HT);
	}
	else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_2  ){
	  NumeratorReco_vsPt_forward_purity_1->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorReco_vsPU_forward_purity_1->Fill(npu->HT);
	}
        else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_2  ){
	  NumeratorReco_vsPt_forward_purity_2->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorReco_vsPU_forward_purity_2->Fill(npu->HT);
	}
      }
      else if(jet->P4().DeltaR(cleanedGenJets.at(1)->P4()) < matchingCone  and !genMatched_2){

        genMatched_2 = true;

	NumeratorReco_vsEta_purity->Fill(cleanedGenJets.at(1)->Eta);

	if(fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_1 ){
	  NumeratorReco_vsPt_central_purity->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorReco_vsPU_central_purity->Fill(npu->HT);
	}
	else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_2  ){
	  NumeratorReco_vsPt_forward_purity_1->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorReco_vsPU_forward_purity_1->Fill(npu->HT);
	}
        else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_2  ){
	  NumeratorReco_vsPt_forward_purity_2->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorReco_vsPU_forward_purity_2->Fill(npu->HT);
	}
      }
    }

    genMatched_1 = false;
    genMatched_2 = false;


    for(size_t iJet = 0; iJet < cleanedPuppiJets.size(); iJet++){

      Jet* jet = cleanedPuppiJets.at(iJet);

      if(jet->P4().DeltaR(cleanedGenJets.at(0)->P4()) < matchingCone and !genMatched_1 ){ // positive match with leading jet

        genMatched_1 = true;

	NumeratorPuppi_vsEta_purity->Fill(cleanedGenJets.at(0)->Eta);

	if(fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_1 ){
  	  NumeratorPuppi_vsPt_central_purity->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorPuppi_vsPU_central_purity->Fill(npu->HT);
	}
	else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(0)->Eta) < etaBinCut_2  ){
	  NumeratorPuppi_vsPt_forward_purity_1->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorPuppi_vsPU_forward_purity_1->Fill(npu->HT);
	}
        else if(fabs(cleanedGenJets.at(0)->Eta) >= etaBinCut_2  ){
	  NumeratorPuppi_vsPt_forward_purity_2->Fill(cleanedGenJets.at(0)->PT);
	  NumeratorPuppi_vsPU_forward_purity_2->Fill(npu->HT);
	}

      }

      else if(jet->P4().DeltaR(cleanedGenJets.at(1)->P4()) < matchingCone and !genMatched_2){

	genMatched_2 = true ;

	NumeratorPuppi_vsEta_purity->Fill(cleanedGenJets.at(1)->Eta);

	if(fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_1 ){
	  NumeratorPuppi_vsPt_central_purity->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorPuppi_vsPU_central_purity->Fill(npu->HT);
	}
	else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_1 and fabs(cleanedGenJets.at(1)->Eta) < etaBinCut_2  ){
	  NumeratorPuppi_vsPt_forward_purity_1->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorPuppi_vsPU_forward_purity_1->Fill(npu->HT);
	}
        else if(fabs(cleanedGenJets.at(1)->Eta) >= etaBinCut_2  ){
	  NumeratorPuppi_vsPt_forward_purity_2->Fill(cleanedGenJets.at(1)->PT);
	  NumeratorPuppi_vsPU_forward_purity_2->Fill(npu->HT);
	}

      }
    }
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

  TLatex * tex = new TLatex(0.94,0.92," 14 TeV");
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

  // leading jet pt
  leadingJetPt->SetLineWidth(2);
  leadingJetPt->SetLineColor(kBlue);
  leadingJetPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  leadingJetPt->GetYaxis()->SetTitle("Entries");
  leadingJetPt->Draw("hist");

  gPad->Update();
  TPaveStats *tps1 = (TPaveStats*) leadingJetPt->FindObject("stats");
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
  
  leadingPuppiJetPt->SetLineWidth(2);
  leadingPuppiJetPt->SetLineColor(kRed);
  leadingPuppiJetPt->SetLineStyle(7);
  leadingPuppiJetPt->Draw("hist");
  
  gPad->Update();
  TPaveStats *tps2 = (TPaveStats*) leadingPuppiJetPt->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  leadingJetPt->GetYaxis()->SetRangeUser(0.001,max(leadingJetPt->GetMaximum(),leadingPuppiJetPt->GetMaximum())*1.25);

  leadingJetPt->Draw("hist");
  leadingPuppiJetPt->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(leadingJetPt,"leading CHS jet","l");
  legend->AddEntry(leadingPuppiJetPt,"leading Puppi jet","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetPt_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // trailing jet pt
  trailingJetPt->SetLineWidth(2);
  trailingJetPt->SetLineColor(kBlue);
  trailingJetPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  trailingJetPt->GetYaxis()->SetTitle("Entries");
  trailingJetPt->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) trailingJetPt->FindObject("stats");
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
  
  trailingPuppiJetPt->SetLineWidth(2);
  trailingPuppiJetPt->SetLineColor(kRed);
  trailingPuppiJetPt->SetLineStyle(7);
  trailingPuppiJetPt->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) trailingPuppiJetPt->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  trailingJetPt->GetYaxis()->SetRangeUser(0.001,max(trailingJetPt->GetMaximum(),trailingPuppiJetPt->GetMaximum())*1.25);

  trailingJetPt->Draw("hist");
  trailingPuppiJetPt->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(trailingJetPt,"trailing CHS jet","l");
  legend->AddEntry(trailingPuppiJetPt,"trailing Puppi jet","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetPt_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  //leading jet eta
  leadingJetEta->SetLineWidth(2);
  leadingJetEta->SetLineColor(kBlue);
  leadingJetEta->GetXaxis()->SetTitle("#eta");
  leadingJetEta->GetYaxis()->SetTitle("Entries");
  leadingJetEta->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) leadingJetEta->FindObject("stats");
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
  
  leadingPuppiJetEta->SetLineWidth(2);
  leadingPuppiJetEta->SetLineColor(kRed);
  leadingPuppiJetEta->SetLineStyle(7);
  leadingPuppiJetEta->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) leadingPuppiJetEta->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  leadingJetEta->GetYaxis()->SetRangeUser(0.001,max(leadingJetEta->GetMaximum(),leadingPuppiJetEta->GetMaximum())*1.25);

  leadingJetEta->Draw("hist");
  leadingPuppiJetEta->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(leadingJetEta,"leading CHS jet","l");
  legend->AddEntry(leadingPuppiJetEta,"leading Puppi jet","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/leadingJetEta_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //leading jet eta
  trailingJetEta->SetLineWidth(2);
  trailingJetEta->SetLineColor(kBlue);
  trailingJetEta->GetXaxis()->SetTitle("#eta");
  trailingJetEta->GetYaxis()->SetTitle("Entries");
  trailingJetEta->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) trailingJetEta->FindObject("stats");
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
  
  trailingPuppiJetEta->SetLineWidth(2);
  trailingPuppiJetEta->SetLineColor(kRed);
  trailingPuppiJetEta->SetLineStyle(7);
  trailingPuppiJetEta->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) trailingPuppiJetEta->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  trailingJetEta->GetYaxis()->SetRangeUser(0.001,max(trailingJetEta->GetMaximum(),trailingPuppiJetEta->GetMaximum())*1.25);

  trailingJetEta->Draw("hist");
  trailingPuppiJetEta->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(trailingJetEta,"trailing CHS jet","l");
  legend->AddEntry(trailingPuppiJetEta,"trailing Puppi jet","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/trailingJetEta_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  ///////////////////////////////////////
  // efficiency plots ///////////////////
  ///////////////////////////////////////

  TGraphAsymmErrors* EfficiencyReco_vsEta = new TGraphAsymmErrors();
  EfficiencyReco_vsEta->SetName("EfficiencyReco_vsEta");

  TGraphAsymmErrors* EfficiencyPuppi_vsEta = new TGraphAsymmErrors();
  EfficiencyPuppi_vsEta->SetName("EfficiencyPuppi_vsEta");

  EfficiencyReco_vsEta->BayesDivide(NumeratorReco_vsEta,DenominatorReco_vsEta);
  EfficiencyPuppi_vsEta->BayesDivide(NumeratorPuppi_vsEta,DenominatorPuppi_vsEta);


  TH2F* frameEfficiency_vsEta = new TH2F("frameEfficiency_vsEta","",50,-5,5,500,0.6,
					 1.1*max(EfficiencyReco_vsEta->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsEta->GetHistogram()->GetMaximum()));
  frameEfficiency_vsEta->SetLineWidth(2);
  frameEfficiency_vsEta->SetMarkerStyle(21);
  frameEfficiency_vsEta->SetMarkerSize(0.3);
  frameEfficiency_vsEta->GetXaxis()->SetTitle("#eta");
  frameEfficiency_vsEta->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsEta->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsEta->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsEta->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsEta->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsEta->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsEta->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsEta->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsEta->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsEta->SetStats(0);
  frameEfficiency_vsEta->Draw();

  EfficiencyReco_vsEta->SetLineWidth(2);
  EfficiencyReco_vsEta->SetLineColor(kBlue);
  EfficiencyReco_vsEta->SetMarkerStyle(20);
  EfficiencyReco_vsEta->SetMarkerColor(kBlue);
  EfficiencyReco_vsEta->GetXaxis()->SetTitle("#eta");
  EfficiencyReco_vsEta->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsEta->Draw("p");

  EfficiencyPuppi_vsEta->SetLineWidth(2);
  EfficiencyPuppi_vsEta->SetLineColor(kRed);
  EfficiencyPuppi_vsEta->SetMarkerStyle(24);
  EfficiencyPuppi_vsEta->SetMarkerColor(kRed);
  EfficiencyPuppi_vsEta->GetXaxis()->SetTitle("#eta");
  EfficiencyPuppi_vsEta->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsEta->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsEta,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsEta,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsEta_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();
  
  // efficiency vs pt
  TGraphAsymmErrors* EfficiencyReco_vsPt_central = new TGraphAsymmErrors();
  EfficiencyReco_vsPt_central->SetName("EfficiencyReco_vsPt_central");

  TGraphAsymmErrors* EfficiencyPuppi_vsPt_central = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPt_central->SetName("EfficiencyPuppi_vsPt_central");

  EfficiencyReco_vsPt_central->BayesDivide(NumeratorReco_vsPt_central,DenominatorReco_vsPt_central);
  EfficiencyPuppi_vsPt_central->BayesDivide(NumeratorPuppi_vsPt_central,DenominatorPuppi_vsPt_central);

  
  TH2F* frameEfficiency_vsPt = new TH2F("frameEfficiency_vsPt","",30,jetPtThreshold,500,500,
					min(EfficiencyReco_vsPt_central->GetHistogram()->GetMinimum(),EfficiencyPuppi_vsPt_central->GetHistogram()->GetMinimum()),
					1.05*max(EfficiencyReco_vsPt_central->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPt_central->GetHistogram()->GetMaximum()));
  frameEfficiency_vsPt->SetLineWidth(2);
  frameEfficiency_vsPt->SetMarkerStyle(21);
  frameEfficiency_vsPt->SetMarkerSize(0.3);
  frameEfficiency_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  frameEfficiency_vsPt->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPt->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPt->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPt->SetStats(0);
  frameEfficiency_vsPt->Draw();

  EfficiencyReco_vsPt_central->SetLineWidth(2);
  EfficiencyReco_vsPt_central->SetLineColor(kBlue);
  EfficiencyReco_vsPt_central->SetMarkerStyle(20);
  EfficiencyReco_vsPt_central->SetMarkerColor(kBlue);
  EfficiencyReco_vsPt_central->GetXaxis()->SetTitle("p_{T} (GeV)");
  EfficiencyReco_vsPt_central->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPt_central->Draw("p");

  EfficiencyPuppi_vsPt_central->SetLineWidth(2);
  EfficiencyPuppi_vsPt_central->SetLineColor(kRed);
  EfficiencyPuppi_vsPt_central->SetMarkerStyle(24);
  EfficiencyPuppi_vsPt_central->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPt_central->GetXaxis()->SetTitle("p_{T} (GeV)");
  EfficiencyPuppi_vsPt_central->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPt_central->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPt_central,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPt_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* EfficiencyReco_vsPt_forward_1 = new TGraphAsymmErrors();
  EfficiencyReco_vsPt_forward_1->SetName("EfficiencyReco_vsPt_forward_1");

  TGraphAsymmErrors* EfficiencyPuppi_vsPt_forward_1 = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPt_forward_1->SetName("EfficiencyPuppi_vsPt_forward_1");

  EfficiencyReco_vsPt_forward_1->BayesDivide(NumeratorReco_vsPt_forward_1,DenominatorReco_vsPt_forward_1);
  EfficiencyPuppi_vsPt_forward_1->BayesDivide(NumeratorPuppi_vsPt_forward_1,DenominatorPuppi_vsPt_forward_1);

  frameEfficiency_vsPt->Delete();
  frameEfficiency_vsPt = new TH2F("frameEfficiency_vsPt","",30,jetPtThreshold,500,500,
				  min(EfficiencyReco_vsPt_forward_1->GetHistogram()->GetMinimum(),EfficiencyPuppi_vsPt_forward_1->GetHistogram()->GetMinimum()),
				  1.05*max(EfficiencyReco_vsPt_forward_1->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPt_forward_1->GetHistogram()->GetMaximum()));

  frameEfficiency_vsPt->SetLineWidth(2);
  frameEfficiency_vsPt->SetMarkerStyle(21);
  frameEfficiency_vsPt->SetMarkerSize(0.3);
  frameEfficiency_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  frameEfficiency_vsPt->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPt->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPt->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPt->SetStats(0);
  frameEfficiency_vsPt->Draw();

  EfficiencyReco_vsPt_forward_1->SetLineWidth(2);
  EfficiencyReco_vsPt_forward_1->SetLineColor(kBlue);
  EfficiencyReco_vsPt_forward_1->SetMarkerStyle(20);
  EfficiencyReco_vsPt_forward_1->SetMarkerColor(kBlue);
  EfficiencyReco_vsPt_forward_1->GetXaxis()->SetTitle("p_{T} (GeV)");
  EfficiencyReco_vsPt_forward_1->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPt_forward_1->Draw("p");

  EfficiencyPuppi_vsPt_forward_1->SetLineWidth(2);
  EfficiencyPuppi_vsPt_forward_1->SetLineColor(kRed);
  EfficiencyPuppi_vsPt_forward_1->SetMarkerStyle(24);
  EfficiencyPuppi_vsPt_forward_1->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPt_forward_1->GetXaxis()->SetTitle("#eta");
  EfficiencyPuppi_vsPt_forward_1->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPt_forward_1->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPt_forward_1,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPt_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* EfficiencyReco_vsPt_forward_2 = new TGraphAsymmErrors();
  EfficiencyReco_vsPt_forward_2->SetName("EfficiencyReco_vsPt_forward_2");

  TGraphAsymmErrors* EfficiencyPuppi_vsPt_forward_2 = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPt_forward_2->SetName("EfficiencyPuppi_vsPt_forward_2");

  EfficiencyReco_vsPt_forward_2->BayesDivide(NumeratorReco_vsPt_forward_2,DenominatorReco_vsPt_forward_2);
  EfficiencyPuppi_vsPt_forward_2->BayesDivide(NumeratorPuppi_vsPt_forward_2,DenominatorPuppi_vsPt_forward_2);

  frameEfficiency_vsPt->Delete();
  frameEfficiency_vsPt = new TH2F("frameEfficiency_vsPt","",30,jetPtThreshold,500,500,
				  min(EfficiencyReco_vsPt_forward_2->GetHistogram()->GetMinimum(),EfficiencyPuppi_vsPt_forward_2->GetHistogram()->GetMinimum()),
				  1.05*max(EfficiencyReco_vsPt_forward_2->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPt_forward_2->GetHistogram()->GetMaximum()));

  frameEfficiency_vsPt->SetLineWidth(2);
  frameEfficiency_vsPt->SetMarkerStyle(21);
  frameEfficiency_vsPt->SetMarkerSize(0.3);
  frameEfficiency_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  frameEfficiency_vsPt->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPt->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPt->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPt->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPt->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPt->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPt->SetStats(0);
  frameEfficiency_vsPt->Draw();

  EfficiencyReco_vsPt_forward_2->SetLineWidth(2);
  EfficiencyReco_vsPt_forward_2->SetLineColor(kBlue);
  EfficiencyReco_vsPt_forward_2->SetMarkerStyle(20);
  EfficiencyReco_vsPt_forward_2->SetMarkerColor(kBlue);
  EfficiencyReco_vsPt_forward_2->GetXaxis()->SetTitle("p_{T} (GeV)");
  EfficiencyReco_vsPt_forward_2->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPt_forward_2->Draw("p");

  EfficiencyPuppi_vsPt_forward_2->SetLineWidth(2);
  EfficiencyPuppi_vsPt_forward_2->SetLineColor(kRed);
  EfficiencyPuppi_vsPt_forward_2->SetMarkerStyle(24);
  EfficiencyPuppi_vsPt_forward_2->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPt_forward_2->GetXaxis()->SetTitle("#eta");
  EfficiencyPuppi_vsPt_forward_2->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPt_forward_2->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPt_forward_2,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPt_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPt_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //vs PU
  TGraphAsymmErrors* EfficiencyReco_vsPU_central = new TGraphAsymmErrors();
  EfficiencyReco_vsPU_central->SetName("EfficiencyReco_vsPU_central");

  TGraphAsymmErrors* EfficiencyPuppi_vsPU_central = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPU_central->SetName("EfficiencyPuppi_vsPU_central");

  EfficiencyReco_vsPU_central->BayesDivide(NumeratorReco_vsPU_central,DenominatorReco_vsPU_central);
  EfficiencyPuppi_vsPU_central->BayesDivide(NumeratorPuppi_vsPU_central,DenominatorPuppi_vsPU_central);


  TH2F* frameEfficiency_vsPU = new TH2F("frameEfficiency_vsPU","",50,PU-40,PU+40,500,0.7,
					1.1*max(EfficiencyReco_vsPU_central->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPU_central->GetHistogram()->GetMaximum()));
  frameEfficiency_vsPU->SetLineWidth(2);
  frameEfficiency_vsPU->SetMarkerStyle(21);
  frameEfficiency_vsPU->SetMarkerSize(0.3);
  frameEfficiency_vsPU->GetXaxis()->SetTitle("N_{PU}");
  frameEfficiency_vsPU->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPU->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPU->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPU->SetStats(0);
  frameEfficiency_vsPU->Draw();

  EfficiencyReco_vsPU_central->SetLineWidth(2);
  EfficiencyReco_vsPU_central->SetLineColor(kBlue);
  EfficiencyReco_vsPU_central->SetMarkerStyle(20);
  EfficiencyReco_vsPU_central->SetMarkerColor(kBlue);
  EfficiencyReco_vsPU_central->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyReco_vsPU_central->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPU_central->Draw("p");

  EfficiencyPuppi_vsPU_central->SetLineWidth(2);
  EfficiencyPuppi_vsPU_central->SetLineColor(kRed);
  EfficiencyPuppi_vsPU_central->SetMarkerStyle(24);
  EfficiencyPuppi_vsPU_central->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPU_central->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyPuppi_vsPU_central->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPU_central->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPU_central,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPU_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* EfficiencyReco_vsPU_forward_1 = new TGraphAsymmErrors();
  EfficiencyReco_vsPU_forward_1->SetName("EfficiencyReco_vsPU_forward_1");

  TGraphAsymmErrors* EfficiencyPuppi_vsPU_forward_1 = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPU_forward_1->SetName("EfficiencyPuppi_vsPU_forward_1");

  EfficiencyReco_vsPU_forward_1->BayesDivide(NumeratorReco_vsPU_forward_1,DenominatorReco_vsPU_forward_1);
  EfficiencyPuppi_vsPU_forward_1->BayesDivide(NumeratorPuppi_vsPU_forward_1,DenominatorPuppi_vsPU_forward_1);

  frameEfficiency_vsPU->Delete();
  frameEfficiency_vsPU = new TH2F("frameEfficiency_vsPU","",50,PU-40,PU+40,500,0.7,
				  1.1*max(EfficiencyReco_vsPU_forward_1->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPU_forward_1->GetHistogram()->GetMaximum()));

  frameEfficiency_vsPU->SetLineWidth(2);
  frameEfficiency_vsPU->SetMarkerStyle(21);
  frameEfficiency_vsPU->SetMarkerSize(0.3);
  frameEfficiency_vsPU->GetXaxis()->SetTitle("N_{PU}");
  frameEfficiency_vsPU->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPU->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPU->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPU->SetStats(0);
  frameEfficiency_vsPU->Draw();

  EfficiencyReco_vsPU_forward_1->SetLineWidth(2);
  EfficiencyReco_vsPU_forward_1->SetLineColor(kBlue);
  EfficiencyReco_vsPU_forward_1->SetMarkerStyle(20);
  EfficiencyReco_vsPU_forward_1->SetMarkerColor(kBlue);
  EfficiencyReco_vsPU_forward_1->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyReco_vsPU_forward_1->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPU_forward_1->Draw("p");

  EfficiencyPuppi_vsPU_forward_1->SetLineWidth(2);
  EfficiencyPuppi_vsPU_forward_1->SetLineColor(kRed);
  EfficiencyPuppi_vsPU_forward_1->SetMarkerStyle(24);
  EfficiencyPuppi_vsPU_forward_1->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPU_forward_1->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyPuppi_vsPU_forward_1->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPU_forward_1->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPU_forward_1,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPU_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* EfficiencyReco_vsPU_forward_2 = new TGraphAsymmErrors();
  EfficiencyReco_vsPU_forward_2->SetName("EfficiencyReco_vsPU_forward_2");

  TGraphAsymmErrors* EfficiencyPuppi_vsPU_forward_2 = new TGraphAsymmErrors();
  EfficiencyPuppi_vsPU_forward_2->SetName("EfficiencyPuppi_vsPU_forward_2");

  EfficiencyReco_vsPU_forward_2->BayesDivide(NumeratorReco_vsPU_forward_2,DenominatorReco_vsPU_forward_2);
  EfficiencyPuppi_vsPU_forward_2->BayesDivide(NumeratorPuppi_vsPU_forward_2,DenominatorPuppi_vsPU_forward_2);

  frameEfficiency_vsPU->Delete();
  frameEfficiency_vsPU = new TH2F("frameEfficiency_vsPU","",50,PU-40,PU+40,500,0.7,
				  1.1*max(EfficiencyReco_vsPU_forward_2->GetHistogram()->GetMaximum(),EfficiencyPuppi_vsPU_forward_2->GetHistogram()->GetMaximum()));

  frameEfficiency_vsPU->SetLineWidth(2);
  frameEfficiency_vsPU->SetMarkerStyle(21);
  frameEfficiency_vsPU->SetMarkerSize(0.3);
  frameEfficiency_vsPU->GetXaxis()->SetTitle("N_{PU}");
  frameEfficiency_vsPU->GetXaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetXaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetXaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetXaxis()->SetTitleOffset(1.10);
  frameEfficiency_vsPU->GetYaxis()->SetTitle("efficiency");
  frameEfficiency_vsPU->GetYaxis()->SetLabelOffset(0.012);
  frameEfficiency_vsPU->GetYaxis()->SetLabelSize(0.038);
  frameEfficiency_vsPU->GetYaxis()->SetTitleSize(0.05);
  frameEfficiency_vsPU->GetYaxis()->SetTitleOffset(1.18);
  frameEfficiency_vsPU->SetStats(0);
  frameEfficiency_vsPU->Draw();

  EfficiencyReco_vsPU_forward_2->SetLineWidth(2);
  EfficiencyReco_vsPU_forward_2->SetLineColor(kBlue);
  EfficiencyReco_vsPU_forward_2->SetMarkerStyle(20);
  EfficiencyReco_vsPU_forward_2->SetMarkerColor(kBlue);
  EfficiencyReco_vsPU_forward_2->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyReco_vsPU_forward_2->GetYaxis()->SetTitle("Efficiency");
  EfficiencyReco_vsPU_forward_2->Draw("p");

  EfficiencyPuppi_vsPU_forward_2->SetLineWidth(2);
  EfficiencyPuppi_vsPU_forward_2->SetLineColor(kRed);
  EfficiencyPuppi_vsPU_forward_2->SetMarkerStyle(24);
  EfficiencyPuppi_vsPU_forward_2->SetMarkerColor(kRed);
  EfficiencyPuppi_vsPU_forward_2->GetXaxis()->SetTitle("N_{PU}");
  EfficiencyPuppi_vsPU_forward_2->GetYaxis()->SetTitle("Efficiency");
  EfficiencyPuppi_vsPU_forward_2->Draw("psame");

  legend->AddEntry(EfficiencyReco_vsPU_forward_2,"CHS jets","l");
  legend->AddEntry(EfficiencyPuppi_vsPU_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Purity_vsPU_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();


  ///////////////////////////////////
  // purity plots ///////////////////
  //////////////////////////////////

  TGraphAsymmErrors* PurityReco_vsEta = new TGraphAsymmErrors();
  PurityReco_vsEta->SetName("PurityReco_vsEta");

  TGraphAsymmErrors* PurityPuppi_vsEta = new TGraphAsymmErrors();
  PurityPuppi_vsEta->SetName("PurityPuppi_vsEta");

  PurityReco_vsEta->BayesDivide(NumeratorReco_vsEta_purity,DenominatorReco_vsEta_purity);
  PurityPuppi_vsEta->BayesDivide(NumeratorPuppi_vsEta_purity,DenominatorPuppi_vsEta_purity);


  TH2F* framePurity_vsEta = new TH2F("framePurity_vsEta","",50,-4.5,4.5,500,0.6,
					 1.1*max(PurityReco_vsEta->GetHistogram()->GetMaximum(),PurityPuppi_vsEta->GetHistogram()->GetMaximum()));
  framePurity_vsEta->SetLineWidth(2);
  framePurity_vsEta->SetMarkerStyle(21);
  framePurity_vsEta->SetMarkerSize(0.3);
  framePurity_vsEta->GetXaxis()->SetTitle("#eta");
  framePurity_vsEta->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsEta->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsEta->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsEta->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsEta->GetYaxis()->SetTitle("purity");
  framePurity_vsEta->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsEta->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsEta->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsEta->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsEta->SetStats(0);
  framePurity_vsEta->Draw();

  PurityReco_vsEta->SetLineWidth(2);
  PurityReco_vsEta->SetLineColor(kBlue);
  PurityReco_vsEta->SetMarkerStyle(20);
  PurityReco_vsEta->SetMarkerColor(kBlue);
  PurityReco_vsEta->GetXaxis()->SetTitle("#eta");
  PurityReco_vsEta->GetYaxis()->SetTitle("Purity");
  PurityReco_vsEta->Draw("p");

  PurityPuppi_vsEta->SetLineWidth(2);
  PurityPuppi_vsEta->SetLineColor(kRed);
  PurityPuppi_vsEta->SetMarkerStyle(24);
  PurityPuppi_vsEta->SetMarkerColor(kRed);
  PurityPuppi_vsEta->GetXaxis()->SetTitle("#eta");
  PurityPuppi_vsEta->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsEta->Draw("psame");

  legend->AddEntry(PurityReco_vsEta,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsEta,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsEta_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* PurityReco_vsPt_central = new TGraphAsymmErrors();
  PurityReco_vsPt_central->SetName("PurityReco_vsPt_central");

  TGraphAsymmErrors* PurityPuppi_vsPt_central = new TGraphAsymmErrors();
  PurityPuppi_vsPt_central->SetName("PurityPuppi_vsPt_central");

  PurityReco_vsPt_central->BayesDivide(NumeratorReco_vsPt_central_purity,DenominatorReco_vsPt_central_purity);
  PurityPuppi_vsPt_central->BayesDivide(NumeratorPuppi_vsPt_central_purity,DenominatorPuppi_vsPt_central_purity);

  
  TH2F* framePurity_vsPt = new TH2F("framePurity_vsPt","",50,jetPtThreshold,500,500,
					min(PurityReco_vsPt_central->GetHistogram()->GetMinimum(),PurityPuppi_vsPt_central->GetHistogram()->GetMinimum()),
					1.05*max(PurityReco_vsPt_central->GetHistogram()->GetMaximum(),PurityPuppi_vsPt_central->GetHistogram()->GetMaximum()));
  framePurity_vsPt->SetLineWidth(2);
  framePurity_vsPt->SetMarkerStyle(21);
  framePurity_vsPt->SetMarkerSize(0.3);
  framePurity_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  framePurity_vsPt->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPt->GetYaxis()->SetTitle("purity");
  framePurity_vsPt->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPt->SetStats(0);
  framePurity_vsPt->Draw();

  PurityReco_vsPt_central->SetLineWidth(2);
  PurityReco_vsPt_central->SetLineColor(kBlue);
  PurityReco_vsPt_central->SetMarkerStyle(20);
  PurityReco_vsPt_central->SetMarkerColor(kBlue);
  PurityReco_vsPt_central->GetXaxis()->SetTitle("p_{T} (GeV)");
  PurityReco_vsPt_central->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPt_central->Draw("p");

  PurityPuppi_vsPt_central->SetLineWidth(2);
  PurityPuppi_vsPt_central->SetLineColor(kRed);
  PurityPuppi_vsPt_central->SetMarkerStyle(24);
  PurityPuppi_vsPt_central->SetMarkerColor(kRed);
  PurityPuppi_vsPt_central->GetXaxis()->SetTitle("p_{T} (GeV)");
  PurityPuppi_vsPt_central->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPt_central->Draw("psame");

  legend->AddEntry(PurityReco_vsPt_central,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPt_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* PurityReco_vsPt_forward_1 = new TGraphAsymmErrors();
  PurityReco_vsPt_forward_1->SetName("PurityReco_vsPt_forward_1");

  TGraphAsymmErrors* PurityPuppi_vsPt_forward_1 = new TGraphAsymmErrors();
  PurityPuppi_vsPt_forward_1->SetName("PurityPuppi_vsPt_forward_1");

  PurityReco_vsPt_forward_1->BayesDivide(NumeratorReco_vsPt_forward_purity_1,DenominatorReco_vsPt_forward_purity_1);
  PurityPuppi_vsPt_forward_1->BayesDivide(NumeratorPuppi_vsPt_forward_purity_1,DenominatorPuppi_vsPt_forward_purity_1);

  framePurity_vsPt->Delete();
  framePurity_vsPt = new TH2F("framePurity_vsPt","",50,jetPtThreshold,500,500,
				  min(PurityReco_vsPt_forward_1->GetHistogram()->GetMinimum(),PurityPuppi_vsPt_forward_1->GetHistogram()->GetMinimum()),
				  1.05*max(PurityReco_vsPt_forward_1->GetHistogram()->GetMaximum(),PurityPuppi_vsPt_forward_1->GetHistogram()->GetMaximum()));

  framePurity_vsPt->SetLineWidth(2);
  framePurity_vsPt->SetMarkerStyle(21);
  framePurity_vsPt->SetMarkerSize(0.3);
  framePurity_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  framePurity_vsPt->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPt->GetYaxis()->SetTitle("purity");
  framePurity_vsPt->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPt->SetStats(0);
  framePurity_vsPt->Draw();

  PurityReco_vsPt_forward_1->SetLineWidth(2);
  PurityReco_vsPt_forward_1->SetLineColor(kBlue);
  PurityReco_vsPt_forward_1->SetMarkerStyle(20);
  PurityReco_vsPt_forward_1->SetMarkerColor(kBlue);
  PurityReco_vsPt_forward_1->GetXaxis()->SetTitle("p_{T} (GeV)");
  PurityReco_vsPt_forward_1->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPt_forward_1->Draw("p");

  PurityPuppi_vsPt_forward_1->SetLineWidth(2);
  PurityPuppi_vsPt_forward_1->SetLineColor(kRed);
  PurityPuppi_vsPt_forward_1->SetMarkerStyle(24);
  PurityPuppi_vsPt_forward_1->SetMarkerColor(kRed);
  PurityPuppi_vsPt_forward_1->GetXaxis()->SetTitle("#eta");
  PurityPuppi_vsPt_forward_1->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPt_forward_1->Draw("psame");

  legend->AddEntry(PurityReco_vsPt_forward_1,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPt_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* PurityReco_vsPt_forward_2 = new TGraphAsymmErrors();
  PurityReco_vsPt_forward_2->SetName("PurityReco_vsPt_forward_2");

  TGraphAsymmErrors* PurityPuppi_vsPt_forward_2 = new TGraphAsymmErrors();
  PurityPuppi_vsPt_forward_2->SetName("PurityPuppi_vsPt_forward_2");

  PurityReco_vsPt_forward_2->BayesDivide(NumeratorReco_vsPt_forward_purity_2,DenominatorReco_vsPt_forward_purity_2);
  PurityPuppi_vsPt_forward_2->BayesDivide(NumeratorPuppi_vsPt_forward_purity_2,DenominatorPuppi_vsPt_forward_purity_2);

  framePurity_vsPt->Delete();
  framePurity_vsPt = new TH2F("framePurity_vsPt","",50,jetPtThreshold,500,500,
				  min(PurityReco_vsPt_forward_2->GetHistogram()->GetMinimum(),PurityPuppi_vsPt_forward_2->GetHistogram()->GetMinimum()),
				  1.05*max(PurityReco_vsPt_forward_2->GetHistogram()->GetMaximum(),PurityPuppi_vsPt_forward_2->GetHistogram()->GetMaximum()));

  framePurity_vsPt->SetLineWidth(2);
  framePurity_vsPt->SetMarkerStyle(21);
  framePurity_vsPt->SetMarkerSize(0.3);
  framePurity_vsPt->GetXaxis()->SetTitle("p_{T} (GeV)");
  framePurity_vsPt->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPt->GetYaxis()->SetTitle("purity");
  framePurity_vsPt->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPt->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPt->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPt->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPt->SetStats(0);
  framePurity_vsPt->Draw();

  PurityReco_vsPt_forward_2->SetLineWidth(2);
  PurityReco_vsPt_forward_2->SetLineColor(kBlue);
  PurityReco_vsPt_forward_2->SetMarkerStyle(20);
  PurityReco_vsPt_forward_2->SetMarkerColor(kBlue);
  PurityReco_vsPt_forward_2->GetXaxis()->SetTitle("p_{T} (GeV)");
  PurityReco_vsPt_forward_2->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPt_forward_2->Draw("p");

  PurityPuppi_vsPt_forward_2->SetLineWidth(2);
  PurityPuppi_vsPt_forward_2->SetLineColor(kRed);
  PurityPuppi_vsPt_forward_2->SetMarkerStyle(24);
  PurityPuppi_vsPt_forward_2->SetMarkerColor(kRed);
  PurityPuppi_vsPt_forward_2->GetXaxis()->SetTitle("#eta");
  PurityPuppi_vsPt_forward_2->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPt_forward_2->Draw("psame");

  legend->AddEntry(PurityReco_vsPt_forward_2,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPt_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPt_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //vs PU
  TGraphAsymmErrors* PurityReco_vsPU_central = new TGraphAsymmErrors();
  PurityReco_vsPU_central->SetName("PurityReco_vsPU_central");

  TGraphAsymmErrors* PurityPuppi_vsPU_central = new TGraphAsymmErrors();
  PurityPuppi_vsPU_central->SetName("PurityPuppi_vsPU_central");

  PurityReco_vsPU_central->BayesDivide(NumeratorReco_vsPU_central_purity,DenominatorReco_vsPU_central_purity);
  PurityPuppi_vsPU_central->BayesDivide(NumeratorPuppi_vsPU_central_purity,DenominatorPuppi_vsPU_central_purity);


  TH2F* framePurity_vsPU = new TH2F("framePurity_vsPU","",50,PU-40,PU+40,500,0.7,
					1.1*max(PurityReco_vsPU_central->GetHistogram()->GetMaximum(),PurityPuppi_vsPU_central->GetHistogram()->GetMaximum()));
  framePurity_vsPU->SetLineWidth(2);
  framePurity_vsPU->SetMarkerStyle(21);
  framePurity_vsPU->SetMarkerSize(0.3);
  framePurity_vsPU->GetXaxis()->SetTitle("N_{PU}");
  framePurity_vsPU->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPU->GetYaxis()->SetTitle("purity");
  framePurity_vsPU->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPU->SetStats(0);
  framePurity_vsPU->Draw();

  PurityReco_vsPU_central->SetLineWidth(2);
  PurityReco_vsPU_central->SetLineColor(kBlue);
  PurityReco_vsPU_central->SetMarkerStyle(20);
  PurityReco_vsPU_central->SetMarkerColor(kBlue);
  PurityReco_vsPU_central->GetXaxis()->SetTitle("N_{PU}");
  PurityReco_vsPU_central->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPU_central->Draw("p");

  PurityPuppi_vsPU_central->SetLineWidth(2);
  PurityPuppi_vsPU_central->SetLineColor(kRed);
  PurityPuppi_vsPU_central->SetMarkerStyle(24);
  PurityPuppi_vsPU_central->SetMarkerColor(kRed);
  PurityPuppi_vsPU_central->GetXaxis()->SetTitle("N_{PU}");
  PurityPuppi_vsPU_central->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPU_central->Draw("psame");

  legend->AddEntry(PurityReco_vsPU_central,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPU_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* PurityReco_vsPU_forward_1 = new TGraphAsymmErrors();
  PurityReco_vsPU_forward_1->SetName("PurityReco_vsPU_forward_1");

  TGraphAsymmErrors* PurityPuppi_vsPU_forward_1 = new TGraphAsymmErrors();
  PurityPuppi_vsPU_forward_1->SetName("PurityPuppi_vsPU_forward_1");

  PurityReco_vsPU_forward_1->BayesDivide(NumeratorReco_vsPU_forward_purity_1,DenominatorReco_vsPU_forward_purity_1);
  PurityPuppi_vsPU_forward_1->BayesDivide(NumeratorPuppi_vsPU_forward_purity_1,DenominatorPuppi_vsPU_forward_purity_1);

  framePurity_vsPU->Delete();
  framePurity_vsPU = new TH2F("framePurity_vsPU","",50,PU-40,PU+40,500,0.7,
				  1.1*max(PurityReco_vsPU_forward_1->GetHistogram()->GetMaximum(),PurityPuppi_vsPU_forward_1->GetHistogram()->GetMaximum()));

  framePurity_vsPU->SetLineWidth(2);
  framePurity_vsPU->SetMarkerStyle(21);
  framePurity_vsPU->SetMarkerSize(0.3);
  framePurity_vsPU->GetXaxis()->SetTitle("N_{PU}");
  framePurity_vsPU->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPU->GetYaxis()->SetTitle("purity");
  framePurity_vsPU->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPU->SetStats(0);
  framePurity_vsPU->Draw();

  PurityReco_vsPU_forward_1->SetLineWidth(2);
  PurityReco_vsPU_forward_1->SetLineColor(kBlue);
  PurityReco_vsPU_forward_1->SetMarkerStyle(20);
  PurityReco_vsPU_forward_1->SetMarkerColor(kBlue);
  PurityReco_vsPU_forward_1->GetXaxis()->SetTitle("N_{PU}");
  PurityReco_vsPU_forward_1->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPU_forward_1->Draw("p");

  PurityPuppi_vsPU_forward_1->SetLineWidth(2);
  PurityPuppi_vsPU_forward_1->SetLineColor(kRed);
  PurityPuppi_vsPU_forward_1->SetMarkerStyle(24);
  PurityPuppi_vsPU_forward_1->SetMarkerColor(kRed);
  PurityPuppi_vsPU_forward_1->GetXaxis()->SetTitle("N_{PU}");
  PurityPuppi_vsPU_forward_1->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPU_forward_1->Draw("psame");

  legend->AddEntry(PurityReco_vsPU_forward_1,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPU_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // efficiency vs pt
  TGraphAsymmErrors* PurityReco_vsPU_forward_2 = new TGraphAsymmErrors();
  PurityReco_vsPU_forward_2->SetName("PurityReco_vsPU_forward_2");

  TGraphAsymmErrors* PurityPuppi_vsPU_forward_2 = new TGraphAsymmErrors();
  PurityPuppi_vsPU_forward_2->SetName("PurityPuppi_vsPU_forward_2");

  PurityReco_vsPU_forward_2->BayesDivide(NumeratorReco_vsPU_forward_purity_2,DenominatorReco_vsPU_forward_purity_2);
  PurityPuppi_vsPU_forward_2->BayesDivide(NumeratorPuppi_vsPU_forward_purity_2,DenominatorPuppi_vsPU_forward_purity_2);

  framePurity_vsPU->Delete();
  framePurity_vsPU = new TH2F("framePurity_vsPU","",50,PU-40,PU+40,500,0.7,
				  1.1*max(PurityReco_vsPU_forward_2->GetHistogram()->GetMaximum(),PurityPuppi_vsPU_forward_2->GetHistogram()->GetMaximum()));

  framePurity_vsPU->SetLineWidth(2);
  framePurity_vsPU->SetMarkerStyle(21);
  framePurity_vsPU->SetMarkerSize(0.3);
  framePurity_vsPU->GetXaxis()->SetTitle("N_{PU}");
  framePurity_vsPU->GetXaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetXaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetXaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetXaxis()->SetTitleOffset(1.10);
  framePurity_vsPU->GetYaxis()->SetTitle("purity");
  framePurity_vsPU->GetYaxis()->SetLabelOffset(0.012);
  framePurity_vsPU->GetYaxis()->SetLabelSize(0.038);
  framePurity_vsPU->GetYaxis()->SetTitleSize(0.05);
  framePurity_vsPU->GetYaxis()->SetTitleOffset(1.18);
  framePurity_vsPU->SetStats(0);
  framePurity_vsPU->Draw();

  PurityReco_vsPU_forward_2->SetLineWidth(2);
  PurityReco_vsPU_forward_2->SetLineColor(kBlue);
  PurityReco_vsPU_forward_2->SetMarkerStyle(20);
  PurityReco_vsPU_forward_2->SetMarkerColor(kBlue);
  PurityReco_vsPU_forward_2->GetXaxis()->SetTitle("N_{PU}");
  PurityReco_vsPU_forward_2->GetYaxis()->SetTitle("Purity");
  PurityReco_vsPU_forward_2->Draw("p");

  PurityPuppi_vsPU_forward_2->SetLineWidth(2);
  PurityPuppi_vsPU_forward_2->SetLineColor(kRed);
  PurityPuppi_vsPU_forward_2->SetMarkerStyle(24);
  PurityPuppi_vsPU_forward_2->SetMarkerColor(kRed);
  PurityPuppi_vsPU_forward_2->GetXaxis()->SetTitle("N_{PU}");
  PurityPuppi_vsPU_forward_2->GetYaxis()->SetTitle("Purity");
  PurityPuppi_vsPU_forward_2->Draw("psame");

  legend->AddEntry(PurityReco_vsPU_forward_2,"CHS jets","l");
  legend->AddEntry(PurityPuppi_vsPU_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/Efficiency_vsPU_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //Response plot
  ResponseReco_lowPt_central->SetLineWidth(2);
  ResponseReco_lowPt_central->SetLineColor(kBlue);
  ResponseReco_lowPt_central->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_lowPt_central->GetYaxis()->SetTitle("Entries");
  int integral = ResponseReco_lowPt_central->Integral();
  ResponseReco_lowPt_central->Scale(1./integral);
  ResponseReco_lowPt_central->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_lowPt_central->FindObject("stats");
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
  
  ResponsePuppi_lowPt_central->SetLineWidth(2);
  ResponsePuppi_lowPt_central->SetLineColor(kRed);
  integral = ResponsePuppi_lowPt_central->Integral();
  ResponsePuppi_lowPt_central->Scale(1./integral);
  ResponsePuppi_lowPt_central->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_lowPt_central->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_lowPt_central->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_lowPt_central->GetMaximum(),ResponsePuppi_lowPt_central->GetMaximum())*1.25);

  ResponseReco_lowPt_central->Draw("hist"); 
  ResponsePuppi_lowPt_central->SetLineStyle(7);  
  ResponsePuppi_lowPt_central->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_lowPt_central,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_lowPt_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  // ------------------
  ResponseReco_lowPt_forward_1->SetLineWidth(2);
  ResponseReco_lowPt_forward_1->SetLineColor(kBlue);
  ResponseReco_lowPt_forward_1->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_lowPt_forward_1->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_lowPt_forward_1->Integral();
  ResponseReco_lowPt_forward_1->Scale(1./integral);
  ResponseReco_lowPt_forward_1->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_lowPt_forward_1->FindObject("stats");
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
  
  ResponsePuppi_lowPt_forward_1->SetLineWidth(2);
  ResponsePuppi_lowPt_forward_1->SetLineColor(kRed);
  integral = ResponsePuppi_lowPt_forward_1->Integral();
  ResponsePuppi_lowPt_forward_1->Scale(1./integral);
  ResponsePuppi_lowPt_forward_1->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_lowPt_forward_1->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_lowPt_forward_1->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_lowPt_forward_1->GetMaximum(),ResponsePuppi_lowPt_forward_1->GetMaximum())*1.25);

  ResponseReco_lowPt_forward_1->Draw("hist");
  ResponsePuppi_lowPt_forward_1->SetLineStyle(7);  
  ResponsePuppi_lowPt_forward_1->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_lowPt_forward_1,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_lowPt_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();

  // ------------------
  ResponseReco_lowPt_forward_2->SetLineWidth(2);
  ResponseReco_lowPt_forward_2->SetLineColor(kBlue);
  ResponseReco_lowPt_forward_2->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_lowPt_forward_2->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_lowPt_forward_2->Integral();
  ResponseReco_lowPt_forward_2->Scale(1./integral);
  ResponseReco_lowPt_forward_2->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_lowPt_forward_2->FindObject("stats");
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
  
  ResponsePuppi_lowPt_forward_2->SetLineWidth(2);
  ResponsePuppi_lowPt_forward_2->SetLineColor(kRed);
  integral = ResponsePuppi_lowPt_forward_2->Integral();
  ResponsePuppi_lowPt_forward_2->Scale(1./integral);
  ResponsePuppi_lowPt_forward_2->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_lowPt_forward_2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_lowPt_forward_2->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_lowPt_forward_2->GetMaximum(),ResponsePuppi_lowPt_forward_2->GetMaximum())*1.25);

  ResponseReco_lowPt_forward_2->Draw("hist");
  ResponsePuppi_lowPt_forward_2->SetLineStyle(7);  
  ResponsePuppi_lowPt_forward_2->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_lowPt_forward_2,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_lowPt_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_lowPt_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();

  legend->Clear();


  // ----------------
  ResponseReco_highPt_central->SetLineWidth(2);
  ResponseReco_highPt_central->SetLineColor(kBlue);
  ResponseReco_highPt_central->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_highPt_central->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_highPt_central->Integral();
  ResponseReco_highPt_central->Scale(1./integral);
  ResponseReco_highPt_central->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_highPt_central->FindObject("stats");
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
  
  ResponsePuppi_highPt_central->SetLineWidth(2);
  ResponsePuppi_highPt_central->SetLineColor(kRed);
  integral = ResponsePuppi_highPt_central->Integral();
  ResponsePuppi_highPt_central->Scale(1./integral);
  ResponsePuppi_highPt_central->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_highPt_central->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_highPt_central->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_highPt_central->GetMaximum(),ResponsePuppi_highPt_central->GetMaximum())*1.25);

  ResponseReco_highPt_central->Draw("hist");
  ResponsePuppi_highPt_central->SetLineStyle(7);  
  ResponsePuppi_highPt_central->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_highPt_central,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_highPt_central,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_central_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //
  ResponseReco_highPt_forward_1->SetLineWidth(2);
  ResponseReco_highPt_forward_1->SetLineColor(kBlue);
  ResponseReco_highPt_forward_1->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_highPt_forward_1->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_highPt_forward_1->Integral();
  ResponseReco_highPt_forward_1->Scale(1./integral);
  ResponseReco_highPt_forward_1->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_highPt_forward_1->FindObject("stats");
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
  
  ResponsePuppi_highPt_forward_1->SetLineWidth(2);
  ResponsePuppi_highPt_forward_1->SetLineColor(kRed);
  integral = ResponsePuppi_highPt_forward_1->Integral();
  ResponsePuppi_highPt_forward_1->Scale(1./integral);
  ResponsePuppi_highPt_forward_1->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_highPt_forward_1->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_highPt_forward_1->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_highPt_forward_1->GetMaximum(),ResponsePuppi_highPt_forward_1->GetMaximum())*1.25);

  ResponseReco_highPt_forward_1->Draw("hist");
  ResponsePuppi_highPt_forward_1->SetLineStyle(7);  
  ResponsePuppi_highPt_forward_1->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_highPt_forward_1,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_highPt_forward_1,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_1_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();


  //
  ResponseReco_highPt_forward_2->SetLineWidth(2);
  ResponseReco_highPt_forward_2->SetLineColor(kBlue);
  ResponseReco_highPt_forward_2->GetXaxis()->SetTitle("(p_{T}^{Reco}-p_{T}^{Gen})/p_{T}^{Gen}");
  ResponseReco_highPt_forward_2->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_highPt_forward_2->Integral();
  ResponseReco_highPt_forward_2->Scale(1./integral);
  ResponseReco_highPt_forward_2->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_highPt_forward_2->FindObject("stats");
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
  
  ResponsePuppi_highPt_forward_2->SetLineWidth(2);
  ResponsePuppi_highPt_forward_2->SetLineColor(kRed);
  integral = ResponsePuppi_highPt_forward_2->Integral();
  ResponsePuppi_highPt_forward_2->Scale(1./integral);
  ResponsePuppi_highPt_forward_2->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_highPt_forward_2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_highPt_forward_2->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_highPt_forward_2->GetMaximum(),ResponsePuppi_highPt_forward_2->GetMaximum())*1.25);

  ResponseReco_highPt_forward_2->Draw("hist");
  ResponsePuppi_highPt_forward_2->SetLineStyle(7);  
  ResponsePuppi_highPt_forward_2->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_highPt_forward_2,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_highPt_forward_2,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_highPt_forward_2_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //
  ResponseReco_deta->SetLineWidth(2);
  ResponseReco_deta->SetLineColor(kBlue);
  ResponseReco_deta->GetXaxis()->SetTitle("(#Delta#eta_{jj}^{Reco}-#Delta#eta_{jj}^{Gen})/#Delta#eta_{jj}^{Gen}");
  ResponseReco_deta->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_deta->Integral();
  ResponseReco_deta->Scale(1./integral);
 
  ResponseReco_deta->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_deta->FindObject("stats");
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
  
  ResponsePuppi_deta->SetLineWidth(2);
  ResponsePuppi_deta->SetLineColor(kRed);
  integral = ResponsePuppi_deta->Integral();
  ResponsePuppi_deta->Scale(1./integral);
  ResponsePuppi_deta->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_deta->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_deta->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_deta->GetMaximum(),ResponsePuppi_deta->GetMaximum())*1.25);

  ResponseReco_deta->Draw("hist");
  ResponsePuppi_deta->SetLineStyle(7);
  ResponsePuppi_deta->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_deta,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_deta,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_deta_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  //
  ResponseReco_mjj->SetLineWidth(2);
  ResponseReco_mjj->SetLineColor(kBlue);
  ResponseReco_mjj->GetXaxis()->SetTitle("(M_{jj}^{Reco}-M_{jj}^{Gen})/M_{jj}^{Gen}");
  ResponseReco_mjj->GetYaxis()->SetTitle("Entries");
  integral = ResponseReco_mjj->Integral();
  ResponseReco_mjj->Scale(1./integral);
  ResponseReco_mjj->Draw("hist");

  gPad->Update();
  tps1 = (TPaveStats*) ResponseReco_mjj->FindObject("stats");
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
  
  ResponsePuppi_mjj->SetLineWidth(2);
  ResponsePuppi_mjj->SetLineColor(kRed);
  integral = ResponsePuppi_mjj->Integral();
  ResponsePuppi_mjj->Scale(1./integral);
  ResponsePuppi_mjj->Draw("hist");
  
  gPad->Update();
  tps2 = (TPaveStats*) ResponsePuppi_mjj->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetLineColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);

  ResponseReco_mjj->GetYaxis()->SetRangeUser(0.001,max(ResponseReco_mjj->GetMaximum(),ResponsePuppi_mjj->GetMaximum())*1.25);

  ResponseReco_mjj->Draw("hist");
  ResponsePuppi_mjj->SetLineStyle(7);
  ResponsePuppi_mjj->Draw("hist same");

  tps1->SetFillStyle(0);
  tps2->SetFillStyle(0);

  tps1->Draw("same");
  tps2->Draw("same");

  legend->AddEntry(ResponseReco_mjj,"CHS jets","l");
  legend->AddEntry(ResponsePuppi_mjj,"Puppi jets","l");

  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  legend->Draw("same");

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj.root").c_str(),"root");

  cCanvas->SetLogy();

  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj_log.pdf").c_str(),"pdf");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj_log.png").c_str(),"png");
  cCanvas->SaveAs(string(outputFileDirectory+"/ResponseReco_mjj_log.root").c_str(),"root");
  cCanvas->SetLogy(0);
  gPad->Update();


  legend->Clear();

  return 0;
    
}
