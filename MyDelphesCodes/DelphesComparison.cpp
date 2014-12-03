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

// g++ -Wall -o DelphesComparison `root-config --glibs --libs --cflags` -lTreePlayer DelphesComparison.cpp

typedef std::pair<TH1F*,TH1F*> histoPair ;
typedef std::pair<std::string, histoPair > mapElement;

float deltaPhi (float a, float b){
  if(fabs(a-b) > TMath::Pi()) return 2*TMath::Pi()-fabs(a-b);
  else return fabs(a-b);

}

void getHistogram(std::map<std::string,histoPair> & map){

  // PileUp distribution <average 50 PU >
  map.insert(mapElement("NPU", histoPair(new TH1F("N_{PU}_1","",40,30,70),new TH1F("N_{PU}_2","",40,30,70))));  

  // Rho evaluation
  map.insert(mapElement("Rho",      histoPair(new TH1F("#rho_1","",30,0,100),new TH1F("#rho_2","",30,0,100))));
  map.insert(mapElement("RhoCentral",histoPair(new TH1F("#rho^{central}_1","",30,0,100),new TH1F("#rho^{central}_2","",30,0,100))));
  map.insert(mapElement("RhoForward", histoPair(new TH1F("#rho^{forward}_1","",30,0,100),new TH1F("#rho^{forward}_2","",30,0,100))));

  // GenParticle
  map.insert(mapElement("genParticlePT",  histoPair(new TH1F("genParticlePT_1","",50,0,50),  new TH1F("genParticlePT_2","",50,0,50))));
  map.insert(mapElement("genParticleEta", histoPair(new TH1F("genParticleEta_1","",50,-5,5), new TH1F("genParticleEta_2","",50,-5,5))));

  // puppi Particle
  map.insert(mapElement("puppiParticlePT",  histoPair(new TH1F("puppiParticlePT_1","",50,0,50),  new TH1F("puppiParticlePT_2","",50,0,50))));
  map.insert(mapElement("puppiParticleEta", histoPair(new TH1F("puppiParticleEta_1","",50,-5,5), new TH1F("puppiParticleEta_2","",50,-5,5))));
  
  // jet occupancy

  map.insert(mapElement("JetOccupancy", histoPair(new TH1F("JetOccupancy_1","",20,0,20),new TH1F("JetOccupancy_2","",20,0,20))));
  map.insert(mapElement("JetOccupancyCentral", histoPair(new TH1F("JetOccupancyCentral_1","",20,0,20),new TH1F("JetOccupancyCentral_2","",20,0,20))));
  map.insert(mapElement("JetOccupancyForward", histoPair(new TH1F("JetOccupancyForward_1","",20,0,20),new TH1F("JetOccupancyForward_2","",20,0,20))));
  
  // leading jet
  map.insert(mapElement("JetPtLead",   histoPair(new TH1F("JetPtLead_1","",50,10,600), new TH1F("JetPtLead_2","",50,10,600))));
  map.insert(mapElement("JetEtaLead",  histoPair(new TH1F("JetEtaLead_1","",50,-5,5),  new TH1F("JetEtaLead_2","",50,-5,5))));
  map.insert(mapElement("JetMassLead", histoPair(new TH1F("JetMassLead_1","",30,0,150),new TH1F("JetMassLead_2","",30,0,150))));

  // second jet
  map.insert(mapElement("JetPtSecond",   histoPair(new TH1F("JetPtSecond_1","",50,10,600), new TH1F("JetPtSecond_2","",50,10,600))));
  map.insert(mapElement("JetEtaSecond",  histoPair(new TH1F("JetEtaSecond_1","",50,-5,5),  new TH1F("JetEtaSecond_2","",50,-5,5))));
  map.insert(mapElement("JetMassSecond", histoPair(new TH1F("JetMassSecond_1","",30,0,150),new TH1F("JetMassSecond_2","",30,0,150))));

  // third jet
  map.insert(mapElement("JetPtThird",  histoPair(new TH1F("JetPtThird_1","",50,10,600),  new TH1F("JetPtThird_2","",50,10,600))));
  map.insert(mapElement("JetEtaThird", histoPair(new TH1F("JetEtaThird_1","",50,-5,5),   new TH1F("JetEtaThird_2","",50,-5,5))));
  map.insert(mapElement("JetMassThird",histoPair(new TH1F("JetMassThird_1","",30,0,150), new TH1F("JetMassThird_2","",30,0,150))));

  // jet occupancy
  map.insert(mapElement("PuppiJetOccupancy",  histoPair(new TH1F("PuppiJetOccupancy_1","",20,0,20),new TH1F("PuppiJetOccupancy_2","",20,0,20))));
  map.insert(mapElement("PuppiJetOccupancyCentral", histoPair(new TH1F("PuppiJetOccupancyCentral_1","",20,0,20),new TH1F("PuppiJetOccupancyCentral_2","",20,0,20))));
  map.insert(mapElement("PuppiJetOccupancyForward", histoPair(new TH1F("PuppiJetOccupancyForward_1","",20,0,20),new TH1F("PuppiJetOccupancyForward_2","",20,0,20))));
  
  // leading puppi jet
  map.insert(mapElement("PuppiJetPtLead",  histoPair(new TH1F("PuppiJetPtLead_1","",30,10,600), new TH1F("PuppiJetPtLead_2","",30,10,600))));
  map.insert(mapElement("PuppiJetEtaLead", histoPair(new TH1F("PuppiJetEtaLead_1","",50,-5,5),  new TH1F("PuppiJetEtaLead_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassLead",histoPair(new TH1F("PuppiJetMassLead_1","",30,0,150),new TH1F("PuppiJetMassLead_2","",30,0,150))));

  // second puppi jet
  map.insert(mapElement("PuppiJetPtSecond",  histoPair(new TH1F("PuppiJetPtSecond_1","",30,10,600), new TH1F("PuppiJetPtSecond_2","",30,10,600))));
  map.insert(mapElement("PuppiJetEtaSecond", histoPair(new TH1F("PuppiJetEtaSecond_1","",50,-5,5),  new TH1F("PuppiJetEtaSecond_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassSecond",histoPair(new TH1F("PuppiJetMassSecond_1","",30,0,150),new TH1F("PuppiJetMassSecond_2","",30,0,150))));

  // third puppi jet
  map.insert(mapElement("PuppiJetPtThird",  histoPair(new TH1F("PuppiJetPtThird_1","",30,10,600), new TH1F("PuppiJetPtThird_2","",30,10,600))));
  map.insert(mapElement("PuppiJetEtaThird", histoPair(new TH1F("PuppiJetEtaThird_1","",50,-5,5),  new TH1F("PuppiJetEtaThird_2","",50,-5,5))));
  map.insert(mapElement("PuppiJetMassThird",histoPair(new TH1F("PuppiJetMassThird_1","",30,0,150),new TH1F("PuppiJetMassThird_2","",30,0,150))));

  // missing et 
  map.insert(mapElement("GenMissingET",histoPair(   new TH1F("GenMissingET_1","",30,0,600),   new TH1F("GenMissingET_2","",30,0,600))));
  map.insert(mapElement("RecoMissingET",histoPair(  new TH1F("RecoMissingET_1","",30,0,600),  new TH1F("RecoMissingET_2","",30,0,600))));
  map.insert(mapElement("PuppiMissingET",histoPair( new TH1F("PuppiMissingET_1","",30,0,600), new TH1F("PuppiMissingET_2","",30,0,600))));

  // muon performance
  map.insert(mapElement("MuonPt",histoPair(  new TH1F("MuonPt_1","",50,10,600),   new TH1F("MuonPt_2","", 50,10,600))));
  map.insert(mapElement("MuonEta",histoPair( new TH1F("MuonEta_1","",35,-2.5,2.5), new TH1F("MuonEta_2","",35,-2.5,2.5))));

  // electron performance
  map.insert(mapElement("ElectronPt", histoPair( new TH1F("ElectronPt_1","",50,10,600),    new TH1F("ElectronPt_2","",50,10,600))));
  map.insert(mapElement("ElectronEta",histoPair( new TH1F("ElectronEta_1","",35,-2.5,2.5), new TH1F("ElectronEta_2","",35,-2.5,2.5))));  
}

void getResponseHistogram( std::map<std::string,histoPair> & map) {

  // met response
  map.insert(mapElement("MetResp",histoPair(new TH1F("MetResp_1","",50,-100,100),new TH1F("MetPtResp_2","",50,-100,100))));
  map.insert(mapElement("puppiMetResp",histoPair(new TH1F("puppiMetResp_1","",50,-100,100),new TH1F("puppiMetPtResp_2","",50,-100,100))));

  // number of jet repsonse
  map.insert(mapElement("JetOccupancyResp",histoPair(new TH1F("JetOccupancyResp_1","",20,0,1),new TH1F("JetOccupancyResp_2","",20,0,1))));
  map.insert(mapElement("JetOccupancyRespCentral",histoPair(new TH1F("JetOccupancyRespCentral_1","",20,0,1),new TH1F("JetOccupancyRespCentral_2","",20,0,1))));
  map.insert(mapElement("JetOccupancyRespForward",histoPair(new TH1F("JetOccupancyRespForward_1","",20,0,1),new TH1F("JetOccupancyRespForward_2","",20,0,1))));

  // leading jet pt response
  map.insert(mapElement("JetPtRespLead",histoPair(new TH1F("JetPtRespLead_1","",30,-60,60),new TH1F("JetPtRespLead_2","",30,-60,60))));
  map.insert(mapElement("JetEtaRespLead",histoPair(new TH1F("JetEtaRespLead_1","",25,-1,1),new TH1F("JetEtaRespLead_2","",25,-1,1))));
  map.insert(mapElement("JetMassRespLead",histoPair(new TH1F("JetMassRespLead_1","",20,-30,30),new TH1F("JetMassRespLead_2","",20,-30,30))));

  // second jet pt response
  map.insert(mapElement("JetPtRespSecond",histoPair(new TH1F("JetPtRespSecond_1","",30,-60,60),new TH1F("JetPtRespSecond_2","",30,-60,60))));
  map.insert(mapElement("JetEtaRespSecond",histoPair(new TH1F("JetEtaRespSecond_1","",25,-1,1),new TH1F("JetEtaRespSecond_2","",25,-1,1))));
  map.insert(mapElement("JetMassRespSecond",histoPair(new TH1F("JetMassRespSecond_1","",20,-30,30),new TH1F("JetMassRespSecond_2","",20,-30,30))));

  // third jet pt response
  map.insert(mapElement("JetPtRespThird",histoPair(new TH1F("JetPtRespThird_1","",30,-60,60),new TH1F("JetPtRespThird_2","",30,-60,60))));
  map.insert(mapElement("JetEtaRespThird",histoPair(new TH1F("JetEtaRespThird_1","",25,-1,1),new TH1F("JetEtaRespThird_2","",25,-1,1))));
  map.insert(mapElement("JetMassRespThird",histoPair(new TH1F("JetMassRespThird_1","",20,-30,30),new TH1F("JetMassRespThird_2","",20,-30,30))));

  // jet pt response (all the jets)
  map.insert(mapElement("JetPtResp",histoPair(new TH1F("JetPtResp_1","",30,-60,60),new TH1F("JetPtResp_2","",30,-60,60))));
  map.insert(mapElement("JetEtaResp",histoPair(new TH1F("JetEtaResp_1","",25,-1,1),new TH1F("JetEtaResp_2","",25,-1,1))));
  map.insert(mapElement("JetMassResp",histoPair(new TH1F("JetMassResp_1","",20,-30,30),new TH1F("JetMassResp_2","",20,-30,30))));

  // jet pt response (all central jets)
  map.insert(mapElement("JetPtRespCentral",histoPair(new TH1F("JetPtRespCentral_1","",30,-60,60),new TH1F("JetPtRespCentral_2","",30,-60,60))));
  map.insert(mapElement("JetEtaRespCentral",histoPair(new TH1F("JetEtaRespCentral_1","",25,-1,1),new TH1F("JetEtaRespCentral_2","",25,-1,1))));
  map.insert(mapElement("JetMassRespCentral",histoPair(new TH1F("JetMassRespCentral_1","",20,-30,30),new TH1F("JetMassRespCentral_2","",20,-30,30))));

  // jet pt response (all forward jets)
  map.insert(mapElement("JetPtRespForward",histoPair(new TH1F("JetPtRespForward_1","",30,-60,60),new TH1F("JetPtRespForward_2","",30,-60,60))));
  map.insert(mapElement("JetEtaRespForward",histoPair(new TH1F("JetEtaRespForward_1","",25,-1,1),new TH1F("JetEtaRespForward_2","",25,-1,1))));
  map.insert(mapElement("JetMassRespForward",histoPair(new TH1F("JetMassRespForward_1","",20,-30,30),new TH1F("JetMassRespForward_2","",20,-30,30))));

  // number of puppi jet repsonse
  map.insert(mapElement("PuppiJetOccupancyResp",histoPair(new TH1F("PuppiJetOccupancyResp_1","",20,0,1),new TH1F("PuppiJetOccupancyResp_2","",20,0,1))));
  map.insert(mapElement("PuppiJetOccupancyRespCentral",histoPair(new TH1F("PuppiJetOccupancyRespCentral_1","",20,0,1),new TH1F("PuppiJetOccupancyRespCentral_2","",20,0,1))));
  map.insert(mapElement("PuppiJetOccupancyRespForward",histoPair(new TH1F("PuppiJetOccupancyRespForward_1","",20,0,1),new TH1F("PuppiJetOccupancyRespForward_2","",20,0,1))));

  map.insert(mapElement("PuppiJetPtRespLead",histoPair(new TH1F("PuppiJetPtRespLead_1","",30,-60,60),  new TH1F("PuppiJetPtRespLead_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaRespLead",histoPair(new TH1F("PuppiJetEtaRespLead_1","",25,-1,1),    new TH1F("PuppiJetEtaRespLead_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassRespLead",histoPair(new TH1F("PuppiJetMassRespLead_1","",20,-30,30),new TH1F("PuppiJetMassRespLead_2","",20,-30,30))));

  map.insert(mapElement("PuppiJetPtRespSecond",histoPair(new TH1F("PuppiJetPtRespSecond_1","",30,-60,60),new TH1F("PuppiJetPtRespSecond_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaRespSecond",histoPair(new TH1F("PuppiJetEtaRespSecond_1","",25,-1,1),new TH1F("PuppiJetEtaRespSecond_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassRespSecond",histoPair(new TH1F("PuppiJetMassRespSecond_1","",20,-30,30),new TH1F("PuppiJetMassRespSecond_2","",20,-30,30))));

  map.insert(mapElement("PuppiJetPtRespThird",histoPair(new TH1F("PuppiJetPtRespThird_1","",30,-60,60),new TH1F("PuppiJetPtRespThird_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaRespThird",histoPair(new TH1F("PuppiJetEtaRespThird_1","",25,-1,1),new TH1F("PuppiJetEtaRespThird_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassRespThird",histoPair(new TH1F("PuppiJetMassRespThird_1","",20,-30,30),new TH1F("PuppiJetMassRespThird_2","",20,-30,30))));

  // jet pt response (all the jets)
  map.insert(mapElement("PuppiJetPtResp",histoPair(new TH1F("PuppiJetPtResp_1","",30,-60,60),new TH1F("PuppiJetPtResp_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaResp",histoPair(new TH1F("PuppiJetEtaResp_1","",25,-1,1),new TH1F("PuppiJetEtaResp_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassResp",histoPair(new TH1F("PuppiJetMassResp_1","",20,-30,30),new TH1F("PuppiJetMassResp_2","",20,-30,30))));

  // jet pt response (all central jets)
  map.insert(mapElement("PuppiJetPtRespCentral",histoPair(new TH1F("PuppiJetPtRespCentral_1","",30,-60,60),new TH1F("PuppiJetPtRespCentral_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaRespCentral",histoPair(new TH1F("PuppiJetEtaRespCentral_1","",25,-1,1),new TH1F("PuppiJetEtaRespCentral_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassRespCentral",histoPair(new TH1F("PuppiJetMassRespCentral_1","",20,-30,30),new TH1F("PuppiJetMassRespCentral_2","",20,-30,30))));

  // jet pt response (all forward jets)
  map.insert(mapElement("PuppiJetPtRespForward",histoPair(new TH1F("PuppiJetPtRespForward_1","",30,-60,60),new TH1F("PuppiJetPtRespForward_2","",30,-60,60))));
  map.insert(mapElement("PuppiJetEtaRespForward",histoPair(new TH1F("PuppiJetEtaRespForward_1","",25,-1,1),new TH1F("PuppiJetEtaRespForward_2","",25,-1,1))));
  map.insert(mapElement("PuppiJetMassRespForward",histoPair(new TH1F("PuppiJetMassRespForward_1","",20,-30,30),new TH1F("PuppiJetMassRespForward_2","",20,-30,30))));
}


/////////////////////////////////////////////////////

int main (int argc, char** argv){

  // Setting for style 
  std::string ROOTStyle;
  if(getenv ("ROOTStyle")!=NULL){
    ROOTStyle = getenv ("ROOTStyle");
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootLogon.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootPalette.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/rootColors.C").c_str());
    gROOT->ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C").c_str());
  }

  gStyle->SetOptStat(111110);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);
 
  float PtCut;
  (argc >=3) ? PtCut = atof(argv[2]) : PtCut = 0 ;  

  // output folders
  std::string outputFileDirectory    = "MyDelphesCodes/outputPlotsDist";
  std::string outputFileDirectoryRes = "MyDelphesCodes/outputPlotsRes";

  system(("mkdir -p "+outputFileDirectory).c_str());
  system(("rm -r "   +outputFileDirectory+"/*").c_str());

  system(("mkdir -p "+outputFileDirectoryRes).c_str());
  system(("rm -r "   +outputFileDirectoryRes+"/*").c_str());

  // input files
  std::vector<std::string> inputFileList_1 ;
  std::vector<std::string> inputFileList_2 ;

  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_0.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_1.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_2.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_3.root");
  //  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_4.root");
  /*  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_5.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_6.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_7.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_8.root");
  */
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_9.root");
  //inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_10.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_12.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_13.root");
  inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_14.root");
  //inputFileList_1.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test/outputtree_15.root");

  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_0.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_1.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_2.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_3.root");
  //  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_4.root");
  /*  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_5.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_6.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_7.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_8.root");
  */
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_9.root");
  //  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_10.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_12.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_13.root");
  inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_14.root");
  //inputFileList_2.push_back("root://eoscms.cern.ch//store/caf/user/rgerosa/TPSAMPLES_14TEV/DELPHES_TREES/SS_EWK_uvev_test_Seth/outputtree_15.root");

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

  TClonesArray* GenParticle_1 = new TClonesArray("GenParticle");
  TClonesArray* GenParticle_2 = new TClonesArray("GenParticle");
  chain_1->SetBranchAddress("GenParticles",&GenParticle_1);    
  chain_2->SetBranchAddress("GenParticles",&GenParticle_2);    

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
  
  // Event Loop
  if(argc >= 2)  maxEvents = maxEvents/atoi(argv[1]);

  for(int iEntry = 0; iEntry < maxEvents ; iEntry++){

    if(iEntry%1000 == 0) std::cout<<" reading entry "<<iEntry<<std::endl;

    // Read the entry
    chain_1->GetEntry(iEntry);
    chain_2->GetEntry(iEntry);

    // fill PU
    for( int i = 0; i < NPU_1->GetEntries() and i < NPU_2->GetEntries(); i++){
      if(histogramSingleVariables.find("NPU") != histogramSingleVariables.end()){
        histogramSingleVariables["NPU"].first->Fill(dynamic_cast<ScalarHT*>(NPU_1->At(i))->HT);
        histogramSingleVariables["NPU"].second->Fill(dynamic_cast<ScalarHT*>(NPU_2->At(i))->HT);
      }
    }

    // Rho Kt4
    for( int i = 0; i < Rho_1->GetEntries() and i < Rho_2->GetEntries() ; i++){
      histogramSingleVariables["Rho"].first->Fill(dynamic_cast<Rho*>(Rho_1->At(i))->Rho);
      histogramSingleVariables["Rho"].second->Fill(dynamic_cast<Rho*>(Rho_2->At(i))->Rho);

      if(fabs(dynamic_cast<Rho*>(Rho_1->At(i))->Edges[0]) <= 2.5 and fabs(dynamic_cast<Rho*>(Rho_1->At(i))->Edges[1]) <= 2.5)
	histogramSingleVariables["RhoCentral"].first->Fill(dynamic_cast<Rho*>(Rho_1->At(i))->Rho);
      else
	histogramSingleVariables["RhoForward"].first->Fill(dynamic_cast<Rho*>(Rho_1->At(i))->Rho);

      if(fabs(dynamic_cast<Rho*>(Rho_2->At(i))->Edges[0]) <= 2.5 and fabs(dynamic_cast<Rho*>(Rho_2->At(i))->Edges[1]) <= 2.5)
	histogramSingleVariables["RhoCentral"].second->Fill(dynamic_cast<Rho*>(Rho_2->At(i))->Rho);
      else
	histogramSingleVariables["RhoForward"].second->Fill(dynamic_cast<Rho*>(Rho_2->At(i))->Rho);

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
      histogramSingleVariables["puppiParticlePT"].first->Fill(dynamic_cast<GenParticle*>(puppiParticle_1->At(i))->PT);
      histogramSingleVariables["puppiParticleEta"].first->Fill(dynamic_cast<GenParticle*>(puppiParticle_1->At(i))->Eta);
    }
    for(int i = 0; i < puppiParticle_2->GetEntries() ; i++){
      histogramSingleVariables["puppiParticlePT"].second->Fill(dynamic_cast<GenParticle*>(puppiParticle_2->At(i))->PT);
      histogramSingleVariables["puppiParticleEta"].second->Fill(dynamic_cast<GenParticle*>(puppiParticle_2->At(i))->Eta);
    }

    // Gen Particles
    for(int i = 0; i < GenParticle_1->GetEntries() ; i++){
      if(dynamic_cast<GenParticle*>(GenParticle_1->At(i))->IsPU != 0) continue;
      histogramSingleVariables["genParticlePT"].first->Fill(dynamic_cast<GenParticle*>(GenParticle_1->At(i))->PT);
      histogramSingleVariables["genParticleEta"].first->Fill(dynamic_cast<GenParticle*>(GenParticle_1->At(i))->Eta);
    }
   
    for(int i = 0; i < GenParticle_2->GetEntries() ; i++){
      if(dynamic_cast<GenParticle*>(GenParticle_2->At(i))->IsPU != 0) continue;
      histogramSingleVariables["genParticlePT"].second->Fill(dynamic_cast<GenParticle*>(GenParticle_2->At(i))->PT);
      histogramSingleVariables["genParticleEta"].second->Fill(dynamic_cast<GenParticle*>(GenParticle_2->At(i))->Eta);
    }

    // JET Occupancy ; central forward
    float nJets_1 = 0,        nJets_2 = 0;
    float nJetsCentral_1 = 0, nJetsCentral_2 = 0;
    float nJetsForward_1 = 0, nJetsForward_2 = 0;
    
    for(int i = 0; i < Jets_1->GetEntries() ; i++){
      nJets_1++;
      if(dynamic_cast<Jet*>(Jets_1->At(i))->PT < PtCut) continue ;
      if(fabs(dynamic_cast<Jet*>(Jets_1->At(i))->Eta) < 2.5) nJetsCentral_1++;
      else nJetsForward_1++;

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
       histogramSingleVariables["JetPtThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass);
      }
     
    }

    histogramSingleVariables["JetOccupancy"].first->Fill(nJets_1);
    histogramSingleVariables["JetOccupancyCentral"].first->Fill(nJetsCentral_1);
    histogramSingleVariables["JetOccupancyForward"].first->Fill(nJetsForward_1);

    for(int i = 0; i < Jets_2->GetEntries() ; i++){
      if(dynamic_cast<Jet*>(Jets_2->At(i))->PT < PtCut) continue ;
      nJets_2++;
      if(fabs(dynamic_cast<Jet*>(Jets_2->At(i))->Eta) < 2.5) nJetsCentral_2++;
      else nJetsForward_2++;

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
       histogramSingleVariables["JetPtThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT);
       histogramSingleVariables["JetEtaThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta);
       histogramSingleVariables["JetMassThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass);
      }
    }

    histogramSingleVariables["JetOccupancy"].second->Fill(nJets_2);
    histogramSingleVariables["JetOccupancyCentral"].second->Fill(nJetsCentral_2);
    histogramSingleVariables["JetOccupancyForward"].second->Fill(nJetsForward_2);

    // PUPPI JET 1 
    // JET Occupancy ; central forward
    nJets_1 = 0,        nJets_2 = 0;
    nJetsCentral_1 = 0, nJetsCentral_2 = 0;
    nJetsForward_1 = 0, nJetsForward_2 = 0;
    
    for(int i = 0; i < puppiJets_1->GetEntries() ; i++){

      if(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT < PtCut) continue ;
      nJets_1++;
      if(fabs(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta) < 2.5) nJetsCentral_1++;
      else nJetsForward_1++;

      if(i==0){
       histogramSingleVariables["PuppiJetPtLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassLead"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["PuppiJetPtSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassSecond"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["PuppiJetPtThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassThird"].first->Fill(dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass);
      }
    }

    histogramSingleVariables["PuppiJetOccupancy"].first->Fill(nJets_1);
    histogramSingleVariables["PuppiJetOccupancyCentral"].first->Fill(nJetsCentral_1);
    histogramSingleVariables["PuppiJetOccupancyForward"].first->Fill(nJetsForward_1);

    // PUPPI JET 2
    for(int i = 0; i < puppiJets_2->GetEntries() ; i++){
      if(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT < PtCut) continue ;
      nJets_2++;
      if(fabs(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta) < 2.5) nJetsCentral_2++;
      else nJetsForward_2++;

      if(i==0){
       histogramSingleVariables["PuppiJetPtLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassLead"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
      else if (i==1){
       histogramSingleVariables["PuppiJetPtSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassSecond"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
      else if (i==2){
       histogramSingleVariables["PuppiJetPtThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT);
       histogramSingleVariables["PuppiJetEtaThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta);
       histogramSingleVariables["PuppiJetMassThird"].second->Fill(dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass);
      }
    }

    histogramSingleVariables["PuppiJetOccupancy"].second->Fill(nJets_2);
    histogramSingleVariables["PuppiJetOccupancyCentral"].second->Fill(nJetsCentral_2);
    histogramSingleVariables["PuppiJetOccupancyForward"].second->Fill(nJetsForward_2);

    ///////////////////////////////////
    // RESPONSE
    /////////////////////////////////////

    //  MET
    for(int i = 0; i < RecoMET_1->GetEntries() and i < GenMET_1->GetEntries(); i++){
      if(histogramResponse.find("MetResp") == histogramResponse.end()) break ;
      histogramResponse["MetResp"].first->Fill(dynamic_cast<MissingET*>(RecoMET_1->At(i))->MET-dynamic_cast<MissingET*>(GenMET_1->At(i))->MET);
    }
    for(int i = 0; i < RecoMET_2->GetEntries() and i < GenMET_2->GetEntries(); i++){
      if(histogramResponse.find("MetResp") == histogramResponse.end()) break ;
      histogramResponse["MetResp"].second->Fill(dynamic_cast<MissingET*>(RecoMET_2->At(i))->MET-dynamic_cast<MissingET*>(GenMET_2->At(i))->MET);
    }

    for(int i = 0; i < PuppiMET_1->GetEntries() and i < GenMET_1->GetEntries(); i++){
      if(histogramResponse.find("puppiMetResp")  == histogramResponse.end()) break ;
      histogramResponse["puppiMetResp"].first->Fill(dynamic_cast<MissingET*>(PuppiMET_1->At(i))->MET-dynamic_cast<MissingET*>(GenMET_1->At(i))->MET);
    }
    for(int i = 0; i < PuppiMET_2->GetEntries() and i < GenMET_2->GetEntries(); i++){
      if(histogramResponse.find("puppiMetResp") == histogramResponse.end()) break ;
      histogramResponse["puppiMetResp"].second->Fill(dynamic_cast<MissingET*>(PuppiMET_2->At(i))->MET-dynamic_cast<MissingET*>(GenMET_2->At(i))->MET);
    }

    ////////// JET 

    nJets_1 = 0,        nJets_2 = 0;
    nJetsCentral_1 = 0, nJetsCentral_2 = 0;
    nJetsForward_1 = 0, nJetsForward_2 = 0;

    for(int i = 0; i < Jets_1->GetEntries() ; i++){
      if(dynamic_cast<Jet*>(Jets_1->At(i))->PT < PtCut) continue ;
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
	  histogramResponse["JetPtRespThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
      
        nJets_1++;
        if(fabs(dynamic_cast<Jet*>(Jets_1->At(i))->Eta) < 2.5){
          nJetsCentral_1++;
          histogramResponse["JetPtRespCentral"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
          histogramResponse["JetEtaRespCentral"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
          histogramResponse["JetMassRespCentral"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else{
          nJetsForward_1++;
          histogramResponse["JetPtRespForward"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
          histogramResponse["JetEtaRespForward"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
          histogramResponse["JetMassRespForward"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);

	}
        histogramResponse["JetPtResp"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->PT-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
        histogramResponse["JetEtaResp"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Eta-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
        histogramResponse["JetMassResp"].first->Fill(dynamic_cast<Jet*>(Jets_1->At(i))->Mass-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);

      }
    }

    histogramResponse["JetOccupancyResp"].first->Fill(float(nJets_1/GenJets_1->GetEntries()));
    histogramResponse["JetOccupancyRespCentral"].first->Fill(float(nJetsCentral_1/GenJets_1->GetEntries()));
    histogramResponse["JetOccupancyRespForward"].first->Fill(float(nJetsForward_1/GenJets_1->GetEntries()));



    for(int i = 0; i < Jets_2->GetEntries() ; i++){
      if(dynamic_cast<Jet*>(Jets_2->At(i))->PT < PtCut) continue ;
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
	  histogramResponse["JetPtRespThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["JetEtaRespThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["JetMassRespThird"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}

        nJets_2++;
        if(fabs(dynamic_cast<Jet*>(Jets_2->At(i))->Eta) < 2.5){ nJetsCentral_2++;
         histogramResponse["JetPtRespCentral"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
         histogramResponse["JetEtaRespCentral"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
         histogramResponse["JetMassRespCentral"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else{
          nJetsForward_2++;
          histogramResponse["JetPtRespForward"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
          histogramResponse["JetEtaRespForward"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
          histogramResponse["JetMassRespForward"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}

        histogramResponse["JetPtResp"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->PT-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
        histogramResponse["JetEtaResp"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Eta-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
        histogramResponse["JetMassResp"].second->Fill(dynamic_cast<Jet*>(Jets_2->At(i))->Mass-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
      }
    }

    histogramResponse["JetOccupancyResp"].second->Fill(float(nJets_2/GenJets_2->GetEntries()));
    histogramResponse["JetOccupancyRespCentral"].second->Fill(float(nJetsCentral_2/GenJets_2->GetEntries()));
    histogramResponse["JetOccupancyRespForward"].second->Fill(float(nJetsForward_2/GenJets_2->GetEntries()));

    /////////// PUPPI JET 
    nJets_1 = 0,        nJets_2 = 0;
    nJetsCentral_1 = 0, nJetsCentral_2 = 0;
    nJetsForward_1 = 0, nJetsForward_2 = 0;

    for(int i = 0; i < puppiJets_1->GetEntries() ; i++){

      TLorentzVector jetVector, jetArea;
      jetVector.SetPtEtaPhiM(dynamic_cast<Jet*>(puppiJets_1->At(i))->PT,dynamic_cast<Jet*>(puppiJets_1->At(i))->Eta,dynamic_cast<Jet*>(puppiJets_1->At(i))->Phi,dynamic_cast<Jet*>(puppiJets_1->At(i))->Mass); 
      jetArea.SetXYZT(dynamic_cast<Jet*>(puppiJets_1->At(i))->AreaX,dynamic_cast<Jet*>(puppiJets_1->At(i))->AreaY,dynamic_cast<Jet*>(puppiJets_1->At(i))->AreaZ,dynamic_cast<Jet*>(puppiJets_1->At(i))->AreaT);

      float RhoValue = 0; 

      for( int i = 0; i < puppiRho_1->GetEntries() ; i++){
	if(fabs(jetVector.Eta()) > fabs(dynamic_cast<Rho*>(puppiRho_1->At(i))->Edges[0]) and fabs(jetVector.Eta()) <= fabs(dynamic_cast<Rho*>(puppiRho_1->At(i))->Edges[1])){
	  RhoValue = dynamic_cast<Rho*>(puppiRho_1->At(i))->Rho;
          break;
	}
      }

      jetVector = jetVector + RhoValue*jetArea; 
      if(jetVector.Pt() < PtCut) continue ;

      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_1->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(j))->Eta),2)+pow(deltaPhi(jetVector.Phi(),dynamic_cast<Jet*>(GenJets_1->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["PuppiJetPtRespLead"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespLead"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespLead"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["PuppiJetPtRespSecond"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespSecond"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespSecond"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["PuppiJetPtRespThird"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespThird"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespThird"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}

        nJets_1++;
        if(fabs(jetVector.Eta()) < 2.5){
         nJetsCentral_1++;
         histogramResponse["PuppiJetPtRespCentral"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
         histogramResponse["PuppiJetEtaRespCentral"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
         histogramResponse["PuppiJetMassRespCentral"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}
        else{
          nJetsForward_1++;
          histogramResponse["PuppiJetPtRespForward"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
          histogramResponse["PuppiJetEtaRespForward"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
          histogramResponse["PuppiJetMassRespForward"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);
	}

        histogramResponse["PuppiJetPtResp"].first->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->PT);
        histogramResponse["PuppiJetEtaResp"].first->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Eta);
        histogramResponse["PuppiJetMassResp"].first->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_1->At(ijetMatched))->Mass);

      }
    }

    histogramResponse["PuppiJetOccupancyResp"].first->Fill(float(nJets_1/GenJets_1->GetEntries()));
    histogramResponse["PuppiJetOccupancyRespCentral"].first->Fill(float(nJetsCentral_1/GenJets_1->GetEntries()));
    histogramResponse["PuppiJetOccupancyRespForward"].first->Fill(float(nJetsForward_1/GenJets_1->GetEntries()));


    ////

    for(int i = 0; i < puppiJets_2->GetEntries() ; i++){

      TLorentzVector jetVector, jetArea;
      jetVector.SetPtEtaPhiM(dynamic_cast<Jet*>(puppiJets_2->At(i))->PT,dynamic_cast<Jet*>(puppiJets_2->At(i))->Eta,dynamic_cast<Jet*>(puppiJets_2->At(i))->Phi,dynamic_cast<Jet*>(puppiJets_2->At(i))->Mass); 
      jetArea.SetXYZT(dynamic_cast<Jet*>(puppiJets_2->At(i))->AreaX,dynamic_cast<Jet*>(puppiJets_2->At(i))->AreaY,dynamic_cast<Jet*>(puppiJets_2->At(i))->AreaZ,dynamic_cast<Jet*>(puppiJets_2->At(i))->AreaT);

      float RhoValue = 0; 

      for( int i = 0; i < puppiRho_2->GetEntries() ; i++){
	if(fabs(jetVector.Eta()) > fabs(dynamic_cast<Rho*>(puppiRho_2->At(i))->Edges[0]) and fabs(jetVector.Eta()) < fabs(dynamic_cast<Rho*>(puppiRho_2->At(i))->Edges[1]))
	  RhoValue = dynamic_cast<Rho*>(puppiRho_2->At(i))->Rho;
      }

      jetVector = jetVector + RhoValue*jetArea; 

      if(jetVector.Pt() < PtCut) continue ;

      float minDr = 9999 ;
      int ijetMatched = -1;
      for(int j = 0; j < GenJets_2->GetEntries(); j++){
        float dR = TMath::Sqrt(pow(fabs(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(j))->Eta),2)+pow(deltaPhi(jetVector.Phi(),dynamic_cast<Jet*>(GenJets_2->At(j))->Phi),2));
        if(dR < 0.3 and dR < minDr){
          minDr = dR ;      
          ijetMatched = j;
	}
      }

      if(minDr != 9999 and ijetMatched !=-1){
        if(i==0){
	  histogramResponse["PuppiJetPtRespLead"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespLead"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespLead"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==1){
	  histogramResponse["PuppiJetPtRespSecond"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespSecond"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespSecond"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else if(i==2){
	  histogramResponse["PuppiJetPtRespThird"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
	  histogramResponse["PuppiJetEtaRespThird"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
	  histogramResponse["PuppiJetMassRespThird"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}

        nJets_2++;
        if(fabs(jetVector.Eta()) < 2.5){
         nJetsCentral_2++;
         histogramResponse["PuppiJetPtRespCentral"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
         histogramResponse["PuppiJetEtaRespCentral"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
         histogramResponse["PuppiJetMassRespCentral"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}
        else{
          nJetsForward_2++;
          histogramResponse["PuppiJetPtRespForward"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
          histogramResponse["PuppiJetEtaRespForward"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
          histogramResponse["PuppiJetMassRespForward"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);
	}

        histogramResponse["PuppiJetPtResp"].second->Fill(jetVector.Pt()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->PT);
        histogramResponse["PuppiJetEtaResp"].second->Fill(jetVector.Eta()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Eta);
        histogramResponse["PuppiJetMassResp"].second->Fill(jetVector.M()-dynamic_cast<Jet*>(GenJets_2->At(ijetMatched))->Mass);

      }
    }

    histogramResponse["PuppiJetOccupancyResp"].second->Fill(float(nJets_2/GenJets_2->GetEntries()));
    histogramResponse["PuppiJetOccupancyRespCentral"].second->Fill(float(nJetsCentral_2/GenJets_2->GetEntries()));
    histogramResponse["PuppiJetOccupancyRespForward"].second->Fill(float(nJetsForward_2/GenJets_2->GetEntries()));
  }
  
  ///////////////////////////////////////////////
  // Plot
  ///////////////////////////////////////////////
  // make the plot vs PT
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

  for(std::map<std::string,histoPair>::const_iterator itMap = histogramSingleVariables.begin(); itMap != histogramSingleVariables.end(); itMap++){
    cCanvas->cd();
    legend->Clear();

    TString VariableName;
    VariableName.Form("%s",itMap->second.first->GetName());
    if(VariableName.Contains("_1"))  VariableName.ReplaceAll("_1","");

    if(VariableName.Contains("genParticle")) continue ;

    if(VariableName.Contains("puppiParticle")){
      itMap->second.first->GetYaxis()->SetRangeUser(0.001,std::max(itMap->second.first->GetMaximum(),itMap->second.second->GetMaximum())*10000);      
    }
    else{
     itMap->second.first->Scale(1/itMap->second.first->Integral());   
     itMap->second.second->Scale(1/itMap->second.second->Integral());   
     itMap->second.first->GetYaxis()->SetRangeUser(0.,itMap->second.first->GetMaximum()*2);
     cCanvas->SetLogy(0);
     gPad->Update();
    }
    
    itMap->second.first->SetLineColor(kBlack);
    itMap->second.first->SetMarkerColor(kBlack);
    itMap->second.first->SetLineWidth(2);

    itMap->second.first->GetXaxis()->SetTitle(VariableName);
    itMap->second.first->GetYaxis()->SetTitle("a.u.");
    itMap->second.first->SetMarkerStyle(20);
    
    itMap->second.first->Draw("hist");
    gPad->Update();
    TPaveStats *tps1 = (TPaveStats*) (itMap)->second.first->FindObject("stats");
    tps1->SetTextColor(kBlack);
    tps1->SetLineColor(kBlack);
    tps1->SetX1NDC(0.72);
    tps1->SetY1NDC(0.68);
    tps1->SetX2NDC(0.93);
    tps1->SetY2NDC(0.89);
    double X1 = tps1->GetX1NDC();
    double Y1 = tps1->GetY1NDC();
    double X2 = tps1->GetX2NDC();
    double Y2 = tps1->GetY2NDC();
    
    itMap->second.second->SetLineColor(kRed);
    itMap->second.second->SetMarkerColor(kRed);
    itMap->second.second->SetLineWidth(2);
    itMap->second.second->SetMarkerStyle(20);
    itMap->second.second->Draw("p");
    gPad->Update();
    itMap->second.first->Draw("hist");
    gPad->Update();
    itMap->second.second->Draw("psame");
    gPad->Update();

    TPaveStats *tps2 = (TPaveStats*) (itMap)->second.second->FindObject("stats");
    tps2->SetTextColor(kRed);
    tps2->SetLineColor(kRed);
    tps2->SetX1NDC(X1);
    tps2->SetX2NDC(X2);
    tps2->SetY1NDC(Y1-(Y2-Y1));
    tps2->SetY2NDC(Y1);
    tps1->SetFillStyle(0);
    tps2->SetFillStyle(0);
    tps1->Draw("same");
    tps2->Draw("same");
    

    if(not VariableName.Contains("puppiParticle")){
     double  chi2    = itMap->second.first->Chi2Test(itMap->second.second,"WW CHI2/NDF");
     double  KS_test = itMap->second.first->KolmogorovTest(itMap->second.second,"UO");
     TString probatext = Form("#chi^{2}/ndf = %0.2f K_{s} = %0.2f",float(chi2),float(KS_test));
     TLatex* tt = new TLatex(0.45,0.86,probatext);
     tt->SetNDC();
     tt->SetTextSize(0.025);
     tt->AppendPad("same");
    }

    legend->AddEntry(itMap->second.first, "New Delphes", "l");
    legend->AddEntry(itMap->second.second,"Old Delphes", "pl");

    if(VariableName.Contains("puppiParticlePT")){
      itMap->second.first->GetYaxis()->SetTitle("number of particles");
      std::map<std::string,histoPair>::iterator itGen = histogramSingleVariables.find("genParticlePT");           
      if(itGen != histogramSingleVariables.end()){
       if(itGen->second.first != NULL and itGen->second.first != 0){
         itGen->second.first->SetLineColor(kBlue);
         itGen->second.first->SetLineWidth(2);    
         itGen->second.first->Draw("histsame");
       }      
       cCanvas->SetLogy();
       gPad->Update();
       legend->AddEntry(itGen->second.first,"Gen Partlces","l");
     }
    }
    
    if(VariableName.Contains("puppiParticleEta")){
      itMap->second.first->GetYaxis()->SetTitle("number of particles");
      std::map<std::string,histoPair>::iterator itGen = histogramSingleVariables.find("genParticleEta");           
      if(itGen != histogramSingleVariables.end()){
       if(itGen->second.first != NULL and itGen->second.first != 0){
         itGen->second.first->SetLineColor(kBlue);
         itGen->second.first->SetLineWidth(2);    
         itGen->second.first->Draw("histsame");
       }      
       cCanvas->SetLogy();
       gPad->Update();
       legend->AddEntry(itGen->second.first,"Gen Partlces","l");
     }
    }
    

    tex->Draw("same");
    tex2->Draw("same");
    tex3->Draw("same");
    legend->Draw("same");

    VariableName.ReplaceAll("{","");
    VariableName.ReplaceAll("}","");
    VariableName.ReplaceAll("#","");
    VariableName.ReplaceAll("^","");

    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".pdf").c_str(),"pdf");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".png").c_str(),"png");
    cCanvas->SaveAs(std::string(outputFileDirectory+"/"+VariableName+".root").c_str(),"root");
    
  }   

  for(std::map<std::string,histoPair>::const_iterator itMap = histogramResponse.begin(); itMap != histogramResponse.end(); itMap++){

    cCanvas->cd();
    cCanvas->SetLogy(0);
    legend->Clear();

    TString VariableName;
    VariableName.Form("%s",itMap->second.first->GetName());
    if(VariableName.Contains("_1"))  VariableName.ReplaceAll("_1","");

    if(VariableName.Contains("puppiParticle")){
      itMap->second.first->GetYaxis()->SetRangeUser(0.001,std::max(itMap->second.first->GetMaximum(),itMap->second.second->GetMaximum())*1000);
    }
    else{
     itMap->second.first->Scale(1/itMap->second.first->Integral());   
     itMap->second.second->Scale(1/itMap->second.second->Integral());   
     itMap->second.first->GetYaxis()->SetRangeUser(0.,itMap->second.first->GetMaximum()*2);
     cCanvas->SetLogy(0);
     gPad->Update();
    }

    itMap->second.first->SetLineColor(kBlack);
    itMap->second.first->SetMarkerColor(kBlack);
    itMap->second.first->SetLineWidth(2);
    itMap->second.first->GetXaxis()->SetTitle(VariableName);
    itMap->second.first->GetYaxis()->SetTitle("a.u.");
    itMap->second.first->SetMarkerStyle(20);
    itMap->second.first->Draw("hist");
    gPad->Update();
    TPaveStats *tps1 = (TPaveStats*) (itMap)->second.first->FindObject("stats");
    tps1->SetTextColor(kBlack);
    tps1->SetLineColor(kBlack);
    tps1->SetX1NDC(0.72);
    tps1->SetY1NDC(0.68);
    tps1->SetX2NDC(0.93);
    tps1->SetY2NDC(0.89);
    double X1 = tps1->GetX1NDC();
    double Y1 = tps1->GetY1NDC();
    double X2 = tps1->GetX2NDC();
    double Y2 = tps1->GetY2NDC();

    itMap->second.second->SetLineColor(kRed);
    itMap->second.second->SetMarkerColor(kRed);
    itMap->second.second->SetLineWidth(2);
    itMap->second.second->SetMarkerStyle(20);
    itMap->second.second->Draw("p");
    gPad->Update();
    TPaveStats *tps2 = (TPaveStats*) (itMap)->second.second->FindObject("stats");
    tps2->SetTextColor(kRed);
    tps2->SetLineColor(kRed);
    tps2->SetX1NDC(X1);
    tps2->SetX2NDC(X2);
    tps2->SetY1NDC(Y1-(Y2-Y1));
    tps2->SetY2NDC(Y1);

    itMap->second.first->Draw("hist");
    gPad->Update();
    itMap->second.second->Draw("psame");
    gPad->Update();

    tps1->SetFillStyle(0);
    tps2->SetFillStyle(0);
    tps1->Draw("same");
    tps2->Draw("same");

    // chi2 and KS test for compatibility
    double  chi2    = itMap->second.first->Chi2Test(itMap->second.second,"WW CHI2/NDF");
    double  KS_test = itMap->second.first->KolmogorovTest(itMap->second.second,"UO");
    TString probatext = Form("#chi^{2}/ndf = %0.2f K_{s} = %0.2f",float(chi2),float(KS_test));
    TLatex* tt = new TLatex(0.45,0.86,probatext);
    tt->SetNDC();
    tt->SetTextSize(0.025);
    tt->AppendPad("same");

    tex->Draw("same");
    tex2->Draw("same");
    tex3->Draw("same");
 
    legend->AddEntry(itMap->second.first,"New Delphes","l");
    legend->AddEntry(itMap->second.second,"Old Delphes","pl");
    legend->Draw("same");
    gPad->Update();

    
    VariableName.ReplaceAll("{","");
    VariableName.ReplaceAll("}","");
    VariableName.ReplaceAll("#","");
    VariableName.ReplaceAll("^","");
    
    cCanvas->SaveAs(std::string(outputFileDirectoryRes+"/"+VariableName+".png").c_str(),"png");
    cCanvas->SaveAs(std::string(outputFileDirectoryRes+"/"+VariableName+".pdf").c_str(),"pdf");
    cCanvas->SaveAs(std::string(outputFileDirectoryRes+"/"+VariableName+".root").c_str(),"root");

  }

  return 0 ;

}
