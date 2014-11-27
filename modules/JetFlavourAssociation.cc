#include "modules/JetFlavourAssociation.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm> 
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

class PartonClassifier : public ExRootClassifier{

public:

  PartonClassifier() {}
  Int_t GetCategory(TObject *object);  
  Double_t fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------
// As done: https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc

Int_t PartonClassifier::GetCategory(TObject *object){ // select parton in the parton list

  Candidate *parton = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1; // inside the eta and momentum range (be a little bit larger that the tracking coverage
  
  pdgCode = TMath::Abs(parton->PID);

  if(parton->Status == -1)         return -1;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton skip
  if(parton->Status == 3 or parton->Status == 2) return 0; // if status 3 return 

  return 0;
}

class LHEPartonClassifier : public ExRootClassifier{

public:

  LHEPartonClassifier() {}
  Int_t GetCategory(TObject *object);  
  Double_t fEtaMax, fPTMin;
};

Int_t LHEPartonClassifier::GetCategory(TObject *object){ // select parton in the parton list

  Candidate *parton = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1; // inside the eta and momentum range (be a little bit larger that the tracking coverage
  
  pdgCode = TMath::Abs(parton->PID);
  if(parton->Status == -1)         return -1 ;
  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton skip
  if(parton->Status != 1)          return -1; // if status 3 return 

  return 0;
}

//------------------------------------------------------------------------------

JetFlavourAssociation::JetFlavourAssociation() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0),fItLHEPartonInputArray(0),fItJetInputArray(0), fItParticleInputArray(0){
  fClassifier    = new PartonClassifier;
  fClassifierLHE = new LHEPartonClassifier;
}

//------------------------------------------------------------------------------

JetFlavourAssociation::~JetFlavourAssociation(){
  if(fClassifier) delete fClassifier;
  if(fClassifierLHE) delete fClassifierLHE;
}

//------------------------------------------------------------------------------

void JetFlavourAssociation::Init(){

  ExRootConfParam param;

  fDeltaR    = GetDouble("DeltaR", 0.5);

  fClassifier->fPTMin  = GetDouble("PartonPTMin", 0.);
  fClassifier->fEtaMax = GetDouble("PartonEtaMax",2.5);

  fClassifierLHE->fPTMin  = GetDouble("PartonPTMin", 0.);
  fClassifierLHE->fEtaMax = GetDouble("PartonEtaMax",2.5);

  // import input array(s)
  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fLHEPartonInputArray = ImportArray(GetString("LHEPartonInputArray", "Delphes/LHEParticles"));
  fItLHEPartonInputArray = fPartonInputArray->MakeIterator();

  fFilter    = new ExRootFilter(fPartonInputArray);
  fFilterLHE = new ExRootFilter(fLHEPartonInputArray);
  
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

}

//------------------------------------------------------------------------------
void JetFlavourAssociation::Finish(){

  if(fFilter)    delete fFilter;
  if(fFilterLHE) delete fFilterLHE;

  if(fItJetInputArray) delete fItJetInputArray;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;
  if(fItLHEPartonInputArray) delete fItLHEPartonInputArray;

}

//------------------------------------------------------------------------------

void JetFlavourAssociation::Process(){

  Candidate *jet;
  TObjArray *partonArray;
  TObjArray *LHEpartonArray;

  // select quark and gluons
  fFilter->Reset();
  partonArray = fFilter->GetSubArray(fClassifier, 0); // get the filtered parton array 

  if(partonArray == 0) return;
  TIter itPartonArray(partonArray);

  fFilterLHE->Reset();
  LHEpartonArray = fFilterLHE->GetSubArray(fClassifierLHE, 0); // get the filtered parton array 

  if(LHEpartonArray == 0) return;
  TIter itLHEPartonArray(LHEpartonArray);
  
  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next()))){
    // get standard flavour
    GetAlgoFlavour(jet,itPartonArray,itLHEPartonArray);
    GetPhysicsFlavour(jet,itPartonArray,itLHEPartonArray);
  }
}

//------------------------------------------------------------------------------
// Standard definition of jet flavour in https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?v=CMSSW_7_3_0_pre1

void JetFlavourAssociation::GetAlgoFlavour(Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray){
  
  Candidate tempParticle ;
  Candidate tempPartonHighestPt;
  Candidate tempNearest ;
  float maxPt = 0;
  float minDr = 1000;
  Candidate* parton, *LHEParton;
  int pdgCode, pdgCodeMax ;

  itPartonArray.Reset();
  while((parton = static_cast<Candidate*>(itPartonArray.Next()))){ 

    // default delphes method
    if(TMath::Abs(parton->PID) == 21) pdgCode = 0;
    if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR){
       if(pdgCodeMax < pdgCode) pdgCodeMax = pdgCode;
    }

    bool isGoodParton = true;
    itLHEPartonArray.Reset();
    while((LHEParton = static_cast<Candidate*>(itLHEPartonArray.Next()))){
      if(parton->Momentum.DeltaR(LHEParton->Momentum) < 0.001 and parton->PID == LHEParton->PID and LHEParton->Charge == parton->Charge){
         isGoodParton = false;
         break;
      }

      if(!isGoodParton) continue; 

      // check the daugheter 
      int ndaughters = 0;
      int nparton_daughters = 0;
      if(parton->D1 !=-1) ndaughters++;
      if(parton->D2 !=-1) ndaughters++;
      if(ndaughters > 0) { // partons are only quarks or gluons            
          int daughterFlavour1 = -1;
          int daughterFlavour2 = -1;          
          if(parton->D1 !=-1) daughterFlavour1 = TMath::Abs(dynamic_cast<Candidate*>(fParticleInputArray->At(parton->D1))->PID);
          if(parton->D2 !=-1) daughterFlavour2 = TMath::Abs(dynamic_cast<Candidate*>(fParticleInputArray->At(parton->D2))->PID);          
          if((daughterFlavour1 == 1 || daughterFlavour1 == 2 || daughterFlavour1 == 3 || daughterFlavour1 == 4 || daughterFlavour1 == 5 || daughterFlavour1 == 21)) nparton_daughters++;
          if((daughterFlavour2 == 1 || daughterFlavour2 == 2 || daughterFlavour2 == 3 || daughterFlavour2 == 4 || daughterFlavour1 == 5 || daughterFlavour2 == 21)) nparton_daughters++;
      }	
      if(nparton_daughters > 0) continue ;        
      if(jet->Momentum.DeltaR(parton->Momentum) <= fDeltaR){
	if(jet->Momentum.DeltaR(parton->Momentum) < minDr){
	    minDr       = jet->Momentum.DeltaR(parton->Momentum);
	    tempNearest = *parton; 
	}
	  
        if(TMath::Abs(parton->PID) == 4) tempParticle = *parton; // if not yet found and pdgId is a c, take as c
	if(TMath::Abs(parton->PID) == 5) tempParticle = *parton;
	if(parton->Momentum.Pt() > maxPt ) {
	     maxPt = parton->Momentum.Pt();
	     tempPartonHighestPt = *parton;
	}	  
      }
    }
  }
 
  jet->flavourHeaviest  = TMath::Abs(tempParticle.PID);
  jet->flavourHighestPt = TMath::Abs(tempPartonHighestPt.PID);
  jet->flavourNearest2  = TMath::Abs(tempNearest.PID);
  if(tempParticle.Momentum.Pt() <= 0) tempParticle = tempPartonHighestPt;
  jet->flavourAlgo      = TMath::Abs(tempParticle.PID);

      
  if(pdgCodeMax == 0)  pdgCodeMax = 21;
  if(pdgCodeMax == -1) pdgCodeMax = 0;

  jet->flavourDefault = pdgCodeMax;

}


void JetFlavourAssociation::GetPhysicsFlavour(Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray){

  Candidate tempParticle ;
  Candidate tempNearest;
  float     minDr = 1000;
  int       nInTheCone = 0;
  float     TheBiggerConeSize = 0.7;
  Candidate* parton, *LHEParton;
  std::vector<Candidate> contaminations;
  contaminations.clear();

  itLHEPartonArray.Reset();
  while((LHEParton = static_cast<Candidate*>(itLHEPartonArray.Next()))){ 

      float dist = jet->Momentum.DeltaR(LHEParton->Momentum); // take the DR
      if(LHEParton->Status == 1 and dist < minDr ){
	    tempNearest = *LHEParton; 
            minDr       = dist ;
      }

      if( LHEParton->Status == 1  and dist <= fDeltaR){
	     tempParticle = *LHEParton;
             nInTheCone++;
      }
  }

  itPartonArray.Reset();
  itLHEPartonArray.Reset();
  while((parton = static_cast<Candidate*>(itPartonArray.Next()))){ 
      float dist = jet->Momentum.DeltaR(parton->Momentum); // take the DR
      bool isGoodCandidate = true;
      while((LHEParton = static_cast<Candidate*>(itLHEPartonArray.Next()))){ 
	if(parton->Momentum.DeltaR(LHEParton->Momentum) < 0.01 and parton->PID == LHEParton->PID and LHEParton->Charge == parton->Charge){
            isGoodCandidate = false ;
	    break;
	}
      }
      if(!isGoodCandidate) continue;
      if(parton->D1 != -1 or parton->D2!=-1){         
        if((TMath::Abs(parton->PID) < 4 or TMath::Abs(parton->PID) == 21)) continue ;
        if(dist < TheBiggerConeSize) contaminations.push_back(*parton);
      }
  }
  

  jet->flavourNearest3 =  TMath::Abs(tempNearest.PID);
  if(nInTheCone != 1) jet->flavourPhysics = 0;
  else if(contaminations.size() == 0) jet->flavourPhysics = TMath::Abs(tempParticle.PID);
  else if(contaminations.size() != 0){ 
    jet->flavourPhysics = TMath::Abs(tempParticle.PID);
    for(size_t iPart = 0; iPart < contaminations.size() ; iPart++){
      int contaminatingFlavour = TMath::Abs(contaminations.at(iPart).PID);
      int numberOfMothers = 0;
      if(contaminations.at(iPart).M1 !=-1) numberOfMothers++;
      if(contaminations.at(iPart).M2 !=-1) numberOfMothers++;
      if( numberOfMothers > 0 && dynamic_cast<Candidate*>(fParticleInputArray->At(contaminations.at(iPart).M1))->Momentum.DeltaR(tempParticle.Momentum) < 0.001 ) continue; // mother is the initialParton --> OK
      if( TMath::Abs(tempParticle.PID) == 4 ) {
	 if( contaminatingFlavour == 4 ) continue; // keep association --> the initialParton is a c --> the contaminated parton is a c
	 jet->flavourPhysics = 0; // all the other cases reject!
         break;
      }
    }
  }
  
}
