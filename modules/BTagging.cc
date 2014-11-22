/** \class BTagging
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags 
 *  $Date: 2013-04-26 12:39:14 +0200 (Fri, 26 Apr 2013) $
 *  $Revision: 1099 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
*/

#include "modules/BTagging.h"

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

class BTaggingPartonClassifier : public ExRootClassifier{

public:

  BTaggingPartonClassifier() {}
  Int_t GetCategory(TObject *object);  
  Double_t fEtaMax, fPTMin;
};

//------------------------------------------------------------------------------
// As done: https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/PartonSelector.cc

Int_t BTaggingPartonClassifier::GetCategory(TObject *object){ // select parton in the parton list

  Candidate *parton = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = parton->Momentum;
  Int_t pdgCode;

  if(momentum.Pt() <= fPTMin || TMath::Abs(momentum.Eta()) > fEtaMax) return -1; // inside the eta and momentum range (be a little bit larger that the tracking coverage
  
  pdgCode = TMath::Abs(parton->PID);

  if(pdgCode != 21 && pdgCode > 5) return -1; // not a parton skip
  if(parton->Status == 3 or parton->Status == 2) return 0 ; // if status 3 return 

  return 0;
}

//------------------------------------------------------------------------------

BTagging::BTagging() :
  fClassifier(0), fFilter(0),
  fItPartonInputArray(0), fItJetInputArray(0), fItParticleInputArray(0){
  fClassifier = new BTaggingPartonClassifier;
}

//------------------------------------------------------------------------------

BTagging::~BTagging(){
  if(fClassifier) delete fClassifier;
}

//------------------------------------------------------------------------------

void BTagging::Init(){

  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  fBitNumber = GetInt("BitNumber", 0);
  fDeltaR    = GetDouble("DeltaR", 0.5);
  fClassifier->fPTMin  = GetDouble("PartonPTMin", 0.);
  fClassifier->fEtaMax = GetDouble("PartonEtaMax",2.5);

  // read efficiency formulas
  param = GetParam("EfficiencyFormula");
  size = param.GetSize();
  
  fEfficiencyMap.clear();
  for(i = 0; i < size/2; ++i){
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMap[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap = fEfficiencyMap.find(0);
  if(itEfficiencyMap == fEfficiencyMap.end()){
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMap[0] = formula;
  }

  // import input array(s)

  fPartonInputArray = ImportArray(GetString("PartonInputArray", "Delphes/partons"));
  fItPartonInputArray = fPartonInputArray->MakeIterator();

  fFilter = new ExRootFilter(fPartonInputArray);
  
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "Delphes/allParticles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

}

//------------------------------------------------------------------------------
void BTagging::Finish(){

  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  if(fFilter) delete fFilter;
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItPartonInputArray) delete fItPartonInputArray;
  for(itEfficiencyMap = fEfficiencyMap.begin(); itEfficiencyMap != fEfficiencyMap.end(); ++itEfficiencyMap){
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }
}

//------------------------------------------------------------------------------

void BTagging::Process(){

  Candidate *jet;
  TObjArray *partonArray;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  // select quark and gluons
  fFilter->Reset();
  partonArray = fFilter->GetSubArray(fClassifier, 0); // get the filtered parton array 
  
  if(partonArray == 0) return;
  TIter itPartonArray(partonArray);
  
  // loop over all input jets
  fItJetInputArray->Reset();
  while((jet = static_cast<Candidate*>(fItJetInputArray->Next()))){

    const TLorentzVector &jetMomentum = jet->Momentum;

    // get standard flavour
    GetAlgoFlavour(jet,itPartonArray);
    GetPhysicsFlavour(jet,itPartonArray);

    float randomNumber = gRandom->Uniform() ; // generate a common random number

    // find haviest flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourHeaviest);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagHeaviest |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;

    // find hisghestpt flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourHighestPt);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagHighestPt |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;

    // find nearest2 flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourNearest2);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagNearest2 |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;

    // find nearest3 flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourNearest3);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagNearest3 |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;

    // find nearest3 flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourAlgo);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagAlgo |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;

    // find nearest3 flavour and Btag
    itEfficiencyMap = fEfficiencyMap.find(jet->flavourPhysics);
    if(itEfficiencyMap == fEfficiencyMap.end()){
      itEfficiencyMap = fEfficiencyMap.find(0);
    }
    formula = itEfficiencyMap->second;
    jet->BTagPhysics |= (randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) << fBitNumber;
  }
}

//------------------------------------------------------------------------------
// Standard definition of jet flavour in https://cmssdt.cern.ch/SDT/lxr/source/PhysicsTools/JetMCAlgos/plugins/JetPartonMatcher.cc?v=CMSSW_7_3_0_pre1

void BTagging::GetAlgoFlavour(Candidate* jet, TIter & itPartonArray){
  
  Candidate tempParticle ;
  Candidate tempPartonHighestPt;
  Candidate tempNearest ;
  float maxPt = 0;
  float minDr = 1000;
  Candidate* parton;

  itPartonArray.Reset();
  while((parton = static_cast<Candidate*>(itPartonArray.Next()))){ 
      if(parton->Status != 3){       
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
}


void BTagging::GetPhysicsFlavour(Candidate* jet, TIter & itPartonArray){

  Candidate tempParticle ;
  Candidate tempNearest;
  float     minDr = 1000;
  int       nInTheCone = 0;
  float     TheBiggerConeSize = 0.7;
  Candidate* parton;
  std::vector<Candidate> contaminations;
  contaminations.clear();

  itPartonArray.Reset();
  while((parton = static_cast<Candidate*>(itPartonArray.Next()))){ 
      float dist = jet->Momentum.DeltaR(parton->Momentum); // take the DR
      if(parton->Status == 3 and dist < minDr ){
	    tempNearest = *parton; 
            minDr       = dist ;
      }

      if( parton->Status == 3  and dist <= fDeltaR){
	     tempParticle = *parton;
             nInTheCone++;
      }

      if(parton->Status != 3 and (parton->D1 != -1 or parton->D2!=-1)){         
        if((TMath::Abs(parton->PID) < 4 or TMath::Abs(parton->PID) == 21) and dist < TheBiggerConeSize)
	  contaminations.push_back(*parton);
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
      if( numberOfMothers > 0 && dynamic_cast<Candidate*>(fParticleInputArray->At(contaminations.at(iPart).M1))->Momentum == tempParticle.Momentum) continue; // mother is the initialParton --> OK
      if( TMath::Abs(tempParticle.PID) == 4 ) {
	 if( contaminatingFlavour == 4 ) continue; // keep association --> the initialParton is a c --> the contaminated parton is a c
	 jet->flavourPhysics = 0; // all the other cases reject!
      }
    }
  }
  
}
