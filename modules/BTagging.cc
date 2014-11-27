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

BTagging::BTagging() :
  fItJetInputArray(0){
}

//------------------------------------------------------------------------------

BTagging::~BTagging(){}

//------------------------------------------------------------------------------

void BTagging::Init(){

  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  ExRootConfParam param;
  DelphesFormula *formula;
  Int_t i, size;

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaLoose");
  size  = param.GetSize();
  
  fEfficiencyMapLoose.clear();
  for(i = 0; i < size/2; ++i){
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapLoose[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap    = fEfficiencyMapLoose.find(0);
  if(itEfficiencyMap == fEfficiencyMapLoose.end()){
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapLoose[0] = formula;
  }

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaMedium");
  size  = param.GetSize();
  
  fEfficiencyMapMedium.clear();
  for(i = 0; i < size/2; ++i){
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapMedium[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap    = fEfficiencyMapMedium.find(0);
  if(itEfficiencyMap == fEfficiencyMapMedium.end()){
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapMedium[0] = formula;
  }

  // read efficiency formulas
  param = GetParam("EfficiencyFormulaTight");
  size  = param.GetSize();
  
  fEfficiencyMapTight.clear();
  for(i = 0; i < size/2; ++i){
    formula = new DelphesFormula;
    formula->Compile(param[i*2 + 1].GetString());
    fEfficiencyMapTight[param[i*2].GetInt()] = formula;
  }

  // set default efficiency formula
  itEfficiencyMap    = fEfficiencyMapTight.find(0);
  if(itEfficiencyMap == fEfficiencyMapTight.end()){
    formula = new DelphesFormula;
    formula->Compile("0.0");
    fEfficiencyMapTight[0] = formula;
  }


  // import input array(s)
  fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

}

//------------------------------------------------------------------------------
void BTagging::Finish(){
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;
  if(fItJetInputArray) delete fItJetInputArray;

  for(itEfficiencyMap = fEfficiencyMapLoose.begin(); itEfficiencyMap != fEfficiencyMapLoose.end(); ++itEfficiencyMap){
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }

  for(itEfficiencyMap = fEfficiencyMapMedium.begin(); itEfficiencyMap != fEfficiencyMapMedium.end(); ++itEfficiencyMap){
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }

  for(itEfficiencyMap = fEfficiencyMapTight.begin(); itEfficiencyMap != fEfficiencyMapTight.end(); ++itEfficiencyMap){
    formula = itEfficiencyMap->second;
    if(formula) delete formula;
  }

}

//------------------------------------------------------------------------------

void BTagging::Process(){

  Candidate *jet;
  map< Int_t, DelphesFormula * >::iterator itEfficiencyMap;
  DelphesFormula *formula;

  
  // loop over all input jets
  fItJetInputArray->Reset();

  while((jet = static_cast<Candidate*>(fItJetInputArray->Next()))){
    const TLorentzVector &jetMomentum = jet->Momentum;
     
    float randomNumber = gRandom->Uniform();

    // --------------------------
    // Loose BTagging
    // --------------------------

    // find haviest flavour and Btag 
    if(fEfficiencyMapLoose.size() != 0) {

     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourHeaviest);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHeaviest = 1 ;
   
     // find hisghestpt flavour and Btag
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourHighestPt);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHighestPt = 1 ;

     // find nearest2 flavour and Btag
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourNearest2);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest2 = 1 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourNearest3);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest3 = 1 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourAlgo);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagAlgo = 1 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourPhysics);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagPhysics = 1 ;

     // default 
     itEfficiencyMap = fEfficiencyMapLoose.find(jet->flavourDefault);
     if(itEfficiencyMap == fEfficiencyMapLoose.end()){
      itEfficiencyMap = fEfficiencyMapLoose.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagDefault = 1 ;

    }

    // --------------------------
    // Medium BTagging
    // --------------------------

    // find haviest flavour and Btag 
    if(fEfficiencyMapMedium.size() != 0){

     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourHeaviest);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHeaviest = 2 ;
   
     // find hisghestpt flavour and Btag
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourHighestPt);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHighestPt = 2 ;

     // find nearest2 flavour and Btag
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourNearest2);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest2 = 2 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourNearest3);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest3 = 2 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourAlgo);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagAlgo = 2 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourPhysics);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
      itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagPhysics = 2 ;

     // default 
     itEfficiencyMap = fEfficiencyMapMedium.find(jet->flavourDefault);
     if(itEfficiencyMap == fEfficiencyMapMedium.end()){
       itEfficiencyMap = fEfficiencyMapMedium.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagDefault = 2 ;

    }

    // --------------------------
    // Tight BTagging
    // --------------------------

    // find haviest flavour and Btag 
    if(fEfficiencyMapTight.size() != 0){

     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourHeaviest);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHeaviest = 3 ;
   
     // find hisghestpt flavour and Btag
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourHighestPt);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagHighestPt = 3 ;

     // find nearest2 flavour and Btag
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourNearest2);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest2 = 3 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourNearest3);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagNearest3 = 3 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourAlgo);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagAlgo = 3 ;

     // find nearest3 flavour and Btag
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourPhysics);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagPhysics = 3 ;

     // default 
     itEfficiencyMap = fEfficiencyMapTight.find(jet->flavourDefault);
     if(itEfficiencyMap == fEfficiencyMapTight.end()){
      itEfficiencyMap = fEfficiencyMapTight.find(0);
     }
     formula = itEfficiencyMap->second;

     if(randomNumber <= formula->Eval(jetMomentum.Pt(),jetMomentum.Eta())) jet->BTagDefault = 3 ;
    }
  }
}
