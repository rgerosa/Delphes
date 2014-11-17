/** \class Isolation
 *  Sums transverse momenta of isolation objects (tracks, calorimeter towers, etc)
 *  within a DeltaR cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *  $Date: 2013-11-04 13:14:33 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1317 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 */

#include "modules/Isolation.h"

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

class IsolationClassifier : public ExRootClassifier{

public:
  IsolationClassifier() {}
  Int_t GetCategory(TObject *object);
  Double_t fPTMin;
};

//------------------------------------------------------------------------------

Int_t IsolationClassifier::GetCategory(TObject *object){

  Candidate *track = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = track->Momentum;
  if(momentum.Pt() < fPTMin) return -1; // pt threshold on the candidate
  return 0;
}

//------------------------------------------------------------------------------

Isolation::Isolation() :
  fClassifier(0), fNeutralFilter(0), fChargedFilter(0),
  fItNeutralIsolationInputArray(0), fItChargedIsolationInputArray(0), fItCandidateInputArray(0),
  fItRhoInputArray(0){
  fClassifier = new IsolationClassifier;
}

//------------------------------------------------------------------------------

Isolation::~Isolation(){}

//------------------------------------------------------------------------------

void Isolation::Init(){

  const char *rhoInputArrayName;

  fDeltaRMax  = GetDouble("DeltaRMax", 0.5);  // cone dimension
  fPTRatioMax = GetDouble("PTRatioMax", 0.1); // in order to cut on isolation
  fPTSumMax = GetDouble("PTSumMax", 5.0);     // maximum sum pt for absolute cut
  fUsePTSum = GetBool("UsePTSum", false);     // relative or absolute isolation cut
  fClassifier->fPTMin = GetDouble("PTMin", 0.5); // min pt (similar to PF one)

  fNeutralIsolationInputArray   = ImportArray(GetString("NeutralIsolationInputArray", "Calorimeter/eflowTowers"));
  fItNeutralIsolationInputArray = fNeutralIsolationInputArray->MakeIterator();
  fNeutralFilter                = new ExRootFilter(fNeutralIsolationInputArray); 

  fChargedIsolationInputArray   = ImportArray(GetString("ChargedIsolationInputArray", "TrackPileUpSubtractor/eflowTracks"));
  fItChargedIsolationInputArray = fChargedIsolationInputArray->MakeIterator();
  fChargedFilter                = new ExRootFilter(fChargedIsolationInputArray);

  fCandidateInputArray   = ImportArray(GetString("CandidateInputArray", "Calorimeter/electrons"));
  fItCandidateInputArray = fCandidateInputArray->MakeIterator();

  rhoInputArrayName = GetString("RhoInputArray", "");
  if(rhoInputArrayName[0] != '\0'){
    fRhoInputArray = ImportArray(rhoInputArrayName);
    fItRhoInputArray = fRhoInputArray->MakeIterator();
  }
  else{
    fRhoInputArray = 0;
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void Isolation::Finish(){
  if(fItRhoInputArray) delete fItRhoInputArray;
  if(fNeutralFilter)   delete fNeutralFilter;
  if(fChargedFilter)   delete fChargedFilter;
  if(fItCandidateInputArray)        delete fItCandidateInputArray;
  if(fItChargedIsolationInputArray) delete fItChargedIsolationInputArray;
  if(fItNeutralIsolationInputArray) delete fItNeutralIsolationInputArray;
}

//------------------------------------------------------------------------------

void Isolation::Process(){
  Candidate *candidate, *isolation, *object;
  TObjArray *chargedIsolationArray;
  TObjArray *neutralIsolationArray;
  Double_t   sumCharged, sumNeutral, sum, ratio;
  Int_t      counter;
  Double_t   eta = 0.0;
  Double_t   rho = 0.0;

  if(fRhoInputArray && fRhoInputArray->GetEntriesFast() > 0){ // if the rho array is found, save the value
    candidate = static_cast<Candidate*>(fRhoInputArray->At(0));
    rho = candidate->Momentum.Pt();
  }

  // select charged isolation objects
  fChargedFilter->Reset();
  chargedIsolationArray = fChargedFilter->GetSubArray(fClassifier,0); // get the charged particle used a filter on the PT

  // select charged isolation objects
  fNeutralFilter->Reset();
  neutralIsolationArray = fNeutralFilter->GetSubArray(fClassifier, 0); // get the neutral particle used a filter on the PT

  if(chargedIsolationArray == 0 || neutralIsolationArray == 0) return;

  TIter itChargedIsolationArray(chargedIsolationArray);
  TIter itNeutralIsolationArray(neutralIsolationArray);

  // loop over all the input lepton or photon collection for which we want to calculate the isolation
  fItCandidateInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItCandidateInputArray->Next()))){
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = TMath::Abs(candidateMomentum.Eta());
    // loop over all input tracks
    sumCharged = 0.0;
    counter    = 0;
    itChargedIsolationArray.Reset();
    while((isolation = static_cast<Candidate*>(itChargedIsolationArray.Next()))){
      const TLorentzVector &isolationMomentum = isolation->Momentum;
      if(candidateMomentum.DeltaR(isolationMomentum) <= fDeltaRMax && !candidate->Overlaps(isolation)){
        sumCharged += isolationMomentum.Pt();
        ++counter;
      }
    }
  

    // loop over all input towers
    sumNeutral = 0.0;
    counter = 0;
    itNeutralIsolationArray.Reset();
    while((isolation = static_cast<Candidate*>(itNeutralIsolationArray.Next()))){
      const TLorentzVector &isolationMomentum = isolation->Momentum;
      if(candidateMomentum.DeltaR(isolationMomentum) <= fDeltaRMax && !candidate->Overlaps(isolation)){
        sumNeutral += isolationMomentum.Pt();
        ++counter;
      }
    }

    // find rho
    rho = 0.0;
    if(fRhoInputArray){
     fItRhoInputArray->Reset();
     while((object = static_cast<Candidate*>(fItRhoInputArray->Next()))){
        if(eta >= object->Edges[0] && eta < object->Edges[1]){
          rho = object->Momentum.Pt();
        }
     }
    }

    // correct sum for pile-up contamination
    sum = sumCharged + TMath::Max(sumNeutral - TMath::Max(rho,0.0)*fDeltaRMax*fDeltaRMax*TMath::Pi(),0.0);

    ratio = sum/candidateMomentum.Pt();

    candidate->IsolationVar = ratio;

    if((fUsePTSum && sum > fPTSumMax) || ratio > fPTRatioMax) continue;
    fOutputArray->Add(candidate);
  }
}


//------------------------------------------------------------------------------
