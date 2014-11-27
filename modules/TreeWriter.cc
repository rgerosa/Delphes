/** \class TreeWriter
 *  Fills ROOT tree branches.
 *  $Date: 2013-05-16 16:28:38 +0200 (Thu, 16 May 2013) $
 *  $Revision: 1115 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 */

#include "modules/TreeWriter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

#include "TROOT.h"
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

TreeWriter::TreeWriter(){}

//------------------------------------------------------------------------------

TreeWriter::~TreeWriter(){}

//------------------------------------------------------------------------------

void TreeWriter::Init() {

  fClassMap[GenParticle::Class()] = &TreeWriter::ProcessParticles;
  fClassMap[Track::Class()]       = &TreeWriter::ProcessTracks;
  fClassMap[Tower::Class()]       = &TreeWriter::ProcessTowers;
  fClassMap[Photon::Class()]      = &TreeWriter::ProcessPhotons;
  fClassMap[Electron::Class()]    = &TreeWriter::ProcessElectrons;
  fClassMap[Muon::Class()]        = &TreeWriter::ProcessMuons;
  fClassMap[Jet::Class()]         = &TreeWriter::ProcessJets;
  fClassMap[MissingET::Class()]   = &TreeWriter::ProcessMissingET;
  fClassMap[ScalarHT::Class()]    = &TreeWriter::ProcessScalarHT;
  fClassMap[Rho::Class()]         = &TreeWriter::ProcessRho;
  fClassMap[IsoTrack::Class()]    = &TreeWriter::ProcessIsoTracks;
  fClassMap[LHEParticle::Class()] = &TreeWriter::ProcessLHEParticles;

  TBranchMap::iterator itBranchMap;
  map< TClass *, TProcessMethod >::iterator itClassMap;

  // read branch configuration and
  // import array with output from filter/classifier/jetfinder modules

  fOffsetFromModifyBeamSpot = GetInt("OffsetFromModifyBeamSpot", 0);

  ExRootConfParam param = GetParam("Branch");
  Long_t i, size;
  TString branchName, branchClassName, branchInputArray;
  TClass *branchClass;
  TObjArray *array;
  ExRootTreeBranch *branch;

  size = param.GetSize();
  for(i = 0; i < size/3; ++i){

    branchInputArray = param[i*3].GetString();
    branchName       = param[i*3 + 1].GetString();
    branchClassName  = param[i*3 + 2].GetString();

    branchClass = gROOT->GetClass(branchClassName);

    if(!branchClass){
      cout << "** ERROR: cannot find class '" << branchClassName << "'" << endl;
      continue;
    }

    itClassMap = fClassMap.find(branchClass);
    if(itClassMap == fClassMap.end()){
      cout << "** ERROR: cannot create branch for class '" << branchClassName << "'" << endl;
      continue;
    }

    array = ImportArray(branchInputArray);
    branch = NewBranch(branchName, branchClass);

    fBranchMap.insert(make_pair(branch, make_pair(itClassMap->second, array)));
  }

}

//------------------------------------------------------------------------------

void TreeWriter::Finish(){}

//------------------------------------------------------------------------------

void TreeWriter::FillParticles(Candidate *candidate, TRefArray *array){

  TIter it1(candidate->GetCandidates());
  it1.Reset();
  array->Clear();
  while((candidate = static_cast<Candidate*>(it1.Next()))){
    TIter it2(candidate->GetCandidates());
    // particle
    if(candidate->GetCandidates()->GetEntriesFast() == 0){
      array->Add(candidate);
      continue;
    }

    candidate = static_cast<Candidate*>(candidate->GetCandidates()->Last()); //?
    if(candidate->GetCandidates()->GetEntriesFast() == 0){
      array->Add(candidate);
      continue;
    }

    // tower
    it2.Reset();
    while((candidate = static_cast<Candidate*>(it2.Next()))){
      array->Add(candidate->GetCandidates()->Last());
    }
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessParticles(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  GenParticle *entry = 0;
  Double_t pt, signPz, cosTheta, eta, rapidity;

  // loop over all particles
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    entry = static_cast<GenParticle*>(branch->NewEntry()); // create the output branch
    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID()); // take the id of the candidate we are looking at 

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    // store the info in the output branch
    entry->PID    = candidate->PID;
    entry->Status = candidate->Status;
    entry->IsPU   = candidate->IsPU;
    entry->M1     = candidate->M1;
    entry->M2     = candidate->M2;
    entry->D1     = candidate->D1;
    entry->D2     = candidate->D2;
    entry->Charge = candidate->Charge;
    entry->Mass   = candidate->Mass;

    entry->E  = momentum.E();
    entry->Px = momentum.Px();
    entry->Py = momentum.Py();
    entry->Pz = momentum.Pz();

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT  = pt;
    entry->Rapidity = rapidity;

    entry->X = position.X();
    entry->Y = position.Y();
    entry->Z = position.Z();
    entry->T = position.T();

  }
}

//------------------------------------------------------------------------------
void TreeWriter::ProcessTracks(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  Candidate *particle = 0;
  Track *entry = 0;
  Double_t pt, signz, cosTheta, eta;

  // loop over all tracks
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){

    const TLorentzVector &position = candidate->Position;
   
    cosTheta = TMath::Abs(position.CosTheta());
    signz    = (position.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signz*999.9 : position.Eta());

    entry = static_cast<Track*>(branch->NewEntry()); // create the output branch
    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->PID    = candidate->PID;
    entry->Status = candidate->Status;
    entry->Charge = candidate->Charge;

    entry->EtaOuter = eta;
    entry->PhiOuter = position.Phi();

    entry->XOuter = position.X();
    entry->YOuter = position.Y();
    entry->ZOuter = position.Z();
    entry->TOuter = position.T();

    const TLorentzVector &momentum = candidate->Momentum;
    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signz    = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signz*999.9 : momentum.Eta());

    entry->Eta  = eta;
    entry->Phi  = momentum.Phi();
    entry->PT   = pt;
    entry->Mass = momentum.M();

    particle = static_cast<Candidate*>(candidate->GetCandidates()->At(fOffsetFromModifyBeamSpot));
    const TLorentzVector &initialPosition = particle->Position;

    entry->X = initialPosition.X();
    entry->Y = initialPosition.Y();
    entry->Z = initialPosition.Z();
    entry->T = initialPosition.T();

    entry->IsPU     = candidate->IsPU;
    entry->IsRecoPU = candidate->IsRecoPU;
    entry->Particle = particle; // save the reference

  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessTowers(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  Tower *entry = 0;
  Double_t pt, signPz, cosTheta, eta;

  // loop over all jets
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){

    const TLorentzVector &momentum = candidate->Momentum;

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<Tower*>(branch->NewEntry());
    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta  = eta;
    entry->Phi  = momentum.Phi();
    entry->ET   = pt;
    entry->E    = momentum.E();
    entry->Eem  = candidate->Eem;
    entry->Ehad = candidate->Ehad;
    entry->Edges[0] = candidate->Edges[0];
    entry->Edges[1] = candidate->Edges[1];
    entry->Edges[2] = candidate->Edges[2];
    entry->Edges[3] = candidate->Edges[3];

    entry->TOuter = candidate->Position.T();
    entry->nTimes = candidate->nTimes;

    FillParticles(candidate,&entry->Particles); // save the reference
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessPhotons(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  Photon *entry = 0;
  Double_t pt, signPz, cosTheta, eta;

  // loop over all photons
  array->Sort();
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){

    TIter it1(candidate->GetCandidates());
    const TLorentzVector &momentum = candidate->Momentum;

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<Photon*>(branch->NewEntry());

    
    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT  = pt;
    entry->E   = momentum.E();
    entry->EhadOverEem = candidate->Eem > 0.0 ? candidate->Ehad/candidate->Eem : 999.9;
    entry->TOuter      = candidate->Position.T();

    entry->Status   = candidate->Status;
    entry->IsRecoPU = candidate->IsRecoPU;
    entry->IsPU     = candidate->IsPU;
    entry->IsFakeObject = candidate->IsFakeObject;

    entry->IsolationVarDBeta   = candidate->IsolationVarDBeta;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->chargedHadronEnergy = candidate->chargedHadronEnergy;
    entry->neutralEnergy       = candidate->neutralEnergy;
    entry->chargedPUEnergy     = candidate->chargedPUEnergy;
    entry->allParticleEnergy   = candidate->allParticleEnergy;

    FillParticles(candidate, &entry->Particles);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessElectrons(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  Electron *entry = 0;
  Double_t pt, signPz, cosTheta, eta;

  array->Sort();
  // loop over all electrons
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    const TLorentzVector &momentum = candidate->Momentum;

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<Electron*>(branch->NewEntry());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT  = pt;

    entry->Charge = candidate->Charge;

    entry->Status   = candidate->Status;
    entry->IsRecoPU = candidate->IsRecoPU;
    entry->IsPU     = candidate->IsPU;
    entry->IsFakeObject = candidate->IsFakeObject;

    entry->IsolationVarDBeta   = candidate->IsolationVarDBeta;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->chargedHadronEnergy = candidate->chargedHadronEnergy;
    entry->neutralEnergy       = candidate->neutralEnergy;
    entry->chargedPUEnergy     = candidate->chargedPUEnergy;
    entry->allParticleEnergy   = candidate->allParticleEnergy;

    entry->EhadOverEem = 0.0;
    entry->TOuter = candidate->Position.T();

    entry->Particle = candidate->GetCandidates()->At(fOffsetFromModifyBeamSpot);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMuons(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  Muon *entry = 0;
  Double_t pt, signPz, cosTheta, eta;

  array->Sort();

  // loop over all muons
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    const TLorentzVector &momentum = candidate->Momentum;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<Muon*>(branch->NewEntry());
    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;

    entry->Charge = candidate->Charge;

    entry->Status   = candidate->Status;
    entry->IsRecoPU = candidate->IsRecoPU;
    entry->IsPU     = candidate->IsPU;
    entry->IsFakeObject = candidate->IsFakeObject;

    entry->IsolationVarDBeta   = candidate->IsolationVarDBeta;
    entry->IsolationVarRhoCorr = candidate->IsolationVarRhoCorr;
    entry->chargedHadronEnergy = candidate->chargedHadronEnergy;
    entry->neutralEnergy       = candidate->neutralEnergy;
    entry->chargedPUEnergy     = candidate->chargedPUEnergy;
    entry->allParticleEnergy   = candidate->allParticleEnergy;

    entry->Particle = candidate->GetCandidates()->At(fOffsetFromModifyBeamSpot);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessIsoTracks(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  IsoTrack *entry = 0;
  Double_t pt, signPz, cosTheta, eta;

  array->Sort();

  // loop over all IsoTracks
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    const TLorentzVector &momentum = candidate->Momentum;

    pt = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<IsoTrack*>(branch->NewEntry());

    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID());

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT = pt;
    entry->IsolationVar = candidate->TrackIsolationVar;

    entry->Charge = candidate->Charge;
    entry->IsEMCand = candidate->IsEMCand;

    entry->Particle = candidate->GetCandidates()->At(fOffsetFromModifyBeamSpot);
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessJets(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0, *constituent = 0;
  Jet *entry = 0;
  Double_t pt, signPz, cosTheta, eta;
  Double_t ecalEnergy, hcalEnergy;

  array->Sort();

  // loop over all jets
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    TIter itConstituents(candidate->GetCandidates());
    const TLorentzVector &momentum = candidate->Momentum;

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());

    entry = static_cast<Jet*>(branch->NewEntry());

    entry->Eta  = eta;
    entry->Phi  = momentum.Phi();
    entry->PT   = pt;
    entry->Mass = momentum.M();

    entry->DeltaEta = candidate->DeltaEta;
    entry->DeltaPhi = candidate->DeltaPhi;


    entry->Tau1 = candidate->Tau1;
    entry->Tau2 = candidate->Tau2;
    entry->Tau3 = candidate->Tau3;

    entry->NSubJetsTrimmed = candidate->NSubJetsTrimmed ;

    entry->TrimmedMass = candidate->TrimmedMass;
    entry->TrimmedPt   = candidate->TrimmedPt;
    entry->TrimmedEta  = candidate->TrimmedEta;
    entry->TrimmedPhi  = candidate->TrimmedPhi;

    entry->TrimmedMassSub1 = candidate->TrimmedMassSub1;
    entry->TrimmedPtSub1   = candidate->TrimmedPtSub1;
    entry->TrimmedEtaSub1  = candidate->TrimmedEtaSub1;
    entry->TrimmedPhiSub1  = candidate->TrimmedPhiSub1;

    entry->TrimmedMassSub2 = candidate->TrimmedMassSub2;
    entry->TrimmedPtSub2   = candidate->TrimmedPtSub2;
    entry->TrimmedEtaSub2  = candidate->TrimmedEtaSub2;
    entry->TrimmedPhiSub2  = candidate->TrimmedPhiSub2;

    entry->TrimmedMassSub3 = candidate->TrimmedMassSub3;
    entry->TrimmedPtSub3   = candidate->TrimmedPtSub3;
    entry->TrimmedEtaSub3  = candidate->TrimmedEtaSub3;
    entry->TrimmedPhiSub3  = candidate->TrimmedPhiSub3;

    entry->NSubJetsPruned = candidate->NSubJetsPruned ;

    entry->PrunedMass = candidate->PrunedMass;
    entry->PrunedPt   = candidate->PrunedPt;
    entry->PrunedEta  = candidate->PrunedEta;
    entry->PrunedPhi  = candidate->PrunedPhi;

    entry->PrunedMassSub1 = candidate->PrunedMassSub1;
    entry->PrunedPtSub1   = candidate->PrunedPtSub1;
    entry->PrunedEtaSub1  = candidate->PrunedEtaSub1;
    entry->PrunedPhiSub1  = candidate->PrunedPhiSub1;

    entry->PrunedMassSub2 = candidate->PrunedMassSub2;
    entry->PrunedPtSub2   = candidate->PrunedPtSub2;
    entry->PrunedEtaSub2  = candidate->PrunedEtaSub2;
    entry->PrunedPhiSub2  = candidate->PrunedPhiSub2;

    entry->PrunedMassSub3 = candidate->PrunedMassSub3;
    entry->PrunedPtSub3   = candidate->PrunedPtSub3;
    entry->PrunedEtaSub3  = candidate->PrunedEtaSub3;
    entry->PrunedPhiSub3  = candidate->PrunedPhiSub3;

    entry->NSubJetsSoftDrop = candidate->NSubJetsSoftDrop ;

    entry->SoftDropMass = candidate->SoftDropMass;
    entry->SoftDropPt   = candidate->SoftDropPt;
    entry->SoftDropEta  = candidate->SoftDropEta;
    entry->SoftDropPhi  = candidate->SoftDropPhi;

    entry->SoftDropMassSub1 = candidate->SoftDropMassSub1;
    entry->SoftDropPtSub1   = candidate->SoftDropPtSub1;
    entry->SoftDropEtaSub1  = candidate->SoftDropEtaSub1;
    entry->SoftDropPhiSub1  = candidate->SoftDropPhiSub1;

    entry->SoftDropMassSub2 = candidate->SoftDropMassSub2;
    entry->SoftDropPtSub2   = candidate->SoftDropPtSub2;
    entry->SoftDropEtaSub2  = candidate->SoftDropEtaSub2;
    entry->SoftDropPhiSub2  = candidate->SoftDropPhiSub2;

    entry->SoftDropMassSub3 = candidate->SoftDropMassSub3;
    entry->SoftDropPtSub3   = candidate->SoftDropPtSub3;
    entry->SoftDropEtaSub3  = candidate->SoftDropEtaSub3;
    entry->SoftDropPhiSub3  = candidate->SoftDropPhiSub3;

    entry->AreaX = candidate->Area.X();
    entry->AreaY = candidate->Area.Y();
    entry->AreaZ = candidate->Area.Z();
    entry->AreaT = candidate->Area.T();

    entry->Charge = candidate->Charge;

    entry->TauTag = candidate->TauTag;

    entry->BTagAlgo      = candidate->BTagAlgo;
    entry->BTagDefault      = candidate->BTagDefault;
    entry->BTagPhysics   = candidate->BTagPhysics;
    entry->BTagNearest2  = candidate->BTagNearest2;
    entry->BTagNearest3  = candidate->BTagNearest3;
    entry->BTagHeaviest  = candidate->BTagHeaviest;
    entry->BTagHighestPt = candidate->BTagHighestPt;

    entry->flavourAlgo      = candidate->flavourAlgo;
    entry->flavourDefault      = candidate->flavourDefault;
    entry->flavourPhysics   = candidate->flavourPhysics;
    entry->flavourNearest2  = candidate->flavourNearest2;
    entry->flavourNearest3  = candidate->flavourNearest3;
    entry->flavourHeaviest  = candidate->flavourHeaviest;
    entry->flavourHighestPt = candidate->flavourHighestPt;

    itConstituents.Reset();
    entry->Constituents.Clear();
    ecalEnergy = 0.0;
    hcalEnergy = 0.0;
    while((constituent = static_cast<Candidate*>(itConstituents.Next()))){
      entry->Constituents.Add(constituent);
      ecalEnergy += constituent->Eem;
      hcalEnergy += constituent->Ehad;
    }

    entry->EhadOverEem = ecalEnergy > 0.0 ? hcalEnergy/ecalEnergy : 999.9;

    // pileup jet ID                                                                                                                                                                          
    entry->dRMean   = candidate->dRMean;
    entry->dR2Mean  = candidate->dR2Mean;
    entry->ptD      = candidate->ptD;
    entry->sumPt    = candidate->sumPt;
    entry->sumPt2   = candidate->sumPt2;

    entry->dRMeanEm  = candidate->dRMeanEm;
    entry->ptDNe     = candidate->ptDNe;
    entry->sumPtNe   = candidate->sumPtNe;
    entry->nNeutral  = candidate->nNeutral;
    entry->neuEMfrac   = candidate->neuEMfrac;
    entry->dRMeanNeut  = candidate->dRMeanNeut;
    entry->neuHadfrac  = candidate->neuHadfrac;

    entry->dRMeanCh  = candidate->dRMeanCh;
    entry->ptDCh     = candidate->ptDCh;
    entry->sumPtCh   = candidate->sumPtCh;
    entry->nCharged  = candidate->nCharged;

    entry->chgEMfrac   = candidate->chgEMfrac;
    entry->chgHadfrac  = candidate->chgHadfrac;

    entry->betaClassic     = candidate->betaClassic;
    entry->betaClassicStar = candidate->betaClassicStar;
    entry->beta          = candidate->beta;
    entry->betaStar      = candidate->betaStar;
    entry->constituents  = candidate->constituents;

    entry->dZ  = candidate->dZ;
    entry->d0  = candidate->d0;

    entry->etaW  = candidate->etaW;
    entry->phiW  = candidate->phiW;
    entry->jetW  = candidate->jetW;

    entry->majW  = candidate->majW;
    entry->minW  = candidate->minW;
    entry->dRLeadCent  = candidate->dRLeadCent;
    entry->dRLead2nd   = candidate->dRLead2nd;

    entry->ptMean  = candidate->ptMean;
    entry->ptRMS   = candidate->ptRMS;
    entry->pt2A    = candidate->pt2A;
    entry->sumChPt  = candidate->sumChPt;
    entry->sumNePt  = candidate->sumNePt;
    entry->axis2    = candidate->axis2;

    entry->leadFrac   = candidate->leadFrac;
    entry->secondFrac = candidate->secondFrac;
    entry->thirdFrac  = candidate->thirdFrac;
    entry->fourthFrac = candidate->fourthFrac;

    entry->leadChFrac   = candidate->leadChFrac;
    entry->secondChFrac = candidate->secondChFrac;
    entry->thirdChFrac  = candidate->thirdChFrac;
    entry->fourthChFrac = candidate->fourthChFrac;

    entry->leadEmFrac   = candidate->leadEmFrac;
    entry->secondEmFrac = candidate->secondEmFrac;
    entry->thirdEmFrac  = candidate->thirdEmFrac;
    entry->fourthEmFrac = candidate->fourthEmFrac;

    entry->leadNeutFrac   = candidate->leadNeutFrac;
    entry->secondNeutFrac = candidate->secondNeutFrac;
    entry->thirdNeutFrac  = candidate->thirdNeutFrac;
    entry->fourthNeutFrac = candidate->fourthNeutFrac;

    entry->pileupIDFlagCutBased = candidate->pileupIDFlagCutBased;

    entry->FracPt     = candidate->FracPt;
    entry->emFracPt   = candidate->emFracPt;
    entry->neutFracPt = candidate->neutFracPt;
    entry->chFracPt   = candidate->chFracPt;

    FillParticles(candidate, &entry->Particles);

        

  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessMissingET(ExRootTreeBranch *branch, TObjArray *array){

  Candidate *candidate = 0;
  MissingET *entry = 0;

  // get the first entry
  if((candidate = static_cast<Candidate*>(array->At(0)))){
    const TLorentzVector &momentum = candidate->Momentum;
    entry = static_cast<MissingET*>(branch->NewEntry());
    entry->Phi = (-momentum).Phi();
    entry->MET = momentum.Pt();
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessScalarHT(ExRootTreeBranch *branch, TObjArray *array){
  Candidate *candidate = 0;
  ScalarHT *entry = 0;

  // get the first entry
  if((candidate = static_cast<Candidate*>(array->At(0)))){
    const TLorentzVector &momentum = candidate->Momentum;
    entry = static_cast<ScalarHT*>(branch->NewEntry());
    entry->HT = momentum.Pt();
  }
}

//------------------------------------------------------------------------------

void TreeWriter::ProcessRho(ExRootTreeBranch *branch, TObjArray *array){
  TIter iterator(array);
  Candidate *candidate = 0;
  Rho *entry = 0;

  // loop over all rho
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
      const TLorentzVector &momentum = candidate->Momentum;
      entry = static_cast<Rho*>(branch->NewEntry());
      entry->Rho = momentum.E();
      entry->Edges[0] = candidate->Edges[0];
      entry->Edges[1] = candidate->Edges[1];
    }
}

// -----------------------------------------------------------------------------
//------------------------------------------------------------------------------

void TreeWriter::ProcessLHEParticles(ExRootTreeBranch *branch, TObjArray *array){

  TIter iterator(array);
  Candidate *candidate = 0;
  LHEParticle *entry = 0;
  Double_t pt, signPz, cosTheta, eta, rapidity;

  // loop over all particles
  iterator.Reset();
  while((candidate = static_cast<Candidate*>(iterator.Next()))){
    const TLorentzVector &momentum = candidate->Momentum;
    const TLorentzVector &position = candidate->Position;

    entry = static_cast<LHEParticle*>(branch->NewEntry()); // create the output branch
    entry->SetBit(kIsReferenced);
    entry->SetUniqueID(candidate->GetUniqueID()); // take the id of the candidate we are looking at 

    pt       = momentum.Pt();
    cosTheta = TMath::Abs(momentum.CosTheta());
    signPz   = (momentum.Pz() >= 0.0) ? 1.0 : -1.0;
    eta      = (cosTheta == 1.0 ? signPz*999.9 : momentum.Eta());
    rapidity = (cosTheta == 1.0 ? signPz*999.9 : momentum.Rapidity());

    // store the info in the output branch
    entry->PID    = candidate->PID; // pdgId 
    entry->Status = candidate->Status; // status written in the LHE file
    entry->IsPU   = candidate->IsPU;
    entry->M1     = candidate->M1;
    entry->M2     = candidate->M2;
    entry->D1     = candidate->D1;
    entry->D2     = candidate->D2;
    entry->Charge = candidate->Charge;
    entry->Mass   = candidate->Mass;

    entry->E  = momentum.E();
    entry->Px = momentum.Px();
    entry->Py = momentum.Py();
    entry->Pz = momentum.Pz();

    entry->Eta = eta;
    entry->Phi = momentum.Phi();
    entry->PT  = pt;
    entry->Rapidity = rapidity;

    entry->X = position.X();
    entry->Y = position.Y();
    entry->Z = position.Z();
    entry->T = position.T();
  }
}


//------------------------------------------------------------------------------

void TreeWriter::Process(){
  TBranchMap::iterator itBranchMap;
  ExRootTreeBranch *branch;
  TProcessMethod method;
  TObjArray *array;

  for(itBranchMap = fBranchMap.begin(); itBranchMap != fBranchMap.end(); ++itBranchMap){
    branch = itBranchMap->first;
    method = itBranchMap->second.first;
    array = itBranchMap->second.second;

    (this->*method)(branch, array);
  }
}

//------------------------------------------------------------------------------
