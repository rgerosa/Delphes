/** \class Calorimeter
 *  Fills calorimeter towers, performs calorimeter resolution smearing,
 *  preselects towers hit by photons and creates energy flow objects.
 *  $Date: 2013-09-04 17:20:22 +0200 (Wed, 04 Sep 2013) $
 *  $Revision: 1280 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
*/

#include "modules/Calorimeter.h"

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

Calorimeter::Calorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0) {

  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();
}

//------------------------------------------------------------------------------

Calorimeter::~Calorimeter(){

  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;

  if(fTowerTrackArray) delete fTowerTrackArray;
  if(fItTowerTrackArray) delete fItTowerTrackArray;
}

//------------------------------------------------------------------------------

void Calorimeter::Init(){

  ExRootConfParam param, paramEtaBins, paramPhiBins, paramFractions;
  Long_t i, j, k, size, sizeEtaBins, sizePhiBins;
  Double_t ecalFraction, hcalFraction;
  TBinMap::iterator itEtaBin;
  set< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  // read eta and phi bins
  param = GetParam("EtaPhiBins");
  size = param.GetSize();
  fBinMap.clear();
  fEtaBins.clear();
  fPhiBins.clear();
  for(i = 0; i < size/2; ++i){

    paramEtaBins = param[i*2];
    sizeEtaBins = paramEtaBins.GetSize();
    paramPhiBins = param[i*2 + 1];
    sizePhiBins = paramPhiBins.GetSize();

    for(j = 0; j < sizeEtaBins; ++j){
      for(k = 0; k < sizePhiBins; ++k){
        fBinMap[paramEtaBins[j].GetDouble()].insert(paramPhiBins[k].GetDouble());
      }
    }
  }

  // for better performance we transform map of sets to parallel vectors:
  // vector< double > and vector< vector< double >* >
  for(itEtaBin = fBinMap.begin(); itEtaBin != fBinMap.end(); ++itEtaBin){
    fEtaBins.push_back(itEtaBin->first);
    phiBins = new vector< double >(itEtaBin->second.size());
    fPhiBins.push_back(phiBins);
    phiBins->clear();
    for(itPhiBin = itEtaBin->second.begin(); itPhiBin != itEtaBin->second.end(); ++itPhiBin){
      phiBins->push_back(*itPhiBin);
    }
  }

  // read energy fractions for different particles
  param = GetParam("EnergyFraction");
  size = param.GetSize();

  // set default energy fractions values
  fFractionMap.clear();
  fFractionMap[0] = make_pair(0.0, 1.0);

  for(i = 0; i < size/2; ++i){
    paramFractions = param[i*2 + 1];
    ecalFraction = paramFractions[0].GetDouble();
    hcalFraction = paramFractions[1].GetDouble();

    fFractionMap[param[i*2].GetInt()] = make_pair(ecalFraction, hcalFraction);
  }

  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  // import array with output from other modules
  fParticleInputArray = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  fTowerOutputArray = ExportArray(GetString("TowerOutputArray", "towers"));
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));

  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "eflowTowers"));

  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  electronsFromTrack = false;

  fTimingEMin = GetDouble("TimingEMin",4.);

}

//------------------------------------------------------------------------------

void Calorimeter::Finish(){

  vector< vector< Double_t >* >::iterator itPhiBin;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  for(itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin){
    delete *itPhiBin;
  }
}

//------------------------------------------------------------------------------

void Calorimeter::Process(){

  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy;
  Int_t pdgCode;

  TFractionMap::iterator itFractionMap;

  vector< Double_t >::iterator itEtaBin;
  vector< Double_t >::iterator itPhiBin;
  vector< Double_t > *phiBins;

  vector< Long64_t >::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerECalFractions.clear();
  fTowerHCalFractions.clear();
  fTrackECalFractions.clear();
  fTrackHCalFractions.clear();


  // loop over all particles of the event
  fItParticleInputArray->Reset();
  number = -1;

  while((particle = static_cast<Candidate*>(fItParticleInputArray->Next()))){

    const TLorentzVector &particlePosition = particle->Position;
    ++number;
    pdgCode = TMath::Abs(particle->PID);
    itFractionMap = fFractionMap.find(pdgCode); // find the particle in the fraction map

    if(itFractionMap == fFractionMap.end()){
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first; // take ecal fraction
    hcalFraction = itFractionMap->second.second; // take Hcal fraction
  
    // fill tower fraction vectors
    fTowerECalFractions.push_back(ecalFraction);
    fTowerHCalFractions.push_back(hcalFraction);

    if(ecalFraction < 1.0E-9 && hcalFraction < 1.0E-9) continue;

    // find eta bin [1, fEtaBins.size - 1] --> find seed position
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), particlePosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), particlePosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags  = 0;
    flags |= (pdgCode == 11 || pdgCode == 22) << 1; // flag for electrons and photons

    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for particle number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);
    fTowerHits.push_back(towerHit); // fill tower hit vectors
  }

  // loop over all tracks
  fItTrackInputArray->Reset();
  number = -1;
  while((track = static_cast<Candidate*>(fItTrackInputArray->Next()))){

    const TLorentzVector &trackPosition = track->Position;
    ++number;

    pdgCode = TMath::Abs(track->PID);

    itFractionMap = fFractionMap.find(pdgCode);
    if(itFractionMap == fFractionMap.end()){
      itFractionMap = fFractionMap.find(0);
    }

    ecalFraction = itFractionMap->second.first;
    hcalFraction = itFractionMap->second.second;

    fTrackECalFractions.push_back(ecalFraction);
    fTrackHCalFractions.push_back(hcalFraction);

    // find eta bin [1, fEtaBins.size - 1]
    itEtaBin = lower_bound(fEtaBins.begin(), fEtaBins.end(), trackPosition.Eta());
    if(itEtaBin == fEtaBins.begin() || itEtaBin == fEtaBins.end()) continue;
    etaBin = distance(fEtaBins.begin(), itEtaBin);

    // phi bins for given eta bin
    phiBins = fPhiBins[etaBin];

    // find phi bin [1, phiBins.size - 1]
    itPhiBin = lower_bound(phiBins->begin(), phiBins->end(), trackPosition.Phi());
    if(itPhiBin == phiBins->begin() || itPhiBin == phiBins->end()) continue;
    phiBin = distance(phiBins->begin(), itPhiBin);

    flags = 1; // charged tracks flags equal one add to the tower hits
 
    // make tower hit {16-bits for eta bin number, 16-bits for phi bin number, 8-bits for flags, 24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  }

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(), fTowerHits.end()); // sort the hit

  // loop over all hits
  towerEtaPhi = 0;
  fTower = 0;
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits) { // loop over the hits
    towerHit  = (*itTowerHits);
    flags     = (towerHit >> 24) & 0x00000000000000FFLL; // take the flag
    number    = (towerHit) & 0x0000000000FFFFFFLL;       // take the number
    hitEtaPhi = towerHit >> 32;
    if(towerEtaPhi != hitEtaPhi){  // first time hit
      // switch to next tower
      towerEtaPhi = hitEtaPhi;
      // finalize previous tower
      FinalizeTower();
      // create new tower
      fTower = factory->NewCandidate();
      // store which type of particle it belongs to
      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;
      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fTowerECalEnergy = 0.0;
      fTowerHCalEnergy = 0.0;

      fTrackECalEnergy = 0.0;
      fTrackHCalEnergy = 0.0;

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTrackArray->Clear();
    }

    // check for track hits
    if(flags & 1){
      ++fTowerTrackHits;
      track = static_cast<Candidate*>(fTrackInputArray->At(number));
      momentum = track->Momentum;

      bool dbg_scz = false;
      if (dbg_scz) {
	cout << "   Calorimeter input track has x y z t " << track->Position.X() << " " << track->Position.Y() << " " << track->Position.Z() << " " << track->Position.T() 
	     << endl;
	Candidate *prt = static_cast<Candidate*>(track->GetCandidates()->Last());
	const TLorentzVector &ini = prt->Position;

	cout << "                and parent has x y z t " << ini.X() << " " << ini.Y() << " " << ini.Z() << " " << ini.T();

      }

      ecalEnergy = momentum.E() * fTrackECalFractions[number];
      hcalEnergy = momentum.E() * fTrackHCalFractions[number];

      fTrackECalEnergy += ecalEnergy;
      fTrackHCalEnergy += hcalEnergy;

      if (ecalEnergy > fTimingEMin && fTower) {
	if (electronsFromTrack) {
	  //	  cout << " SCZ Debug pushing back track hit E=" << ecalEnergy << " T=" << track->Position.T() << " isPU=" << track->IsPU << " isRecoPU=" << track->IsRecoPU 
	  //	       << " PID=" << track->PID << endl;
	  fTower->ecal_E_t.push_back(std::make_pair<float,float>(ecalEnergy,track->Position.T()));
	} else {
	  //	  cout << " Skipping track hit E=" << ecalEnergy << " T=" << track->Position.T() << " isPU=" << track->IsPU << " isRecoPU=" << track->IsRecoPU 
	  //	       << " PID=" << track->PID << endl;
	}
      }

      fTowerTrackArray->Add(track);

      continue;
    }

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;

    // fill current tower
    ecalEnergy = momentum.E() * fTowerECalFractions[number];
    hcalEnergy = momentum.E() * fTowerHCalFractions[number];

    fTowerECalEnergy += ecalEnergy;
    fTowerHCalEnergy += hcalEnergy;

    if (ecalEnergy > fTimingEMin && fTower) {
      if (abs(particle->PID) != 11 || !electronsFromTrack) {
	//	cout << " SCZ Debug About to push back particle hit E=" << ecalEnergy << " T=" << particle->Position.T() << " isPU=" << particle->IsPU 
	//   << " PID=" << particle->PID << endl;
	fTower->ecal_E_t.push_back(std::make_pair<float,float>(ecalEnergy,particle->Position.T()));
      } else {
	
	// N.B. Only charged particle set to leave ecal energy is the electrons
	//	cout << " SCZ Debug To avoid double-counting, skipping particle hit E=" << ecalEnergy << " T=" << particle->Position.T() << " isPU=" << particle->IsPU 
	//	     << " PID=" << particle->PID << endl;
	
      }
    }

    fTower->AddCandidate(particle);
  }
  
  // finalize last tower
  FinalizeTower();
}

//------------------------------------------------------------------------------

void Calorimeter::FinalizeTower(){

  Candidate *track, *tower;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;

  if(!fTower) return;

  ecalSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy); // take ecal resolution
  ecalEnergy = LogNormal(fTowerECalEnergy, ecalSigma);

  hcalSigma  = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy);
  hcalEnergy = LogNormal(fTowerHCalEnergy, hcalSigma);

  energy     = ecalEnergy + hcalEnergy;

  eta = fTowerEta;
  phi = fTowerPhi;

  pt = energy / TMath::CosH(eta);

  // Time calculation for tower
  fTower->nTimes = 0;
  float tow_sumT = 0;
  float tow_sumW = 0;

  for (unsigned int i = 0 ; i < fTower->ecal_E_t.size() ; i++) { 
    float w = TMath::Sqrt(fTower->ecal_E_t[i].first);
    tow_sumT += w*fTower->ecal_E_t[i].second;
    tow_sumW += w;
    fTower->nTimes++;
  }
  
  if (tow_sumW > 0.) {
    fTower->Position.SetPtEtaPhiE(1.0, eta, phi,tow_sumT/tow_sumW);
  } else {
    fTower->Position.SetPtEtaPhiE(1.0,eta,phi,999999.);
  }

  //  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, 0.);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
  fTower->Ehad = hcalEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];

  // fill calorimeter towers and photon candidates
  if(energy > 0.0){
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0){
      fPhotonOutputArray->Add(fTower);
    }
    fTowerOutputArray->Add(fTower);
  }

  // save all the tracks as energy flow tracks
  fItTowerTrackArray->Reset();
  while((track = static_cast<Candidate*>(fItTowerTrackArray->Next()))){
    fEFlowTrackOutputArray->Add(track);
  }

  ecalEnergy -= fTrackECalEnergy;
  if(ecalEnergy < 0.0) ecalEnergy = 0.0;

  hcalEnergy -= fTrackHCalEnergy;
  if(hcalEnergy < 0.0) hcalEnergy = 0.0;

  energy = ecalEnergy + hcalEnergy;

  // save ECAL and/or HCAL energy excess as an energy flow tower
  if(energy > 0.0){
    // create new tower
    tower = static_cast<Candidate*>(fTower->Clone());
    pt = energy / TMath::CosH(eta);
    tower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    tower->Eem = ecalEnergy;
    tower->Ehad = hcalEnergy;
    fEFlowTowerOutputArray->Add(tower);
  }
}

//------------------------------------------------------------------------------

Double_t Calorimeter::LogNormal(Double_t mean, Double_t sigma)
{
  Double_t a, b;

  if(mean > 0.0)
  {
    b = TMath::Sqrt(TMath::Log((1.0 + (sigma*sigma)/(mean*mean))));
    a = TMath::Log(mean) - 0.5*b*b;

    return TMath::Exp(a + b*gRandom->Gaus(0, 1));
  }
  else
  {
    return 0.0;
  }
}
