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
  fECalEfficiencyFormula(0), fHCalEfficiencyFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0) {

  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;
  fECalEfficiencyFormula = new DelphesFormula;
  fHCalEfficiencyFormula = new DelphesFormula;

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();

}

//------------------------------------------------------------------------------

Calorimeter::~Calorimeter(){

  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;
  if(fECalEfficiencyFormula) delete fECalEfficiencyFormula;
  if(fHCalEfficiencyFormula) delete fHCalEfficiencyFormula;
  if(fTowerTrackArray)   delete fTowerTrackArray;
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

  // some basic info
  fECalEnergyMin = GetDouble("ECalEnergyMin", 0.0);
  fHCalEnergyMin = GetDouble("HCalEnergyMin", 0.0);

  fECalEnergySignificanceMin = GetDouble("ECalEnergySignificanceMin", 0.0);
  fHCalEnergySignificanceMin = GetDouble("HCalEnergySignificanceMin", 0.0);

  // switch on or off the dithering of the center of calorimeter towers
  fSmearTowerCenter = GetBool("SmearTowerCenter", true);

  // read resolution formulas
  fECalResolutionFormula->Compile(GetString("ECalResolutionFormula", "0"));
  fHCalResolutionFormula->Compile(GetString("HCalResolutionFormula", "0"));

  fECalEfficiencyFormula->Compile(GetString("ECalEfficiencyFormula", "0"));
  fHCalEfficiencyFormula->Compile(GetString("HCalEfficiencyFormula", "0"));

  // import array with output from other modules
  fParticleInputArray   = ImportArray(GetString("ParticleInputArray", "ParticlePropagator/particles"));
  fItParticleInputArray = fParticleInputArray->MakeIterator();

  fTrackInputArray   = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fItTrackInputArray = fTrackInputArray->MakeIterator();

  // create output arrays
  // calo tower output
  fTowerOutputArray  = ExportArray(GetString("TowerOutputArray", "towers"));
  // photons 
  fPhotonOutputArray = ExportArray(GetString("PhotonOutputArray", "photons"));
  // track for particle flow
  fEFlowTrackOutputArray = ExportArray(GetString("EFlowTrackOutputArray", "eflowTracks"));
  // tower for particle flow
  fEFlowTowerOutputArray = ExportArray(GetString("EFlowTowerOutputArray", "eflowTowers"));

  // For timing
  // So far this flag needs to be false
  // Curved extrapolation not supported
  electronsFromTrack = false;

}

//------------------------------------------------------------------------------

void Calorimeter::Finish(){

  vector< vector< Double_t >* >::iterator itPhiBin;
  if(fItParticleInputArray) delete fItParticleInputArray;
  if(fItTrackInputArray)    delete fItTrackInputArray;
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

  vector<Long64_t>::iterator itTowerHits;

  DelphesFactory *factory = GetFactory();
  fTowerHits.clear();
  fTowerECalFractions.clear();
  fTowerHCalFractions.clear();
  fTrackECalFractions.clear();
  fTrackHCalFractions.clear();
  fParticlePDGId.clear();
  fTrackPDGId.clear();


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
    fTrackPDGId.push_back(pdgCode);
 
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
    fParticlePDGId.push_back(pdgCode);

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


  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(),fTowerHits.end()); // sort the hit

  // loop over all hits
  towerEtaPhi = 0;
  fTower = 0;
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits) { // loop over the hits
    towerHit  = (*itTowerHits);
    flags     = (towerHit >> 24) & 0x00000000000000FFLL; // take the flag
    number    = (towerHit) & 0x0000000000FFFFFFLL;       // take the number
    hitEtaPhi = towerHit >> 32;

    if(towerEtaPhi != hitEtaPhi){  // first time hit, no  more than once since we have sorted tower as a function of eta and phi hits

      // switch to next tower
      towerEtaPhi = hitEtaPhi;

      // finalize previous tower
      FinalizeTower();

      // create new tower using the calorimeter information
      fTower = factory->NewCandidate();

      // store which type of particle it belongs to
      phiBin = (towerHit >> 32) & 0x000000000000FFFFLL;
      etaBin = (towerHit >> 48) & 0x000000000000FFFFLL;

      // phi bins for given eta bin
      phiBins = fPhiBins[etaBin];

      // calculate eta and phi of the tower's center
      fTowerEta = 0.5*(fEtaBins[etaBin - 1] + fEtaBins[etaBin]);
      fTowerPhi = 0.5*((*phiBins)[phiBin - 1] + (*phiBins)[phiBin]);

      // take the edges of the tower
      fTowerEdges[0] = fEtaBins[etaBin - 1];
      fTowerEdges[1] = fEtaBins[etaBin];
      fTowerEdges[2] = (*phiBins)[phiBin - 1];
      fTowerEdges[3] = (*phiBins)[phiBin];

      fTowerECalEnergy = 0.0;
      fTowerHCalEnergy = 0.0;

      fTrackECalEnergy = 0.0;
      fTrackHCalEnergy = 0.0;

      fTowerECalTime = 0.0;
      fTowerHCalTime = 0.0;

      fTrackECalTime = 0.0;
      fTrackHCalTime = 0.0;

      fTowerECalTimeWeight = 0.0;
      fTowerHCalTimeWeight = 0.0;

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

      fTrackECalTime += TMath::Sqrt(ecalEnergy)*position.T();
      fTrackHCalTime += TMath::Sqrt(hcalEnergy)*position.T();

      fTrackECalTimeWeight += TMath::Sqrt(ecalEnergy);
      fTrackHCalTimeWeight += TMath::Sqrt(hcalEnergy);

      fTowerTrackArray->Add(track);

      continue; // go to the next hit
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

    fTowerECalTime += TMath::Sqrt(ecalEnergy)*position.T();
    fTowerHCalTime += TMath::Sqrt(hcalEnergy)*position.T();

    fTowerECalTimeWeight += TMath::Sqrt(ecalEnergy);
    fTowerHCalTimeWeight += TMath::Sqrt(hcalEnergy);


    fTower->PID = fParticlePDGId.at(number);

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
  Double_t ecalEfficiency, hcalEfficiency;
  Double_t ecalTime, hcalTime, time;

  if(!fTower) return;

  energy = 0;

  // compute the new ecal energy
  ecalSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy); // take ecal resolution
  ecalEnergy = LogNormal(fTowerECalEnergy, ecalSigma);

  // compute the ecal tower efficiency (in case of non gaussian resolution, for aged scenario as example)
  ecalEfficiency = fECalEfficiencyFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy);
  bool isEcalGood = gRandom->Uniform(0.,1.) < ecalEfficiency ;

  // compute the new hcal energy
  hcalSigma  = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy);
  hcalEnergy = LogNormal(fTowerHCalEnergy, hcalSigma);

  // compute the hcal tower efficiency (in case of non gaussian resolution, for aged scenario as example)
  hcalEfficiency = fHCalEfficiencyFormula->Eval(0.0, fTowerEta, 0.0, fTowerHCalEnergy);
  bool isHcalGood = gRandom->Uniform(0.,1.) < hcalEfficiency ;

  ecalTime = (fTowerECalTimeWeight < 1.0E-09 ) ? 0.0 : fTowerECalTime/fTowerECalTimeWeight;
  hcalTime = (fTowerHCalTimeWeight < 1.0E-09 ) ? 0.0 : fTowerHCalTime/fTowerHCalTimeWeight;

  // take the resoluton for the new real ecal energy
  ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, ecalEnergy);
  hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, hcalEnergy);

  if(ecalEnergy < fECalEnergyMin || ecalEnergy < fECalEnergySignificanceMin*ecalSigma) ecalEnergy = 0.0;
  if(hcalEnergy < fHCalEnergyMin || hcalEnergy < fHCalEnergySignificanceMin*hcalSigma) hcalEnergy = 0.0;

  time = (TMath::Sqrt(ecalEnergy)*ecalTime + TMath::Sqrt(hcalEnergy)*hcalTime)/(TMath::Sqrt(ecalEnergy) + TMath::Sqrt(hcalEnergy));

  // add ecal and hcal to the total energy and save the tower information
  if( isEcalGood ) 
    energy     += ecalEnergy;
  if( isHcalGood ) 
    energy += hcalEnergy;

  eta = fTowerEta;
  phi = fTowerPhi;

  if(fSmearTowerCenter){
      eta = gRandom->Uniform(fTowerEdges[0], fTowerEdges[1]);
      phi = gRandom->Uniform(fTowerEdges[2], fTowerEdges[3]);
  }
  else{
      eta = fTowerEta;
      phi = fTowerPhi;
  }

  pt = energy / TMath::CosH(eta);

  fTower->Position.SetPtEtaPhiE(1.0, eta, phi, time);
  fTower->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
  fTower->Eem = ecalEnergy;
  fTower->Ehad = hcalEnergy;

  fTower->Edges[0] = fTowerEdges[0];
  fTower->Edges[1] = fTowerEdges[1];
  fTower->Edges[2] = fTowerEdges[2];
  fTower->Edges[3] = fTowerEdges[3];


  // fill calorimeter towers and photon candidates
  if(energy > 0.0){
    if(fTowerPhotonHits > 0 && fTowerTrackHits == 0)
      fPhotonOutputArray->Add(fTower);            
    fTowerOutputArray->Add(fTower);
  }

  // save all the tracks as energy flow tracks
  fItTowerTrackArray->Reset();
  while((track = static_cast<Candidate*>(fItTowerTrackArray->Next()))){
    fEFlowTrackOutputArray->Add(track);
  }

  

  // apply selective readout threshold before building the input set of particles to the Delphes PF algo
  if(ecalEnergy < fECalEnergyMin) ecalEnergy = 0.0;
  if(hcalEnergy < fHCalEnergyMin) hcalEnergy = 0.0;

  // subtract the contribution of charge particles measured by tracking system in order to split correctly
  ecalEnergy -= fTrackECalEnergy;
  if(ecalEnergy < 0.0) ecalEnergy = 0.0;

  hcalEnergy -= fTrackHCalEnergy;
  if(hcalEnergy < 0.0) hcalEnergy = 0.0;

  // recompute the energy resolution of the "tower"
  ecalSigma = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, ecalEnergy);
  hcalSigma = fHCalResolutionFormula->Eval(0.0, fTowerEta, 0.0, hcalEnergy);

  if(ecalEnergy < fECalEnergySignificanceMin*ecalSigma)
    ecalEnergy = 0.0;
  
  if(hcalEnergy < fHCalEnergySignificanceMin*hcalSigma)
    hcalEnergy = 0.0;
  
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
