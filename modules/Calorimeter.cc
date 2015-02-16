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

inline bool isMatching (TLorentzVector & part, vector<TLorentzVector> & set) 
{
  for (unsigned int i = 0 ; i < set.size () ; ++i)
    {
      if (part.DeltaR (set.at (i)) < 0.1) return true ;
    }
  return false ;
}

//------------------------------------------------------------------------------

Calorimeter::Calorimeter() :
  fECalResolutionFormula(0), fHCalResolutionFormula(0),
  fItParticleInputArray(0), fItTrackInputArray(0),
  fTowerTrackArray(0), fItTowerTrackArray(0),
  fItLHEPartonInputArray (0),
  fDelayBarrel ("fDelayBarrel", "pol5", 0, 5),
  fDelayEndcap ("fDelayEndcap", "pol5", 0, 5)
  {

  fECalResolutionFormula = new DelphesFormula;
  fHCalResolutionFormula = new DelphesFormula;

  fTowerTrackArray = new TObjArray;
  fItTowerTrackArray = fTowerTrackArray->MakeIterator();
}

//------------------------------------------------------------------------------

Calorimeter::~Calorimeter(){

  if(fECalResolutionFormula) delete fECalResolutionFormula;
  if(fHCalResolutionFormula) delete fHCalResolutionFormula;
  if(fTowerTrackArray)       delete fTowerTrackArray;
  if(fItTowerTrackArray)     delete fItTowerTrackArray;

  if(fItParticleInputArray)  delete fItParticleInputArray;
  if(fItTrackInputArray)     delete fItTrackInputArray;
  if(fItLHEPartonInputArray) delete fItLHEPartonInputArray;

  vector< vector< Double_t >* >::iterator itPhiBin;
  for (itPhiBin = fPhiBins.begin(); itPhiBin != fPhiBins.end(); ++itPhiBin){
    delete *itPhiBin;
  }
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

  //PG for matching to LHE particles, for timing studies
  fLHEPartonInputArray = ImportArray(GetString("LHEPartonInputArray", "Delphes/LHEParticles"));
  fItLHEPartonInputArray = fLHEPartonInputArray->MakeIterator();

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
  // since the curved extrapolation not supported
  // if this value is true, electron timing is taken from the track,
  // otherwise is taken from the particle collection
  fElectronsFromTrack  = false;

  Double_t delayBarrelParams[6] ;
  delayBarrelParams[0] = GetDouble ("DelayBarrel_0",  1257.22) ;
  delayBarrelParams[1] = GetDouble ("DelayBarrel_1",  436.655) ;
  delayBarrelParams[2] = GetDouble ("DelayBarrel_2", -444.181) ;
  delayBarrelParams[3] = GetDouble ("DelayBarrel_3",  961.054) ;
  delayBarrelParams[4] = GetDouble ("DelayBarrel_4", -235.284) ;
  delayBarrelParams[5] = GetDouble ("DelayBarrel_5", -235.284) ;
  fDelayBarrel.SetParameters (delayBarrelParams) ;

  Double_t delayEndcapParams[6] ;
  delayEndcapParams[0] = GetDouble ("DelayEndcap_0",  4494.45) ;
  delayEndcapParams[1] = GetDouble ("DelayEndcap_1", -1525.16) ;
  delayEndcapParams[2] = GetDouble ("DelayEndcap_2",  619.747) ;
  delayEndcapParams[3] = GetDouble ("DelayEndcap_3", -114.044) ;
  delayEndcapParams[4] = GetDouble ("DelayEndcap_4",  7.84058) ;
  delayEndcapParams[5] = GetDouble ("DelayEndcap_5",  7.84058) ;
  fDelayEndcap.SetParameters (delayEndcapParams) ;

  //PG FIXME where does this come from?
  // suggested from A. Bornheim, reasonable according to him
  fTimingEMin = GetDouble ("TimingEMin", 4.) ;

  //simple outputs during running
  fDebugOutputCollector.addVariable ("DR") ;
  fDebugOutputCollector.addVariable ("partType") ;
  fDebugOutputCollector.addVariable ("m_partType") ;
  fDebugOutputCollector.addVariable ("Nm_partType") ;
  fDebugOutputCollector.addVariable3D ("eta:time:pt") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time:pt") ;
  fDebugOutputCollector.addVariable3D ("m_eta:time_nc:pt") ;
  fDebugOutputCollector.addVariable3D ("Nm_eta:time:pt") ;
  fDebugOutputCollector.addNtuple ("eta:time:pt:PID") ;
  fEventCounter = 0 ;

}

//------------------------------------------------------------------------------

void Calorimeter::Finish(){

  std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_Ca.root") ;
  fDebugOutputCollector.save (outfile) ;
}

//------------------------------------------------------------------------------

void Calorimeter::Process(){

  Candidate *particle, *track;
  TLorentzVector position, momentum;
  Short_t etaBin, phiBin, flags;
  Int_t number;
  Long64_t towerHit, towerEtaPhi, hitEtaPhi;
  Double_t ecalFraction, hcalFraction;
  Double_t ecalEnergy, hcalEnergy, totalEnergy;
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

  // get the eta, phi, pt of the LHE partons that could generate
  // a Delphes jet
  vector<TLorentzVector> LHEParticles ;
  fItLHEPartonInputArray->Reset () ;
  while (Candidate * LHEparticle = static_cast<Candidate*> (fItLHEPartonInputArray->Next ()))
    {
      if (LHEparticle->Status != 1) continue ;
      if (fabs (LHEparticle->PID) == 11 ||   // electron
//          1) 
          fabs (LHEparticle->PID) < 7   ||   // quarks
          fabs (LHEparticle->PID) == 21 ||   // gluon
          fabs (LHEparticle->PID) == 21)     // photon
        LHEParticles.push_back (LHEparticle->Momentum) ;
    }

  // loop over all tracks to get the deposited energy due to the
  // charged hadrons and electrons. 
  // THis energy is added to the tower as additional information,
  // that will be used in FinalizeTower to mimic the particle flow
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
 
    // make tower hit:
    //    {16-bits for eta bin number, 
    //     16-bits for phi bin number, 
    //     8-bits for flags, 
    //     24-bits for track number}
    towerHit = (Long64_t(etaBin) << 48) | (Long64_t(phiBin) << 32) | (Long64_t(flags) << 24) | Long64_t(number);

    fTowerHits.push_back(towerHit);
  } // loop over all tracks

  // loop over all particles of the event,
  // to get the energy of all the particles that hit the calorimeter
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
  } // loop over all particles of the event

  // all hits are sorted first by eta bin number, then by phi bin number,
  // then by flags and then by particle or track number
  sort(fTowerHits.begin(),fTowerHits.end()); // sort the hit

  // loop over all hits (tracks and calo deposits)
  towerEtaPhi = 0;
  fTower = 0;
  
  // loop over the hits
  for(itTowerHits = fTowerHits.begin(); itTowerHits != fTowerHits.end(); ++itTowerHits) { 
    towerHit  = (*itTowerHits);
    flags     = (towerHit >> 24) & 0x00000000000000FFLL; // take the flag
    number    = (towerHit) & 0x0000000000FFFFFFLL;       // take the number
    hitEtaPhi = towerHit >> 32;
    if(towerEtaPhi != hitEtaPhi){  
      // first time hit, no  more than once since we have sorted tower 
      // as a function of eta and phi hits

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

      fTowerTrackHits = 0;
      fTowerPhotonHits = 0;

      fTowerTrackArray->Clear();
    }  // first time hit

    // check for track hits
    if(flags & 1){
      ++fTowerTrackHits;
      track = static_cast<Candidate*>(fTrackInputArray->At(number));
      momentum = track->Momentum;

      bool dbg_scz = false;
      if (dbg_scz) {
        cout << "   Calorimeter input track has x y z t " << track->Position.X() 
             << " " << track->Position.Y() << " " << track->Position.Z() << " " << track->Position.T() 
             << endl;
        Candidate *prt = static_cast<Candidate*>(track->GetCandidates()->Last());
        const TLorentzVector &ini = prt->Position;
        cout << "                and parent has x y z t " << ini.X() << " " << ini.Y() 
             << " " << ini.Z() << " " << ini.T();
      }

      ecalEnergy = momentum.E() * fTrackECalFractions[number];
      hcalEnergy = momentum.E() * fTrackHCalFractions[number];

      fTrackECalEnergy += ecalEnergy;
      fTrackHCalEnergy += hcalEnergy;

      totalEnergy = ecalEnergy + hcalEnergy ;

      if ( totalEnergy > fTimingEMin && fTower && fElectronsFromTrack) {
        float delay = 0. ;
        float feta = fabs (track->Position.Eta ()) ;
        if (feta < 1.6) delay = fDelayBarrel.Eval (feta) ;
        else            delay = fDelayEndcap.Eval (feta) ;
        
        fTower->ecal_E_t.push_back(
          std::make_pair<float,float>(totalEnergy, track->Position.T() - delay));
      }
      
      fTowerTrackArray->Add(track);

      continue; // go to the next hit
    } // check for track hits

    // check for photon and electron hits in current tower
    if(flags & 2) ++fTowerPhotonHits;

    particle = static_cast<Candidate*>(fParticleInputArray->At(number));
    momentum = particle->Momentum;

    // fill current tower
    ecalEnergy = momentum.E() * fTowerECalFractions[number];
    hcalEnergy = momentum.E() * fTowerHCalFractions[number];

    fTowerECalEnergy += ecalEnergy;
    fTowerHCalEnergy += hcalEnergy;

    totalEnergy = ecalEnergy + hcalEnergy ;

    if ( (totalEnergy > fTimingEMin && fTower) &&
         (abs(particle->PID) != 11 || !fElectronsFromTrack) ) {

        float delay = 0. ;
        float feta = fabs (particle->Position.Eta ()) ;
        if (feta < 1.6) delay = fDelayBarrel.Eval (feta) ;
        else            delay = fDelayEndcap.Eval (feta) ;
        
        fTower->ecal_E_t.push_back(
          std::make_pair<float,float>(totalEnergy, particle->Position.T() - delay));

        double Rmin = 100. ; 
        for (unsigned int i = 0 ; i < LHEParticles.size () ; ++i)
          {
            if (particle->Position.DeltaR (LHEParticles.at (i)) < Rmin) 
              Rmin = particle->Position.DeltaR (LHEParticles.at (i)) ;
          }
        fDebugOutputCollector.fillContainer ("DR", Rmin) ;
        fDebugOutputCollector.fillContainer ("partType", particle->PID) ;

        fDebugOutputCollector.fillContanier3D ("eta:time:pt", 
           feta, particle->Position.T () - delay, particle->Momentum.Pt ()) ;

        fDebugOutputCollector.getNtuple ("eta:time:pt:PID")->Fill 
          (
            feta,
            particle->Position.T () - delay,
            particle->Momentum.Pt (),
            particle->PID
          ) ;

        if (isMatching (particle->Position, LHEParticles)) 
          {
            fDebugOutputCollector.fillContainer ("m_partType", particle->PID) ;
            fDebugOutputCollector.fillContanier3D ("m_eta:time:pt", 
               feta, particle->Position.T () - delay, particle->Momentum.Pt ()) ;
            fDebugOutputCollector.fillContanier3D ("m_eta:time_nc:pt", 
               feta, particle->Position.T (), particle->Momentum.Pt ()) ;
          }
        else
          {
            fDebugOutputCollector.fillContainer ("Nm_partType", particle->PID) ;
            fDebugOutputCollector.fillContanier3D ("Nm_eta:time:pt", 
               feta, particle->Position.T () - delay, particle->Momentum.Pt ()) ;
          }
    }

    fTower->PID = fParticlePDGId.at(number);
    fTower->AddCandidate(particle);
  } // loop over the hits

  // finalize last tower
  FinalizeTower();
  ++fEventCounter ;
}

//------------------------------------------------------------------------------

void Calorimeter::FinalizeTower(){

  Candidate *track, *tower;
  Double_t energy, pt, eta, phi;
  Double_t ecalEnergy, hcalEnergy;
  Double_t ecalSigma, hcalSigma;

  if(!fTower) return;

  // ECAL resolution
  ecalSigma  = fECalResolutionFormula->Eval(0.0, fTowerEta, 0.0, fTowerECalEnergy);
  ecalEnergy = LogNormal(fTowerECalEnergy, ecalSigma);

  // HCAL resolution
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
