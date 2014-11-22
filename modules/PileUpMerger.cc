/** \class PileUpMerger
 *
 *  Merges particles from pile-up sample into event
 *
 *
 *  $Date: 2013-02-12 15:13:59 +0100 (Tue, 12 Feb 2013) $
 *  $Revision: 907 $
 *
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/PileUpMerger.h"

//#include "CLHEP/Units/GlobalSystemOfUnits.h"
//#include "CLHEP/Units/GlobalPhysicalConstants.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"
#include "classes/DelphesPileUpReader.h"

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

static const double mm  = 1.;
static const double m = 1000.*mm;
static const double ns  = 1.;
static const double s = 1.e+9 *ns;
static const double c_light   = 2.99792458e+8 * m/s;


using namespace std;

//------------------------------------------------------------------------------

PileUpMerger::PileUpMerger() :
  fReader(0), fItInputArray(0)
{
}

//------------------------------------------------------------------------------

PileUpMerger::~PileUpMerger()
{
}

//------------------------------------------------------------------------------

void PileUpMerger::Init()
{
  const char *fileName;

  fMeanPileUp  = GetDouble("MeanPileUp", 10);
  fZVertexSpread = GetDouble("ZVertexSpread", 0.05)*1.0E3;

  fInputBSX = GetDouble("InputBSX",0.);
  fInputBSY = GetDouble("InputBSY",0.);
  fOutputBSX = GetDouble("OutputBSX",0.);
  fOutputBSY = GetDouble("OutputBSY",0.);


  fileName = GetString("PileUpFile", "MinBias.pileup");
  fReader = new DelphesPileUpReader(fileName);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fNPUOutputArray = ExportArray(GetString("NPUOutputArray", "NPU"));
}

//------------------------------------------------------------------------------

void PileUpMerger::Finish()
{
  if(fReader) delete fReader;
}

//------------------------------------------------------------------------------

void PileUpMerger::Process()
{
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  TParticlePDG *pdgParticle;
  Int_t pid;
  Float_t x, y, z, t;
  Float_t px, py, pz, e;
  Double_t dz, dt, dphi;
  Int_t poisson, event;
  Long64_t allEntries, entry;
  Candidate *candidate;
  DelphesFactory *factory;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    fOutputArray->Add(candidate);
  }

  factory = GetFactory();

  poisson = gRandom->Poisson(fMeanPileUp);

  allEntries = fReader->GetEntries();

  for(event = 0; event < poisson; ++event){
    do{
      entry = TMath::Nint(gRandom->Rndm()*allEntries);
    }
    while(entry >= allEntries);

    fReader->ReadEntry(entry);

    dz   = gRandom->Gaus(0.0, fZVertexSpread);
    dphi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
    dt   = gRandom->Gaus(0., fZVertexSpread*(mm/ns)/c_light);
 
    while(fReader->ReadParticle(pid, x, y, z, t, px, py, pz, e)){

      candidate = factory->NewCandidate();

      // Get rid of BS position in PU
      // To deal with http://red-gridftp11.unl.edu/Snowmass/MinBias100K_14TeV.pileup fInputBSX = 2.44 and fInputBSY = 3.93
      x = x - fInputBSX;
      y = y - fInputBSY;

      candidate->PID    = pid;
      candidate->Status = 1;

      pdgParticle = pdg->GetParticle(pid);
      candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge()/3.0) : -999;
      candidate->Mass = pdgParticle ? pdgParticle->Mass() : -999.9;

      candidate->IsPU = event+1; // might as well store which PU vertex this comes from so they can be separated

      candidate->Momentum.SetPxPyPzE(px, py, pz, e);
      candidate->Momentum.RotateZ(dphi);

      candidate->Position.SetXYZT(x,y,dz,dt); // Use CMSSW dz,dt (ignoring old values), keep x y
      // Kept for backward compatibility - ModifyBeamSpot will change them again, it's needed to do the primary event also

      candidate->Position.RotateZ(dphi);

      // CMSSW uses fOutputBSX = 0.24 and fOutputBSY = 0.39 in the files I looked at
      candidate->Position += TLorentzVector(fOutputBSX,fOutputBSY,0.,0.);

      fOutputArray->Add(candidate);
    }
  }

  // Store true number of pileup vertices
  candidate = factory->NewCandidate();
  candidate->Momentum.SetPtEtaPhiE((float)poisson, 0.0, 0.0, (float)poisson); // cheating and storing NPU as a float
  fNPUOutputArray->Add(candidate);

}

//------------------------------------------------------------------------------

