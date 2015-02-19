/** \class ModifyBeamSpot
 *  \author S. Zenz
*/

#include "modules/ModifyBeamSpot.h"

//#include "CLHEP/Units/GlobalSystemOfUnits.h"
//#include "CLHEP/Units/GlobalPhysicalConstants.h"

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

// static const double mm  = 1.;
// static const double m = 1000.*mm;
// static const double ns  = 1.;
// static const double s = 1.e+9 *ns;
// static const double c_light   = 2.99792458e+8 * m/s;

static const double c_light     = 0.299792458 ; // mm/ps , or m/ns
static const double inv_c_light = 3.335640952 ; // ps/mm or ns/m


using namespace std;

//------------------------------------------------------------------------------

ModifyBeamSpot::ModifyBeamSpot() :
  fFormula(0), fItInputArray(0)
{
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

ModifyBeamSpot::~ModifyBeamSpot()
{
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Init()
{
  // read resolution formula

  fZVertexSpread = GetDouble ("ZVertexSpread", 53)  ;
  // mm in the cfg file, mm in the code
  fTVertexSpread = GetDouble ("TVertexSpread", 160) ;
  // ps in the cfg file, ps in the code

  fOutputBSX = GetDouble ("OutputBSX", 0.) ;
  fOutputBSY = GetDouble ("OutputBSY", 0.) ;
  fOutputBSZ = GetDouble ("OutputBSZ", 0.) ;

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
  fPVOutputArray = ExportArray(GetString("PVOutputArray", "PV"));
  
  fDebugOutputCollector.addVariable ("PUindex") ;
  
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Finish()
{  
  std::string outfile = GetString ("simpleOutputFileName", "simpleOutput_Ca.root") ;
  fDebugOutputCollector.save (outfile) ;

  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Process()
{
  Candidate *candidate, *mother;
  // Average position of primary particles
  Double_t PVX = 0., PVY = 0., PVZ = 0., PVT = 0.; 

  PVX = fOutputBSX ;
  PVY = fOutputBSY ;
  PVZ = gRandom->Gaus (0., fZVertexSpread) + fOutputBSZ ;
  PVT = gRandom->Gaus (0., fTVertexSpread * c_light) ;

  fItInputArray->Reset () ;
  // loop over particles in the event
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;

    double X = candidatePosition.X () ;
    double Y = candidatePosition.Y () ;
    double Z = candidatePosition.Z () ;
    double T = candidatePosition.T () ;

    // LHE files are centered in (0,0,0,0)
    // while PU has been generated with smearing already, therefore I smear only 
    // the hard scattering.
    if (candidate->IsPU == 0)
      {
        X += PVX ;
        Y += PVY ;
        Z += PVZ ;
        T += PVT ;
      }

    fDebugOutputCollector.fillContainer ("PUindex", candidate->IsPU) ;

    mother = candidate ;
    candidate = static_cast<Candidate*> (candidate->Clone ()) ;
    candidate->Position.SetXYZT (X, Y, Z, T) ;
    candidate->AddCandidate (mother) ;
    fOutputArray->Add (candidate) ;

  } // loop over particles in the event

  // Store the PV "beam spot"
  DelphesFactory *factory = GetFactory () ;
  candidate               = factory->NewCandidate () ;
  candidate->Position.SetXYZT (PVX, PVY, PVZ, PVT) ;
  fPVOutputArray->Add (candidate) ;

  ++fEventCounter ;                                                                                                                                             

}

//------------------------------------------------------------------------------
