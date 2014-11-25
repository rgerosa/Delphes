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

static const double mm  = 1.;
static const double m = 1000.*mm;
static const double ns  = 1.;
static const double s = 1.e+9 *ns;
static const double c_light   = 2.99792458e+8 * m/s;


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

  fZVertexSpread = GetDouble("ZVertexSpread", 0.05)*1.0E3;

  currentPU = -1;

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));

  fPVOutputArray = ExportArray(GetString("PVOutputArray", "PV"));
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Finish()
{  
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void ModifyBeamSpot::Process()
{
  Candidate *candidate, *mother;
  Double_t PVX = 0., PVY = 0., PVZ = 0., PVT = 0.; // Average position of primary particles
  Int_t PVN = 0;
  currentPU = -1;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    if (currentPU < 0 || currentPU != candidate->IsPU) {
      //cout << "SCZ Debug: currentPU=" << currentPU << " candidate->IsPU=" << candidate->IsPU << " so throwing new numbers" << endl;

      if (currentPU == 0) {
	// Done with PV
	PVX /= PVN;
	PVY /= PVN;
	PVZ /= PVN;
	PVT /= PVN;
      }

      currentPU = candidate->IsPU; // Depends on generated particle IsPU being set to a different value for each vertex

      // N.B. Z and T are not correlated
      // Not exactly right but it seems to be what CMSSW does
      currentZ = gRandom->Gaus(0., fZVertexSpread);
      currentT = gRandom->Gaus(0., fZVertexSpread*(mm/ns)/c_light);
    }

    // For now simply keep X, Y - they looked okish
    double X = candidatePosition.X();
    double Y = candidatePosition.Y();
    double Z = currentZ;
    double T = currentT;

    if (currentPU == 0) {
      PVN++;
      PVX += X;
      PVY += Y;
      PVZ += Z;
      PVT += T;
    }
 
    mother = candidate;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Position.SetXYZT(X,Y,Z,T);
    candidate->AddCandidate(mother);
    fOutputArray->Add(candidate);

  }

  // If PV is somehow last (i.e. no pileup, still gotta divide out
  if (currentPU == 0) {
    // Done with PV
    PVX /= PVN;
    PVY /= PVN;
    PVZ /= PVN;
    PVT /= PVN;
  }

  //  cout << " SCZ ModifyBeamSpot PVZ=" << PVZ << endl;

  // Store the PV "beam spot"
  DelphesFactory *factory = GetFactory();
  candidate               = factory->NewCandidate();
  candidate->Position.SetXYZT(PVX,PVY,PVZ,PVT);
  fPVOutputArray->Add(candidate);

}

//------------------------------------------------------------------------------
