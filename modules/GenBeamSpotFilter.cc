/** \class GenBeamSpotFilter
 *
 *  Removes all generated particles except electrons, muons, taus,
 *  and particles with status == 3.
 *
 *  $Date: 2013-05-20 22:22:07 +0200 (Mon, 20 May 2013) $
 *  $Revision: 1118 $
 *
 *
 *  \author J. Hirschauer - FNAL
 *
 */

#include "modules/GenBeamSpotFilter.h"

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

GenBeamSpotFilter::GenBeamSpotFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

GenBeamSpotFilter::~GenBeamSpotFilter()
{
}

//------------------------------------------------------------------------------

void GenBeamSpotFilter::Init()
{
  // PT threshold

  fPTMin = GetDouble("PTMin", 0.5);

  // import input array
  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void GenBeamSpotFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void GenBeamSpotFilter::Process()
{
  Candidate *candidate;
  Bool_t pass;
  TLorentzVector origin, firstOrigin;
  int N_same = 0;
  int N_different = 0 ;
  int N_pileup = 0;
  fPassedOne = false;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    origin = candidate->Position;
    pass = kFALSE;

    if (candidate->IsPU == 0 && !fPassedOne) {
      firstOrigin = origin;
      pass = true;
      fPassedOne = true;
      N_same = 0;
      N_different = 0;
      N_pileup = 0;
    }

    //    cout << candidate->IsPU << " " << origin.X() << " "<< origin.Y() << " " << origin.Z() << endl;
    if (candidate->IsPU == 0) {
      if (origin == firstOrigin) {
	//	cout << "From firstOrigin: " << candidate->IsPU << " " << origin.X() << " "<< origin.Y() << " " << origin.Z() << endl;
	N_same++;
      } else {
	//	cout << "NOT from firstOrigin: " << candidate->IsPU << " " << origin.X() << " "<< origin.Y() << " " << origin.Z() << endl;
	N_different++;
      }
    } else {
      N_pileup++;
    }

    if (pass) fOutputArray->Add(candidate);
  }

  //  cout << "GenBeamSpotFilter Total same/different/pileup in this processing: " << N_same << "/" << N_different << "/" << N_pileup << endl;
}

