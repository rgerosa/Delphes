/** \class NeutrinoFilter
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

#include "modules/NeutrinoFilter.h"

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
//#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

NeutrinoFilter::NeutrinoFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

NeutrinoFilter::~NeutrinoFilter()
{
}

//------------------------------------------------------------------------------

void NeutrinoFilter::Init()
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

void NeutrinoFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void NeutrinoFilter::Process()
{
  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    pdgCode = TMath::Abs(candidate->PID);

    pass = kTRUE;

    if (pdgCode == 12 || pdgCode == 14 || pdgCode == 16) pass = kFALSE;

    if (pass) fOutputArray->Add(candidate);
  }
}

