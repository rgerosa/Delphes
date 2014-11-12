
/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
 *
 *  $Date: 2013-11-04 11:59:27 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1315 $
 *
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \with modifications by J. Stupak and J. Dolen
 *
 */

#include "modules/FastJetFinder.h"

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
#include <vector>

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

#include "fastjet/SISConePlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"

using namespace std;
using namespace fastjet;
using namespace contrib;
//------------------------------------------------------------------------------

FastJetFinder::FastJetFinder() :
  fPlugin(0), fDefinition(0), fAreaDefinition(0), fItInputArray(0)
{

}

//------------------------------------------------------------------------------

FastJetFinder::~FastJetFinder()
{

}

//------------------------------------------------------------------------------

void FastJetFinder::Init()
{
  JetDefinition::Plugin *plugin = NULL;

  // read eta ranges

  ExRootConfParam param = GetParam("RhoEtaRange");
  Long_t i, size;

  fEtaRangeMap.clear();
  size = param.GetSize();
  for(i = 0; i < size/2; ++i)
  {
    fEtaRangeMap[param[i*2].GetDouble()] = param[i*2 + 1].GetDouble();
  }

  fKeepPileUp = GetInt("KeepPileUp",1);
  //  cout << "set fKeepPileUp=" << fKeepPileUp << endl;

  // define algorithm

  fJetAlgorithm = GetInt("JetAlgorithm", 6);
  fParameterR = GetDouble("ParameterR", 0.5);

  fConeRadius = GetDouble("ConeRadius", 0.5);
  fSeedThreshold = GetDouble("SeedThreshold", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fMaxIterations = GetInt("MaxIterations", 100);
  fMaxPairSize = GetInt("MaxPairSize", 2);
  fIratch = GetInt("Iratch", 1);
  fAdjacencyCut = GetDouble("AdjacencyCut", 2.0);
  fOverlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetPTMin = GetDouble("JetPTMin", 10.0);

  // ---  Jet Area Parameters ---
  fAreaAlgorithm = GetInt("AreaAlgorithm", 0);
  fComputeRho = GetBool("ComputeRho", false);
  // - ghost based areas -
  fGhostEtaMax = GetDouble("GhostEtaMax", 5.0);
  fRepeat = GetInt("Repeat", 1);
  fGhostArea = GetDouble("GhostArea", 0.01);
  fGridScatter = GetDouble("GridScatter", 1.0);
  fPtScatter = GetDouble("PtScatter", 0.1);
  fMeanGhostPt = GetDouble("MeanGhostPt", 1.0E-100);
  // - voronoi based areas -
  fEffectiveRfact = GetDouble("EffectiveRfact", 1.0);

  switch(fAreaAlgorithm)
  {
    case 1:
      fAreaDefinition = new fastjet::AreaDefinition(active_area_explicit_ghosts, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 2:
      fAreaDefinition = new fastjet::AreaDefinition(one_ghost_passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 3:
      fAreaDefinition = new fastjet::AreaDefinition(passive_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    case 4:
      fAreaDefinition = new fastjet::AreaDefinition(VoronoiAreaSpec(fEffectiveRfact));
      break;
    case 5:
      fAreaDefinition = new fastjet::AreaDefinition(active_area, GhostedAreaSpec(fGhostEtaMax, fRepeat, fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
    default:
    case 0:
      fAreaDefinition = 0;
      break;
  }

  switch(fJetAlgorithm)
  {
    case 1:
      plugin = new fastjet::CDFJetCluPlugin(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 2:
      plugin = new fastjet::CDFMidPointPlugin(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 3:
      plugin = new fastjet::SISConePlugin(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 4:
      fDefinition = new fastjet::JetDefinition(fastjet::kt_algorithm, fParameterR);
      break;
    case 5:
      fDefinition = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fParameterR);
      break;
    default:
    case 6:
      fDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fParameterR);
      break;
  }

  fPlugin = plugin;

  ClusterSequence::print_banner();

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays

  fOutputArray = ExportArray(GetString("OutputArray", "jets"));
  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
}

//------------------------------------------------------------------------------

void FastJetFinder::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fDefinition) delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fPlugin) delete static_cast<JetDefinition::Plugin*>(fPlugin);
}

//------------------------------------------------------------------------------

void FastJetFinder::Process()
{
  // cout << " Starting FastJetFinder::Process() " << GetString("InputArray", "Calorimeter/towers") << endl;

  Candidate *candidate, *constituent;
  TLorentzVector momentum;
  Double_t deta, dphi, detaMax, dphiMax;
  Int_t number;
  Double_t rho = 0;
  PseudoJet jet, area;
  vector<PseudoJet> inputList, outputList;
  ClusterSequence *sequence;
  map< Double_t, Double_t >::iterator itEtaRangeMap;

  DelphesFactory *factory = GetFactory();

  inputList.clear();

  // loop over input objects
  fItInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {

    // for genjets mostly
    if (fKeepPileUp == 0 && candidate->IsPU > 0) {
    //      cout << " skipping pile up object" << endl;
      continue;
    }

    momentum = candidate->Momentum;
    jet = PseudoJet(momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;

    // cout << " fastJet pt eta phi m " << jet.pt() << " " << jet.eta() << " " << jet.phi() << " " << jet.m() << endl;
  }

  // construct jets
  if(fAreaDefinition)
  {
    sequence = new ClusterSequenceArea(inputList, *fDefinition, *fAreaDefinition);
  }
  else
  {
    sequence = new ClusterSequence(inputList, *fDefinition);
  }

  // compute rho and store it
  if(fComputeRho && fAreaDefinition)
  {
    for(itEtaRangeMap = fEtaRangeMap.begin(); itEtaRangeMap != fEtaRangeMap.end(); ++itEtaRangeMap)
    {
      Selector select_rapidity = SelectorAbsRapRange(itEtaRangeMap->first, itEtaRangeMap->second);
      JetMedianBackgroundEstimator estimator(select_rapidity, *fDefinition, *fAreaDefinition);
      estimator.set_particles(inputList);
      rho = estimator.rho();

      candidate = factory->NewCandidate();
      candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
      candidate->Edges[0] = itEtaRangeMap->first;
      candidate->Edges[1] = itEtaRangeMap->second;
      fRhoOutputArray->Add(candidate);
    }
  }

  outputList.clear();
  outputList = sorted_by_pt(sequence->inclusive_jets(fJetPTMin));

  // loop over all jets and export them
  detaMax = 0.0;
  dphiMax = 0.0;
  vector<PseudoJet>::iterator itInputList, itOutputList;
  for(itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList)
  {
    momentum.SetPxPyPzE(itOutputList->px(), itOutputList->py(), itOutputList->pz(), itOutputList->E());
    area.reset(0.0, 0.0, 0.0, 0.0);
    if(fAreaDefinition) area = itOutputList->area_4vector();

    candidate = factory->NewCandidate();

    inputList.clear();
    inputList = sequence->constituents(*itOutputList);
    for(itInputList = inputList.begin(); itInputList != inputList.end(); ++itInputList)
    {
      constituent = static_cast<Candidate*>(fInputArray->At(itInputList->user_index()));

      //      cout << "fKeepPileUp=" << fKeepPileUp << " constituent->IsPU=" << constituent->IsPU << endl;

      deta = TMath::Abs(momentum.Eta() - constituent->Momentum.Eta());
      dphi = TMath::Abs(momentum.DeltaPhi(constituent->Momentum));
      if(deta > detaMax) detaMax = deta;
      if(dphi > dphiMax) dphiMax = dphi;

      candidate->AddCandidate(constituent);
    }

    candidate->Momentum = momentum;
    candidate->Area.SetPxPyPzE(area.px(), area.py(), area.pz(), area.E());

    candidate->DeltaEta = detaMax;
    candidate->DeltaPhi = dphiMax;

    //------------------------------------
    // SubStructure
    //------------------------------------

    if (itOutputList->perp()>200){

      //------------------------------------
      // Trimming
      //------------------------------------

      double Rtrim = 0.2;
      double ptfrac = 0.05;
      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::cambridge_algorithm, Rtrim), fastjet::SelectorPtFractionMin(ptfrac) );
      fastjet::PseudoJet trimmed_jet = trimmer(*itOutputList);

      double trimmedMass=trimmed_jet.m();
      if (trimmedMass<0) trimmedMass=0;
      candidate->TrimmedMass=trimmedMass;
  
      //------------------------------------
      // Subjet Quantities
      //------------------------------------

      vector<PseudoJet> kept_subjets = trimmed_jet.pieces();
      candidate->NSubJets=kept_subjets.size();

      double largest_mass_subjet = 0;
      for (size_t i = 0; i < kept_subjets.size(); i++)
	{
	  if ( kept_subjets[i].m() > largest_mass_subjet )
	    {
	      largest_mass_subjet = kept_subjets[i].m();
	    }
	}

      double massdrop=1;
      if ( trimmedMass!=0 ) massdrop = largest_mass_subjet / trimmedMass;
      candidate->MassDrop=massdrop;

      //------------------------------------
      // NSubJettiness
      //------------------------------------
      double beta = 1.0; // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
      double R0   = 0.8; // Characteristic jet radius for normalization
      double Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff)
      
      fastjet::contrib::Nsubjettiness nSub1(1,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      fastjet::contrib::Nsubjettiness nSub2(2,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      fastjet::contrib::Nsubjettiness nSub3(3,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      candidate->Tau1 = nSub1(*itOutputList);
      candidate->Tau2 = nSub2(*itOutputList);
      candidate->Tau3 = nSub3(*itOutputList);

      //------------------------------------
      // W-tag - Mass drop from trimmed jets
      //------------------------------------

      if (massdrop < 0.4 && 60 < trimmedMass && trimmedMass < 120)
	candidate->WTag=1;

      //------------------------------------
      // Top-tag
      //------------------------------------

      if (kept_subjets.size()>=3 && 140 < trimmedMass && trimmedMass < 230)
	candidate->TopTag=1;

      //------------------------------------
      // h-tag - Mass drop from trimmed jets
      //------------------------------------

      if (massdrop < 0.4 && 100 < trimmedMass && trimmedMass < 140)
	candidate->HTag=1;
    }

    fOutputArray->Add(candidate);
  }
  delete sequence;

  // cout << " Ending FastJetFinder::Process() " << endl;
}
