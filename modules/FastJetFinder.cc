/** \class FastJetFinder
 *  Finds jets using FastJet library.
 *  $Date: 2013-11-04 11:59:27 +0100 (Mon, 04 Nov 2013) $
 *  $Revision: 1315 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *  \with modifications by J. Stupak and J. Dolen
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
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/SoftDrop.hh"

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
  fPlugin(0), fDefinition(0), fAreaDefinition(0), fItInputArray(0){}

//------------------------------------------------------------------------------

FastJetFinder::~FastJetFinder(){}

//------------------------------------------------------------------------------
void FastJetFinder::Init(){

  JetDefinition::Plugin *plugin = NULL;

  // read eta ranges
  ExRootConfParam param = GetParam("RhoEtaRange");

  Long_t i, size;
  fEtaRangeMap.clear();
  size = param.GetSize();

  for(i = 0; i < size/2; ++i) fEtaRangeMap[param[i*2].GetDouble()] = param[i*2 + 1].GetDouble();  

  fKeepPileUp = GetInt("KeepPileUp",1);

  // define algorithm
  fJetAlgorithm = GetInt("JetAlgorithm", 6);

  fParameterR   = GetDouble("ParameterR", 0.5);

  fConeRadius       = GetDouble("ConeRadius", 0.5);
  fSeedThreshold    = GetDouble("SeedThreshold", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fMaxIterations    = GetInt("MaxIterations", 100);
  fMaxPairSize      = GetInt("MaxPairSize", 2);
  fIratch           = GetInt("Iratch", 1);
  fAdjacencyCut     = GetDouble("AdjacencyCut", 2.0);
  fOverlapThreshold = GetDouble("OverlapThreshold", 0.75);

  fJetPTMin = GetDouble("JetPTMin", 10.0);

  // ---  Jet Area Parameters ---
  fAreaAlgorithm  = GetInt("AreaAlgorithm", 0);
  fComputeRho     = GetBool("ComputeRho", false);
  fComputeRhoGrid = GetBool("ComputeRhoGrid", false);

  // - ghost based areas -
  fGhostEtaMax = GetDouble("GhostEtaMax", 5.0);
  fRepeat      = GetInt("Repeat", 1);
  fGhostArea   = GetDouble("GhostArea", 0.01);
  fGridScatter = GetDouble("GridScatter", 1.0);
  fPtScatter   = GetDouble("PtScatter", 0.1);
  fMeanGhostPt = GetDouble("MeanGhostPt", 1.0E-100);

  // - voronoi based areas -
  fEffectiveRfact = GetDouble("EffectiveRfact", 1.0);

  switch(fAreaAlgorithm){
    case 1:
      fAreaDefinition = new fastjet::AreaDefinition(active_area_explicit_ghosts,GhostedAreaSpec(fGhostEtaMax,fRepeat,fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
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
      fAreaDefinition = new fastjet::AreaDefinition(active_area_explicit_ghosts,GhostedAreaSpec(fGhostEtaMax,fRepeat,fGhostArea, fGridScatter, fPtScatter, fMeanGhostPt));
      break;
  }

  switch(fJetAlgorithm){
    case 1:
      plugin      = new fastjet::CDFJetCluPlugin(fSeedThreshold, fConeRadius, fAdjacencyCut, fMaxIterations, fIratch, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 2:
      plugin      = new fastjet::CDFMidPointPlugin(fSeedThreshold, fConeRadius, fConeAreaFraction, fMaxPairSize, fMaxIterations, fOverlapThreshold);
      fDefinition = new fastjet::JetDefinition(plugin);
      break;
    case 3:
      plugin      = new fastjet::SISConePlugin(fConeRadius, fOverlapThreshold, fMaxIterations, fJetPTMin);
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
  fInputArray   = ImportArray(GetString("InputArray", "Calorimeter/towers"));
  fItInputArray = fInputArray->MakeIterator();

  // create output arrays
  fOutputArray    = ExportArray(GetString("OutputArray", "jets"));
  fRhoOutputArray = ExportArray(GetString("RhoOutputArray", "rho"));
}

//------------------------------------------------------------------------------

void FastJetFinder::Finish()
{
  if(fItInputArray)   delete fItInputArray;
  if(fDefinition)     delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fPlugin) delete static_cast<JetDefinition::Plugin*>(fPlugin);
}

//------------------------------------------------------------------------------

void FastJetFinder::Process(){

  Candidate  *candidate, *constituent;
  TLorentzVector momentum;
  Double_t deta, dphi, detaMax, dphiMax;
  Int_t number;
  Double_t rho = 0;

  std::vector<fastjet::PseudoJet> inputList, outputList;
  std::map<Double_t,Double_t >::iterator itEtaRangeMap;

  DelphesFactory *factory = GetFactory();

  // loop over input objects
  inputList.clear();
  outputList.clear();

  fItInputArray->Reset();
  number = 0;
  while((candidate = static_cast<Candidate*>(fItInputArray->Next()))){
    // for genjets mostly
    if (fKeepPileUp == 0 && candidate->IsPU > 0) {
      continue;
    }

    momentum = candidate->Momentum;
    fastjet::PseudoJet jet (momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
    jet.set_user_index(number);
    inputList.push_back(jet);
    ++number;
  }

  // construct jets
  fastjet::ClusterSequenceArea *sequence = new ClusterSequenceArea(inputList, *fDefinition, *fAreaDefinition);

  // compute rho and store it
  if(fComputeRho && fAreaDefinition){
    for(itEtaRangeMap = fEtaRangeMap.begin(); itEtaRangeMap != fEtaRangeMap.end(); ++itEtaRangeMap){
      Selector select_rapidity = SelectorAbsRapRange(itEtaRangeMap->first, itEtaRangeMap->second); // define an eta region
      fastjet::JetMedianBackgroundEstimator estimator(select_rapidity,*fDefinition,*fAreaDefinition);
      estimator.set_particles(inputList);
      rho = estimator.rho();
      //store rho
      candidate = factory->NewCandidate();
      candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
      candidate->Edges[0] = itEtaRangeMap->first;
      candidate->Edges[1] = itEtaRangeMap->second;
      fRhoOutputArray->Add(candidate);
    }
  }

  else if(fComputeRhoGrid && fAreaDefinition){
    for(itEtaRangeMap = fEtaRangeMap.begin(); itEtaRangeMap != fEtaRangeMap.end(); ++itEtaRangeMap){
      // numbers taken from here https://cmssdt.cern.ch/SDT/lxr/source/RecoJets/JetAlgorithms/src/FixedGridEnergyDensity.cc?v=CMSSW_7_2_0_pre8
      double etaSpacing  = 0;
      if(fabs(itEtaRangeMap->first) < 2.1 and fabs(itEtaRangeMap->second) < 2.1) etaSpacing = 0.6;
      else if(fabs(itEtaRangeMap->first) > 2.1 and fabs(itEtaRangeMap->second) > 2.1) etaSpacing = 0.6;
      else etaSpacing = 0.6;
     
      fastjet::RectangularGrid grid (itEtaRangeMap->first,itEtaRangeMap->second,etaSpacing,TMath::TwoPi()/20.);
      fastjet::GridMedianBackgroundEstimator estimator(grid);
      estimator.set_particles(inputList);
      rho = estimator.rho();
      //store rho
      candidate = factory->NewCandidate();
      candidate->Momentum.SetPtEtaPhiE(rho, 0.0, 0.0, rho);
      candidate->Edges[0] = itEtaRangeMap->first;
      candidate->Edges[1] = itEtaRangeMap->second;
      fRhoOutputArray->Add(candidate);
    }    
  }

  else{

   // take inclusive jets with a pt cut applied  
   outputList.clear();
   outputList = sorted_by_pt(sequence->inclusive_jets(fJetPTMin));

   // loop over all jets and export them
   detaMax = 0.0;
   dphiMax = 0.0;
  
   std::vector<fastjet::PseudoJet>::const_iterator itInputList, itOutputList;
   PseudoJet area;

   // Loop on the outputjets after clustering
   for(itOutputList = outputList.begin(); itOutputList != outputList.end(); ++itOutputList){

    // set momentum
    momentum.SetPxPyPzE(itOutputList->px(), itOutputList->py(), itOutputList->pz(), itOutputList->E());
    area.reset(0.0, 0.0, 0.0, 0.0);
    if(fAreaDefinition) area = itOutputList->area_4vector(); // take the jet aarea

    candidate = factory->NewCandidate();

    inputList.clear();

    // filter away the ghosts
    std::vector<fastjet::PseudoJet> ghosts,jetParticles;
    SelectorIsPureGhost().sift(itOutputList->constituents(), ghosts, jetParticles);

    for(itInputList = jetParticles.begin(); itInputList != jetParticles.end(); ++itInputList){

      // Take the original constistuen from delphes particle array
      constituent = static_cast<Candidate*>(fInputArray->At(itInputList->user_index()));
      deta = TMath::Abs(momentum.Eta()-constituent->Momentum.Eta());
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
    //-----------------------------------
    if (itOutputList->perp() > 200){

      //------------------------------------
      // Trimming
      //------------------------------------

      double Rtrim   = 0.2;
      double ptfrac = 0.05;
      fastjet::Filter    trimmer(fastjet::JetDefinition(fastjet::kt_algorithm,Rtrim),fastjet::SelectorPtFractionMin(ptfrac));
      fastjet::PseudoJet trimmed_jet = trimmer(*itOutputList);

      candidate->TrimmedMass = trimmed_jet.m();
      candidate->TrimmedPt   = trimmed_jet.pt();
      candidate->TrimmedEta  = trimmed_jet.eta();
      candidate->TrimmedPhi  = trimmed_jet.phi();
      
      // three hardest subjet 
      std::vector<fastjet::PseudoJet>  subjets = trimmed_jet.pieces();
      subjets    = sorted_by_pt(subjets);
      candidate->NSubJetsTrimmed = subjets.size();

      for (size_t i = 0; i < subjets.size() and i < 3; i++){
	if(subjets.at(i).pt() < 0) continue ; 
        if(i == 1){ 	
         candidate->TrimmedMassSub1 = subjets.at(i).m();
         candidate->TrimmedPtSub1   = subjets.at(i).pt();
         candidate->TrimmedEtaSub1  = subjets.at(i).eta();
         candidate->TrimmedPhiSub1  = subjets.at(i).phi();
	}
        else if(i == 2){ 	
         candidate->TrimmedMassSub2 = subjets.at(i).m();
         candidate->TrimmedPtSub2   = subjets.at(i).pt();
         candidate->TrimmedEtaSub2  = subjets.at(i).eta();
         candidate->TrimmedPhiSub2  = subjets.at(i).phi();
	}
        else if(i == 3){ 	
         candidate->TrimmedMassSub3 = subjets.at(i).m();
         candidate->TrimmedPtSub3   = subjets.at(i).pt();
         candidate->TrimmedEtaSub3  = subjets.at(i).eta();
         candidate->TrimmedPhiSub3  = subjets.at(i).phi();
	}
      }

      //------------------------------------
      // Pruning
      //------------------------------------

      double Zcut   = 0.1;
      double Rcut   = 0.5;
      double Rprun  = 0.8;

      fastjet::Pruner    pruner(fastjet::JetDefinition(fastjet::cambridge_algorithm,Rprun),Zcut,Rcut);
      fastjet::PseudoJet pruned_jet = pruner(*itOutputList);

      candidate->PrunedMass = pruned_jet.m();
      candidate->PrunedPt   = pruned_jet.pt();
      candidate->PrunedEta  = pruned_jet.eta();
      candidate->PrunedPhi  = pruned_jet.phi();
      subjets.clear();
      
      // three hardest subjet 
      subjets    = pruned_jet.pieces();
      subjets    = sorted_by_pt(subjets);
      candidate->NSubJetsPruned = subjets.size();

      for (size_t i = 0; i < subjets.size() and i < 3; i++){
	if(subjets.at(i).pt() < 0) continue ; 
        if(i == 1){ 	
         candidate->PrunedMassSub1 = subjets.at(i).m();
         candidate->PrunedPtSub1   = subjets.at(i).pt();
         candidate->PrunedEtaSub1  = subjets.at(i).eta();
         candidate->PrunedPhiSub1  = subjets.at(i).phi();
	}
        else if(i == 2){ 	
         candidate->PrunedMassSub2 = subjets.at(i).m();
         candidate->PrunedPtSub2   = subjets.at(i).pt();
         candidate->PrunedEtaSub2  = subjets.at(i).eta();
         candidate->PrunedPhiSub2  = subjets.at(i).phi();
	}
        else if(i == 3){ 	
         candidate->PrunedMassSub3 = subjets.at(i).m();
         candidate->PrunedPtSub3   = subjets.at(i).pt();
         candidate->PrunedEtaSub3  = subjets.at(i).eta();
         candidate->PrunedPhiSub3  = subjets.at(i).phi();
	}
      }

      //------------------------------------
      // SoftDrop
      //------------------------------------

      double beta   = 0.;
      double symmetry_cut   = 0.1;
      double R0     = 0.8;

      contrib::SoftDrop  softDrop(beta,symmetry_cut,R0);
      fastjet::PseudoJet softdrop_jet = softDrop(*itOutputList);

      candidate->SoftDropMass = softdrop_jet.m();
      candidate->SoftDropPt   = softdrop_jet.pt();
      candidate->SoftDropEta  = softdrop_jet.eta();
      candidate->SoftDropPhi  = softdrop_jet.phi();

      // three hardest subjet 
      subjets.clear();
      subjets    = softdrop_jet.pieces();
      subjets    = sorted_by_pt(subjets);
      candidate->NSubJetsSoftDrop = softdrop_jet.pieces().size();

      for (size_t i = 0; i < subjets.size()  and i < 3; i++){
	if(subjets.at(i).pt() < 0) continue ; 
        if(i == 1){ 	
         candidate->SoftDropMassSub1 = subjets.at(i).m();
         candidate->SoftDropPtSub1   = subjets.at(i).pt();
         candidate->SoftDropEtaSub1  = subjets.at(i).eta();
         candidate->SoftDropPhiSub1  = subjets.at(i).phi();
	}
        else if(i == 2){ 	
         candidate->SoftDropMassSub2 = subjets.at(i).m();
         candidate->SoftDropPtSub2   = subjets.at(i).pt();
         candidate->SoftDropEtaSub2  = subjets.at(i).eta();
         candidate->SoftDropPhiSub2  = subjets.at(i).phi();
	}
        else if(i == 3){ 	
         candidate->SoftDropMassSub3 = subjets.at(i).m();
         candidate->SoftDropPtSub3   = subjets.at(i).pt();
         candidate->SoftDropEtaSub3  = subjets.at(i).eta();
         candidate->SoftDropPhiSub3  = subjets.at(i).phi();
	}
      }

      //------------------------------------
      // NSubJettiness
      //------------------------------------
      beta = 1.0;     // power for angular dependence, e.g. beta = 1 --> linear k-means, beta = 2 --> quadratic/classic k-means
      R0   = 0.8;     // Characteristic jet radius for normalization
      Rcut = 10000.0; // maximum R particles can be from axis to be included in jet (large value for no cutoff)
      
      fastjet::contrib::Nsubjettiness nSub1(1,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      fastjet::contrib::Nsubjettiness nSub2(2,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      fastjet::contrib::Nsubjettiness nSub3(3,fastjet::contrib::Njettiness::onepass_kt_axes,beta,R0,Rcut);
      candidate->Tau1 = nSub1(*itOutputList);
      candidate->Tau2 = nSub2(*itOutputList);
      candidate->Tau3 = nSub3(*itOutputList);

    }

    fOutputArray->Add(candidate);
   }
  }

  delete sequence;
}
