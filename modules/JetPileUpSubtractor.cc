/** \class JetPileUpSubtractor
 *  Subtract pile-up contribution from jets using the fastjet area method
 *  $Date: 2012-11-18 15:57:08 +0100 (Sun, 18 Nov 2012) $
 *  $Revision: 814 $
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
*/

#include "modules/JetPileUpSubtractor.h"

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

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/contrib/SafeSubtractor.hh"

#include "fastjet/SISConePlugin.hh"
#include "fastjet/CDFMidPointPlugin.hh"
#include "fastjet/CDFJetCluPlugin.hh"

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"


using namespace std;
using namespace fastjet;
using namespace contrib;

//------------------------------------------------------------------------------

JetPileUpSubtractor::JetPileUpSubtractor() :
  fItJetInputArray(0), fItRhoInputArray(0), fPlugin(0), fDefinition(0), fAreaDefinition(0), fItInputArray(0){

}

//------------------------------------------------------------------------------

JetPileUpSubtractor::~JetPileUpSubtractor(){

}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Init(){

  // import input array(s)
  fJetInputArray   = ImportArray(GetString("JetInputArray","FastJetFinder/jets"));
  fItJetInputArray = fJetInputArray->MakeIterator();

  fRhoInputArray   = ImportArray(GetString("RhoInputArray", "Rho/rho"));
  fItRhoInputArray = fRhoInputArray->MakeIterator();

  fInputArray   = ImportArray(GetString("InputArray","EFlowMerger/eflow"));
  fItInputArray = fInputArray->MakeIterator();

  fSafe4VAreaSubtraction = GetBool("doSafe4VAreaSubtraction", false);
  
  JetDefinition::Plugin *plugin    = NULL;
  JetDefinition::Plugin *pluginRho = NULL;

  fJetPTMin = GetDouble("JetPTMin", 20.0);

  // read eta ranges                                                                                                                                                                       
  ExRootConfParam param = GetParam("RhoEtaRange");
  Long_t i, size;
  fEtaRangeMap.clear();
  size = param.GetSize();
  for(i = 0; i < size/2; ++i) fEtaRangeMap[param[i*2].GetDouble()] = param[i*2 + 1].GetDouble();

  // define algorithm                                                                                                                                                                
  fJetAlgorithm     = GetInt("JetAlgorithm", 6);
  fJetAlgorithmRho  = GetInt("JetAlgorithmRho", 4);
  fParameterR       = GetDouble("ParameterR", 0.5);
  fParameterRRho    = GetDouble("ParameterRRho", 0.4);
  fConeRadius       = GetDouble("ConeRadius", 0.5);
  fConeRadiusRho    = GetDouble("ConeRadiusRho", 0.5);
  fSeedThreshold    = GetDouble("SeedThreshold", 1.0);
  fSeedThresholdRho = GetDouble("SeedThresholdRho", 1.0);
  fConeAreaFraction = GetDouble("ConeAreaFraction", 1.0);
  fConeAreaFractionRho = GetDouble("ConeAreaFractionRho", 1.0);
  fMaxIterations    = GetInt("MaxIterations", 100);
  fMaxIterationsRho = GetInt("MaxIterationsRho", 100);
  fMaxPairSize      = GetInt("MaxPairSize", 2);
  fMaxPairSizeRho   = GetInt("MaxPairSizeRho", 2);
  fIratch           = GetInt("Iratch", 1);
  fIratchRho        = GetInt("IratchRho", 1);
  fAdjacencyCut     = GetDouble("AdjacencyCut", 2.0);
  fAdjacencyCutRho  = GetDouble("AdjacencyCutRho", 2.0);
  fOverlapThresholdRho = GetDouble("OverlapThresholdRho", 0.75);

  // ---  Jet Area Parameters ---                                                                                                                                                        
  fAreaAlgorithm    = GetInt("AreaAlgorithm", 0);
  fAreaAlgorithmRho = GetInt("AreaAlgorithmRho", 0);

  // - ghost based areas -                                                                                                                                                              
  fGhostEtaMax    = GetDouble("GhostEtaMax", 5.0);
  fGhostEtaMaxRho = GetDouble("GhostEtaMaxRho", 5.0);
  fRepeat         = GetInt("Repeat", 1);
  fRepeatRho      = GetInt("RepeatRho", 1);
  fGhostArea      = GetDouble("GhostArea", 0.01);
  fGhostAreaRho   = GetDouble("GhostAreaRho", 0.01);
  fGridScatter    = GetDouble("GridScatter", 1.0);
  fGridScatterRho = GetDouble("GridScatterRho", 1.0);
  fPtScatter      = GetDouble("PtScatter", 0.1);
  fPtScatterRho   = GetDouble("PtScatterRho", 0.1);
  fMeanGhostPt    = GetDouble("MeanGhostPt", 1.0E-100);
  fMeanGhostPtRho = GetDouble("MeanGhostPtRho", 1.0E-100);

  // - voronoi based areas -                                                                                                                                                               
  fEffectiveRfactRho = GetDouble("EffectiveRfactRho", 1.0);
  fEffectiveRfact    = GetDouble("EffectiveRfact", 1.0);

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

  switch(fAreaAlgorithmRho){
  case 1:
    fAreaDefinitionRho = new fastjet::AreaDefinition(active_area_explicit_ghosts,GhostedAreaSpec(fGhostEtaMaxRho,fRepeatRho,fGhostAreaRho,fGridScatterRho,fPtScatterRho,fMeanGhostPtRho));
    break;
  case 2:
    fAreaDefinitionRho = new fastjet::AreaDefinition(one_ghost_passive_area,GhostedAreaSpec(fGhostEtaMaxRho,fRepeatRho,fGhostAreaRho,fGridScatterRho,fPtScatterRho,fMeanGhostPtRho));
    break;
  case 3:
    fAreaDefinitionRho = new fastjet::AreaDefinition(passive_area,GhostedAreaSpec(fGhostEtaMaxRho,fRepeatRho,fGhostAreaRho,fGridScatterRho,fPtScatterRho,fMeanGhostPtRho));
    break;
  case 4:
    fAreaDefinitionRho = new fastjet::AreaDefinition(VoronoiAreaSpec(fEffectiveRfactRho));
    break;
  case 5:
    fAreaDefinitionRho = new fastjet::AreaDefinition(active_area, GhostedAreaSpec(fGhostEtaMax,fRepeat,fGhostArea,fGridScatter,fPtScatter,fMeanGhostPt));
    break;
  default:
  case 0:
    fAreaDefinitionRho = new fastjet::AreaDefinition(active_area_explicit_ghosts,GhostedAreaSpec(fGhostEtaMax,fRepeatRho,fGhostAreaRho,fGridScatterRho,fPtScatterRho,fMeanGhostPtRho));
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

  switch(fJetAlgorithmRho){
  case 1:
    pluginRho      = new fastjet::CDFJetCluPlugin(fSeedThresholdRho,fConeRadiusRho,fAdjacencyCutRho,fMaxIterationsRho,fIratchRho,fOverlapThresholdRho);
    fDefinitionRho = new fastjet::JetDefinition(pluginRho);
    break;
  case 2:
    pluginRho      = new fastjet::CDFMidPointPlugin(fSeedThresholdRho,fConeRadiusRho,fConeAreaFractionRho,fMaxPairSizeRho,fMaxIterationsRho,fOverlapThresholdRho);
    fDefinitionRho = new fastjet::JetDefinition(pluginRho);
    break;
  case 3:
    pluginRho      = new fastjet::SISConePlugin(fConeRadiusRho,fOverlapThresholdRho,fMaxIterationsRho,fJetPTMin);
    fDefinitionRho = new fastjet::JetDefinition(pluginRho);
    break;
  case 4:
    fDefinitionRho = new fastjet::JetDefinition(fastjet::kt_algorithm,fParameterRRho);
    break;
  case 5:
    fDefinitionRho = new fastjet::JetDefinition(fastjet::cambridge_algorithm,fParameterRRho);
    break;
  default:
  case 6:
    fDefinitionRho = new fastjet::JetDefinition(fastjet::antikt_algorithm, fParameterRRho);
    break;
  }

  fPluginRho = pluginRho;
  fPlugin = plugin;
  ClusterSequence::print_banner();

  // create output array(s)
  fOutputArray = ExportArray(GetString("OutputArray", "jets"));

  
}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Finish(){
  if(fItRhoInputArray) delete fItRhoInputArray;
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItInputArray)   delete fItInputArray;
  if(fDefinition)     delete fDefinition;
  if(fAreaDefinition) delete fAreaDefinition;
  if(fDefinitionRho)     delete fDefinitionRho;
  if(fAreaDefinitionRho) delete fAreaDefinitionRho;
  if(fPlugin) delete static_cast<JetDefinition::Plugin*>(fPlugin);
  if(fPluginRho) delete static_cast<JetDefinition::Plugin*>(fPluginRho);
}

//------------------------------------------------------------------------------

void JetPileUpSubtractor::Process(){

  Candidate *candidate, *object;
  TLorentzVector momentum;
  Double_t eta = 0.0;
  Double_t rho = 0.0;

  Double_t deta, dphi, detaMax, dphiMax;
  Int_t number;

  DelphesFactory *factory = GetFactory();

  // loop over all input candidates
  if(not fSafe4VAreaSubtraction){

   TLorentzVector area ;
   fItJetInputArray->Reset();
   while((candidate = static_cast<Candidate*>(fItJetInputArray->Next()))){

    momentum = candidate->Momentum;
    area     = candidate->Area; // take the jet area
    eta      = TMath::Abs(momentum.Eta());
    // find rho
    rho = 0.0;

    if(fRhoInputArray){
      fItRhoInputArray->Reset();
      while((object = static_cast<Candidate*>(fItRhoInputArray->Next()))){
        if(eta >= object->Edges[0] && eta < object->Edges[1]){
          rho = object->Momentum.Pt();
        }
      }
    }

    // apply pile-up correction
    if(momentum.Pt() <= rho * area.Pt()) continue;
    momentum -= rho * area;
    
    if(momentum.Pt() <= fJetPTMin) continue;
    candidate = static_cast<Candidate*>(candidate->Clone());
    candidate->Momentum = momentum;
    fOutputArray->Add(candidate);

   }
  }
  else{

    // calculate the correction 4V safe subtraction
    // loop over input particle objects                                                                                                                                                    
    fItInputArray->Reset();
    number = 0;
    std::map< Double_t, Double_t >::iterator itEtaRangeMap;
    std::vector<PseudoJet> inputList; 
    while((candidate = static_cast<Candidate*>(fItInputArray->Next()))){
      momentum = candidate->Momentum;
      PseudoJet jet (momentum.Px(), momentum.Py(), momentum.Pz(), momentum.E());
      jet.set_user_index(number);
      inputList.push_back(jet);
      ++number;
    }

    // construct background estimators                                                                                                                                                  
    std::vector<fastjet::JetMedianBackgroundEstimator*> bge_rho ;
    std::vector<fastjet::JetMedianBackgroundEstimator*> bge_rhom ;
    std::vector<BackgroundJetPtMDensity*> m_density;
    bge_rho.clear();
    bge_rhom.clear();
    m_density.clear();    
    
    fastjet::ClusterSequenceArea *sequenceRho = new ClusterSequenceArea(inputList, *fDefinitionRho, *fAreaDefinitionRho);

    for(itEtaRangeMap = fEtaRangeMap.begin(); itEtaRangeMap != fEtaRangeMap.end(); ++itEtaRangeMap){
      Selector select_rapidity = SelectorAbsRapRange(itEtaRangeMap->first, itEtaRangeMap->second); // define an eta region                                                                 
      bge_rho.push_back(new fastjet::JetMedianBackgroundEstimator(select_rapidity, *sequenceRho));
      bge_rhom.push_back(new fastjet::JetMedianBackgroundEstimator(select_rapidity,*sequenceRho));
      m_density.push_back(new BackgroundJetPtMDensity);
      bge_rhom.back()->set_jet_density_class(m_density.back());
    }

    fastjet::ClusterSequenceArea *sequence = new ClusterSequenceArea(inputList, *fDefinition, *fAreaDefinition);
    std::vector<fastjet::PseudoJet> jetList = sequence->inclusive_jets(fJetPTMin);

    contrib::SafeAreaSubtractor* area_subtractor = 0;
    std::vector<PseudoJet> subtractedJets ;

    for( size_t iJet = 0; iJet < jetList.size(); iJet++){
      // select the right eta region
      int EtaRegion = 0 ; 
      for(itEtaRangeMap = fEtaRangeMap.begin(); itEtaRangeMap != fEtaRangeMap.end(); ++itEtaRangeMap){
       Selector select_rapidity = SelectorAbsRapRange(itEtaRangeMap->first, itEtaRangeMap->second); // define an eta region                                                                
       if(select_rapidity.pass(jetList.at(iJet))){
	 area_subtractor = new contrib::SafeAreaSubtractor(bge_rho.at(EtaRegion),bge_rhom.at(EtaRegion));
         subtractedJets.push_back((*area_subtractor)(jetList.at(iJet)));
         break;
       }
       ++EtaRegion;
      }
    } 

    // loop over all jets and export them                                                                                                                                          
    detaMax = 0.0;
    dphiMax = 0.0;

    // Loop on the outputjets after clustering                                                                                                                                     
    PseudoJet area ;    
    for(std::vector<PseudoJet>::const_iterator itOutputList = subtractedJets.begin(); itOutputList != subtractedJets.end(); ++itOutputList){

      // set momentum                                                                                               
      if((*itOutputList).pt() < fJetPTMin) continue ;                                                                       
      momentum.SetPxPyPzE(itOutputList->px(), itOutputList->py(), itOutputList->pz(), itOutputList->E());
      area.reset(0.0, 0.0, 0.0, 0.0);
      if(fAreaDefinition) area = itOutputList->area_4vector(); // take the jet aarea                                                                                                       

      candidate = factory->NewCandidate();

      // filter away the ghosts                                                                                                                                                            
      std::vector<fastjet::PseudoJet> ghosts,jetParticles;
      SelectorIsPureGhost().sift(itOutputList->constituents(), ghosts, jetParticles);
 
      for(std::vector<PseudoJet>::const_iterator itInputList = jetParticles.begin(); itInputList != jetParticles.end(); ++itInputList){

	// Take the original constistuen from delphes particle array                                                                                                                       
	object = static_cast<Candidate*>(fInputArray->At(itInputList->user_index()));
	deta = TMath::Abs(momentum.Eta()-object->Momentum.Eta());
	dphi = TMath::Abs(momentum.DeltaPhi(object->Momentum));
	if(deta > detaMax) detaMax = deta;
	if(dphi > dphiMax) dphiMax = dphi;
	candidate->AddCandidate(object);
      }

      candidate->Momentum = momentum;
      candidate->Area.SetPxPyPzE(area.px(), area.py(), area.pz(), area.E());

      candidate->DeltaEta = detaMax;
      candidate->DeltaPhi = dphiMax;

      //------------------------------------                                                                                                                         
      // SubStructure                                                                                                                                                             
      //------------------------------------                                                                                                                                           
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
	// SoftDrop                                                                                                                                                                                //------------------------------------                                                                                                                                                
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

    if(sequence) delete sequence ;
    if(sequenceRho) delete sequenceRho ;
    if(area_subtractor) delete  area_subtractor;
  }
}
