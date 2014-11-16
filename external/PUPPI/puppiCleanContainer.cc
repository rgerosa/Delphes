#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "fastjet/Selector.hh"

#include <algorithm>

#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"

using namespace std;

// ------------- Constructor
puppiCleanContainer::puppiCleanContainer(std::vector<RecoObj> inParticles, 
                                         std::vector<puppiAlgoBin> puppiAlgo,
                                         float minPuppiWeight,
                                         bool  fUseExp){

    // take the input particles
    fRecoParticles_.clear();
    fRecoParticles_ = inParticles;

    // puppi algo 
    puppiAlgo_.clear();
    puppiAlgo_ = puppiAlgo; 

    // min puppi weight
    fMinPuppiWeight_ = minPuppiWeight;

    //Clear everything
    fPFParticles_.clear();
    fPFchsParticles_.clear();
    fChargedPV_.clear();
    fChargedNoPV_.clear();
    fPuppiWeights_.clear();

    fNPV_    = 1 ;
    fPVFrac_ = 0.;
    fUseExp_ = fUseExp;

    //Link to the RecoObjects --> loop on the input particles
    for (unsigned int i = 0; i < fRecoParticles_.size(); i++){
        fastjet::PseudoJet curPseudoJet;
        curPseudoJet.reset_PtYPhiM (fRecoParticles_[i].pt,fRecoParticles_[i].eta,fRecoParticles_[i].phi,fRecoParticles_[i].m);
        curPseudoJet.set_user_index(fRecoParticles_[i].id);  
        // fill vector of pseudojets for internal references
        fPFParticles_.push_back(curPseudoJet);
        if(fRecoParticles_[i].id <= 2) fPFchsParticles_.push_back(curPseudoJet);    //Remove Charged particles associated to other vertex
        if(fRecoParticles_[i].id == 2) fChargedPV_.push_back(curPseudoJet);         //Take Charged particles associated to PV
        if(fRecoParticles_[i].id == 3) fChargedNoPV_.push_back(curPseudoJet);
        if(fRecoParticles_[i].id >= 1) fPVFrac_++ ;
	if(fNPV_ < fRecoParticles_[i].vtxId) fNPV_ = fRecoParticles_[i].vtxId;
    }

    fPVFrac_ = double(fChargedPV_.size())/fPVFrac_;

}

// ------------- De-Constructor
puppiCleanContainer::~puppiCleanContainer(){}

// main function to compute puppi Event
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent(){

    // output particles
    std::vector<fastjet::PseudoJet> particles;
    particles.clear();

    // calculate puppi metric, RMS and mean value for all the algorithms
    for(size_t iPuppiAlgo = 0; iPuppiAlgo < puppiAlgo_.size(); iPuppiAlgo++){
      getRMSAvg(iPuppiAlgo,fPFParticles_,fChargedPV_); // give all the particles in the event and the charged one
    }

    for(size_t iPart = 0; iPart < fPFParticles_.size(); iPart++) {

      double pWeight = 1;
      //Get the Puppi Id and if move on or not .. Pt criteria should be taken into account by itself
      int pPupId = getPuppiId(fPFParticles_[iPart].pt(),fPFParticles_[iPart].eta(),puppiAlgo_);
      if(pPupId == -1) { // out acceptance should be kept as they are
        fPuppiWeights_.push_back(pWeight);
        fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());    
        curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
        particles.push_back(curjet);                                                                                                                                                
	continue;
      }

      if(fPFParticles_.at(iPart).pt() < puppiAlgo_.at(pPupId).fPtMin_){ // low momentum particles
        fPuppiWeights_.push_back(0);
	continue;
      }
   
      // fill the p-values
      double pChi2 = 0;   
      if(fUseExp_){
       //Compute an Experimental Puppi Weight with delta Z info (very simple example)
       if(iPart <= fRecoParticles_.size() and fRecoParticles_[iPart].id == fPFParticles_.at(iPart).user_index()){
        pChi2 = getChi2FromdZ(fRecoParticles_[iPart].dZ);
        //Now make sure Neutrals are not set
        if(fRecoParticles_[iPart].pfType > 3) pChi2 = 0;
       }
      }
      
      puppiParticle partTmp ;
      int found = 0;
      if(fabs(fPFParticles_[iPart].user_index()) <= 2 and puppiAlgo_.at(pPupId).fUseCharged_){
	for(size_t puppiIt = 0 ; puppiIt < puppiAlgo_.at(pPupId).fPuppiParticlesPV_.size(); puppiIt++){
          if(puppiAlgo_.at(pPupId).fPuppiParticlesPV_.at(puppiIt).fPosition_ == int(iPart)){
            partTmp = puppiAlgo_.at(pPupId).fPuppiParticlesPV_.at(puppiIt);           
            found = 1 ;
            break;
	  }
	}
      }
      else if ((fabs(fPFParticles_[iPart].user_index()) <= 2 and !puppiAlgo_.at(pPupId).fUseCharged_) or fabs(fPFParticles_[iPart].user_index()) >= 3){
	for(size_t puppiIt = 0; puppiIt < puppiAlgo_.at(pPupId).fPuppiParticlesPU_.size(); puppiIt++){
          if(puppiAlgo_.at(pPupId).fPuppiParticlesPU_.at(puppiIt).fPosition_ == int(iPart)){
            partTmp = puppiAlgo_.at(pPupId).fPuppiParticlesPU_.at(puppiIt);
            found = 1;
            break;
	  }
	}
      }

      // means that is inside the NULL vector for some reasons
      if(found == 0){
	for(size_t puppiIt = 0; puppiIt < puppiAlgo_.at(pPupId).fPuppiParticlesNULL_.size(); puppiIt++){
          if(puppiAlgo_.at(pPupId).fPuppiParticlesNULL_.at(puppiIt).fPosition_ == int(iPart)){
            partTmp = puppiAlgo_.at(pPupId).fPuppiParticlesNULL_.at(puppiIt);
            found = 1 ;
            break;
	  }
	}
      }

      if(found == 0){
	fPuppiWeights_.push_back(pWeight);
        fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());    
        curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
        particles.push_back(curjet);                                                                                                                                                
        continue;
      }
      else if(found == 1 and partTmp.fPval_ == -999){
	fPuppiWeights_.push_back(pWeight);
        fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());    
        curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                      
        particles.push_back(curjet);                                                                                                                                                
        continue;
      }

      
      pWeight = compute(partTmp.fPval_,pChi2,puppiAlgo_.at(pPupId));      
      
      if(fPFParticles_[iPart].user_index() == 1 && puppiAlgo_.at(pPupId).fApplyCHS_ ) pWeight = 1;
      if(fPFParticles_[iPart].user_index() == 2 && puppiAlgo_.at(pPupId).fApplyCHS_ ) pWeight = 0;

      //Basic Weight Checks
      if( std::isinf(pWeight) || std::isnan(pWeight)){
	std::cerr << "====> Weight is nan : pt " << fPFParticles_[iPart].pt() << " -- eta : " << fPFParticles_[iPart].eta() << " -- Value" << partTmp.fPval_ << " -- id : " << fPFParticles_[iPart].user_index() << std::endl;
         pWeight = 1;
      }

      //Basic Cuts
      if(pWeight < fMinPuppiWeight_) pWeight = 0; //==> Elminate the low Weight stuff
      
      //threshold cut on the neutral Pt
      if(pWeight*fPFParticles_[iPart].pt() < getNeutralPtCut(puppiAlgo_.at(pPupId).fNeutralMinE_,puppiAlgo_.at(pPupId).fNeutralPtSlope_,fNPV_) && fPFParticles_[iPart].user_index() == 1 ) 
       pWeight = 0; 

      fPuppiWeights_.push_back(pWeight);

      //Now get rid of the thrown out weights for the particle collection
      if(pWeight == 0) continue;

      //Produce
      fastjet::PseudoJet curjet( pWeight*fPFParticles_[iPart].px(), pWeight*fPFParticles_[iPart].py(), pWeight*fPFParticles_[iPart].pz(), pWeight*fPFParticles_[iPart].e());           
      curjet.set_user_index(fPFParticles_[iPart].user_index());                                                                                                                            
      particles.push_back(curjet);                                                                                                                                                
    }
    
    return particles;
}

void puppiCleanContainer::getRMSAvg(int iPuppiAlgo, std::vector<fastjet::PseudoJet> & particlesAll, std::vector<fastjet::PseudoJet> &chargedPV) { 

  std::vector<puppiParticle> puppiParticles; // puppi particles to be set for a specific algo
  puppiParticles.clear();

  // Loop on all the particles of the event  
  for(size_t iPart = 0; iPart < particlesAll.size(); iPart++ ) { 

    float pVal    = -999;
    int   pPupId  = getPuppiId(particlesAll[iPart].pt(),particlesAll[iPart].eta(),puppiAlgo_); // get the puppi id algo asaf of eta and phi of the particle

    // does not exsist and algorithm for this particle, store -1    
    if(pPupId == -1){ 
      puppiParticles.push_back(puppiParticle(particlesAll.at(iPart).pt(),particlesAll.at(iPart).eta(),pVal,particlesAll.at(iPart).user_index(),iPart));
      continue;
    }

    if(pPupId != iPuppiAlgo) continue; // mis-match of the index, it will call later.
    
    // apply CHS -> use only LV hadrons to compute the metric for each particle
    if(puppiAlgo_.at(pPupId).fUseCharged_)  
       pVal = goodVar(particlesAll[iPart], chargedPV,    puppiAlgo_.at(pPupId).fMetricId_,puppiAlgo_.at(pPupId).fConeSize_);
    else if(!puppiAlgo_.at(pPupId).fUseCharged_) 
       pVal = goodVar(particlesAll[iPart], particlesAll, puppiAlgo_.at(pPupId).fMetricId_,puppiAlgo_.at(pPupId).fConeSize_);
    
    // fill the value
    if(std::isnan(pVal) || std::isinf(pVal)) std::cout << "====>  Value is Nan " << pVal << " == " << particlesAll[iPart].pt() << " -- " << particlesAll[iPart].eta() << std::endl;
    if(std::isnan(pVal) || std::isinf(pVal)) continue;
    
    puppiParticles.push_back(puppiParticle(particlesAll.at(iPart).pt(),particlesAll.at(iPart).eta(),pVal,particlesAll.at(iPart).user_index(),iPart));
  }

  // set the puppi particles for the algorithm
  puppiAlgo_.at(iPuppiAlgo).setPuppiParticles(puppiParticles);

  // compute RMS, median and mean value  
  computeMedRMS(iPuppiAlgo);
  
}


float puppiCleanContainer::goodVar(const fastjet::PseudoJet & particle, const std::vector<fastjet::PseudoJet> & particleAll, const int & pPupId, const float & coneSize) {
  float lPup = 0;
  lPup = var_within_R(pPupId,particleAll,particle,coneSize);
  return lPup;
}


// ----------------------
float puppiCleanContainer::compute(const float & val, const float & chi2, const puppiAlgoBin & puppiAlgo) {

  if(puppiAlgo.fMetricId_ == -1) return 1;
  if(puppiAlgo.fPuppiParticlesPU_.size() + puppiAlgo.fPuppiParticlesPV_.size()  == 0) return 1.; 

  float pVal  = val;
  float lVal  = 0.;
  float lPVal = 1.;


  if(puppiAlgo.fMetricId_ == 0 && val == 0) pVal = puppiAlgo.fMedian_;
  if(puppiAlgo.fMetricId_ == 3 && val == 0) pVal = puppiAlgo.fMedian_;
  if(puppiAlgo.fMetricId_ == 5 && val == 0) pVal = puppiAlgo.fMedian_;

  int lNDOF = 0 ;

  lVal += (pVal-puppiAlgo.fMedian_)*(fabs(pVal-puppiAlgo.fMedian_))/puppiAlgo.fRMS_/puppiAlgo.fRMS_;
  lNDOF++;
  if(chi2 != 0) lNDOF++; 
  if(chi2 != 0) lVal+=chi2; //Add external Chi2 to first element

  lPVal *= ROOT::Math::chisquared_cdf(lVal,lNDOF);

  return lPVal;

}

float puppiCleanContainer::getNeutralPtCut(const float & fNeutralMinE, const float & fNeutralPtSlope, const int & fNPV) {
  return fNeutralMinE + fNPV * fNeutralPtSlope;
}

// take the type of algorithm in the vector index
int puppiCleanContainer::getPuppiId(const float & pt, const float & eta, const std::vector<puppiAlgoBin> & puppiAlgos){
  int PuppiId = -1;
  for(size_t iPuppiAlgo = 0; iPuppiAlgo < puppiAlgos.size() ; iPuppiAlgo++){
    if(fabs(eta) < puppiAlgos[iPuppiAlgo].fEtaMin_) continue;
    if(fabs(eta) > puppiAlgos[iPuppiAlgo].fEtaMax_) continue;
    PuppiId = iPuppiAlgo ;
    break;
  }
  return PuppiId;  
}


float puppiCleanContainer::var_within_R(const int & pPupId, const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R){

  if(pPupId == -1) return 1;
  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  float var = 0;

  for(size_t iPart = 0; iPart < near_particles.size(); iPart++){

    double pDEta = near_particles[iPart].eta()-centre.eta();
    double pDPhi = fabs(near_particles[iPart].phi()-centre.phi());
    if(pDPhi > 2.*3.14159265-pDPhi) pDPhi = 2.*3.14159265-pDPhi;
    double pDR = sqrt(pDEta*pDEta+pDPhi*pDPhi);

    if(pDR < 0.001) continue;
    if(pDR < 0.01)  continue;//pDR = 0.01;
    if(pDR == 0)    continue;

    if(pPupId == 0) var += (near_particles[iPart].pt()/pDR/pDR);
    if(pPupId == 1) var += near_particles[iPart].pt();
    if(pPupId == 2) var += (1./pDR)*(1./pDR);
    if(pPupId == 3) var += (1./pDR)*(1./pDR);
    if(pPupId == 4) var += near_particles[iPart].pt();
    if(pPupId == 5) var += (near_particles[iPart].pt()/pDR)*(near_particles[iPart].pt()/pDR);
  }

  if(pPupId == 0 && var != 0) var = log(var);
  if(pPupId == 3 && var != 0) var = log(var);
  if(pPupId == 5 && var != 0) var = log(var);


  return var;
 
}


float puppiCleanContainer::pt_within_R(const std::vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet & centre, const float & R){

  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  double answer = 0.0;
  for(size_t iPart = 0; iPart<near_particles.size(); iPart++){
    answer += near_particles[iPart].pt();
  }

  return answer;
}

fastjet::PseudoJet puppiCleanContainer::flow_within_R(const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R){

  fastjet::Selector sel = fastjet::SelectorCircle(R);
  sel.set_reference(centre);
  std::vector<fastjet::PseudoJet> near_particles = sel(particles);
  fastjet::PseudoJet flow;
  for(unsigned int i=0; i<near_particles.size(); i++){
    flow += near_particles[i];
  }
  return flow;

}


void puppiCleanContainer::computeMedRMS(const int & puppiAlgo) {

  if(puppiAlgo > int(puppiAlgo_.size())  ) return;
  if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size() == 0) return;

  // sort in pt increasing order
  std::sort(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.begin(),puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.end(),puppiValSort());

  // if apply correction
  float lCorr = 1.;
  if(puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_) lCorr *= 1.-fPVFrac_;

  // count the position of the last particle with pval zero coming from PU
  int lNum0 = 0;
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size(); i0++) {
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_[i0].fPval_ == 0) lNum0 = i0;
  }

  // take the median value on PU particles
  int lNHalfway = lNum0 + int(float(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size()-lNum0)*0.50*lCorr);
  puppiAlgo_.at(puppiAlgo).fMedian_ = puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(lNHalfway).fPval_;
  float lMed = puppiAlgo_.at(puppiAlgo).fMedian_; //Just to make the readability easier

  // take the RMS
  int lNRMS = 0;
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size(); i0++) {
    puppiAlgo_.at(puppiAlgo).fMean_ += puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_;
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_ == 0) continue;
    if(!puppiAlgo_.at(puppiAlgo).fUseCharged_ && puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_ && puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_ > lMed) continue;
    lNRMS++;
    puppiAlgo_.at(puppiAlgo).fRMS_ += (puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_-lMed)*( puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.at(i0).fPval_-lMed);
  }

  puppiAlgo_.at(puppiAlgo).fMean_ /= puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size();
  if(lNRMS > 0) puppiAlgo_.at(puppiAlgo).fRMS_/=lNRMS;
  if(puppiAlgo_.at(puppiAlgo).fRMS_ == 0) puppiAlgo_.at(puppiAlgo).fRMS_ = 1e-5;
  puppiAlgo_.at(puppiAlgo).fRMS_  = sqrt(puppiAlgo_.at(puppiAlgo).fRMS_);
  puppiAlgo_.at(puppiAlgo).fRMS_ *= puppiAlgo_.at(puppiAlgo).fRMSScaleFactor_;

  if(!puppiAlgo_.at(puppiAlgo).fApplyLowPUCorr_) return;

  //Adjust the p-value to correspond to the median
  std::sort(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.begin(),puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.end(),puppiValSort());
  int lNPV = 0; 
  for(size_t i0 = 0; i0 < puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.size(); i0++){ 
    if(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_[i0].fPval_ <= lMed ) lNPV++;
  }

  // it helps in puppi the median value close to the mean one when a lot of pval 0 are present
  float lAdjust = 1.5*float(lNPV)/float(puppiAlgo_.at(puppiAlgo).fPuppiParticlesPV_.size()+puppiAlgo_.at(puppiAlgo).fPuppiParticlesPU_.size());
  if(lAdjust > 0) puppiAlgo_.at(puppiAlgo).fMedian_ -= sqrt(ROOT::Math::chisquared_quantile(lAdjust,1.)*puppiAlgo_.at(puppiAlgo).fRMS_);

}


float puppiCleanContainer::getChi2FromdZ(float iDZ) {
   //We need to obtain prob of PU + (1-Prob of LV)
   // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm (its really more like 1mm)
   //double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.2)*2.; //*2 is to do it double sided
   //Take iDZ to be corrected by sigma already
   double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),1.)*2.; //*2 is to do it double sided
   double lProbPU = 1-lProbLV;
   if(lProbPU <= 0) lProbPU = 1e-16; //Quick Trick to through out infs
   if(lProbPU >= 0) lProbPU = 1-1e-16; //Ditto
   double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
   lChi2PU*=lChi2PU;
   return lChi2PU;
}
