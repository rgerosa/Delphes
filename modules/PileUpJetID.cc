#include "modules/PileUpJetID.h"

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
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <assert.h>   

using namespace std;

//------------------------------------------------------------------------------

PileUpJetID::PileUpJetID() :
  fItJetInputArray(0),fTrackInputArray(0),fNeutralInputArray(0){}

//------------------------------------------------------------------------------

PileUpJetID::~PileUpJetID(){}

//------------------------------------------------------------------------------

void PileUpJetID::Init(){

  fJetPTMin        = GetDouble("JetPTMin",     10.0);
  fParameterR      = GetDouble("ParameterR",   0.5);
  fUseConstituents = GetInt("UseConstituents", 0);

  // import input array(s)
  fJetInputArray     = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
  fTrackInputArray   = ImportArray(GetString("TrackInputArray", "ParticlePropagator/tracks"));
  fNeutralInputArray = ImportArray(GetString("NeutralInputArray", "ParticlePropagator/tracks"));
  fPVInputArray      = ImportArray(GetString("PVInputArray", "ModifyBeamSpot/PV"));

  fItJetInputArray     = fJetInputArray->MakeIterator();
  fItTrackInputArray   = fTrackInputArray->MakeIterator();
  fItNeutralInputArray = fNeutralInputArray->MakeIterator();
  fItPVInputArray      = fPVInputArray->MakeIterator();

  // read eta min ranges                                                                                                                                                                   
  ExRootConfParam param = GetParam("Cones");
  fCones.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) fCones.push_back(param[iMap].GetFloat());

  // create output array(s)
  fOutputArray           = ExportArray(GetString("OutputArray", "jets"));

  // Tight Id
  param = GetParam("Pt010_Tight_betaStar");
  Pt010_Tight_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Tight_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Tight_betaStar");
  Pt1020_Tight_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Tight_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Tight_betaStar");
  Pt2030_Tight_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Tight_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Tight_betaStar");
  Pt3050_Tight_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Tight_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt010_Tight_RMS");
  Pt010_Tight_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Tight_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Tight_RMS");
  Pt1020_Tight_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Tight_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Tight_RMS");
  Pt2030_Tight_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Tight_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Tight_RMS");
  Pt3050_Tight_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Tight_RMS.push_back(param[iMap].GetFloat());

  // Medium Id
  param = GetParam("Pt010_Medium_betaStar");
  Pt010_Medium_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Medium_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Medium_betaStar");
  Pt1020_Medium_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Medium_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Medium_betaStar");
  Pt2030_Medium_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Medium_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Medium_betaStar");
  Pt3050_Medium_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Medium_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt010_Medium_RMS");
  Pt010_Medium_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Medium_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Medium_RMS");
  Pt1020_Medium_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Medium_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Medium_RMS");
  Pt2030_Medium_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Medium_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Medium_RMS");
  Pt3050_Medium_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Medium_RMS.push_back(param[iMap].GetFloat());

  // Loose Id
  param = GetParam("Pt010_Loose_betaStar");
  Pt010_Loose_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Loose_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Loose_betaStar");
  Pt1020_Loose_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Loose_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Loose_betaStar");
  Pt2030_Loose_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Loose_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Loose_betaStar");
  Pt3050_Loose_betaStar.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Loose_betaStar.push_back(param[iMap].GetFloat());

  param = GetParam("Pt010_Loose_RMS");
  Pt010_Loose_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt010_Loose_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt1020_Loose_RMS");
  Pt1020_Loose_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt1020_Loose_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt2030_Loose_RMS");
  Pt2030_Loose_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt2030_Loose_RMS.push_back(param[iMap].GetFloat());

  param = GetParam("Pt3050_Loose_RMS");
  Pt3050_Loose_RMS.clear();
  for(int iMap = 0; iMap < param.GetSize(); ++iMap) Pt3050_Loose_RMS.push_back(param[iMap].GetFloat());

  //-------------------------
  for(int i0 = 0; i0 < 3; i0++) {  // number of WP
    if(i0 == 0){

     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][0][i2] = Pt010_Tight_betaStar[i2];     
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][1][i2] = Pt1020_Tight_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][2][i2] = Pt2030_Tight_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][3][i2] = Pt3050_Tight_betaStar[i2];

     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][0][i2] = Pt010_Tight_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][1][i2] = Pt1020_Tight_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][2][i2] = Pt2030_Tight_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][3][i2] = Pt3050_Tight_RMS[i2];
    }
    else if(i0 == 1){
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][0][i2] = Pt010_Medium_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][1][i2] = Pt1020_Medium_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][2][i2] = Pt2030_Medium_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][3][i2] = Pt3050_Medium_betaStar[i2];

     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][0][i2] = Pt010_Medium_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][1][i2] = Pt1020_Medium_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][2][i2] = Pt2030_Medium_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][3][i2] = Pt3050_Medium_RMS[i2];
    }
    else if(i0 == 2){
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][0][i2] = Pt010_Loose_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][1][i2] = Pt1020_Loose_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][2][i2] = Pt2030_Loose_betaStar[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_betaStar[i0][3][i2] = Pt3050_Loose_betaStar[i2];

     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][0][i2] = Pt010_Loose_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][1][i2] = Pt1020_Loose_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][2][i2] = Pt2030_Loose_RMS[i2];
     for(int i2 = 0; i2 < 4; i2++) pileUpIDCut_RMS[i0][3][i2] = Pt3050_Loose_RMS[i2];
    }
  }
}

//------------------------------------------------------------------------------

void PileUpJetID::Finish(){
  if(fItJetInputArray) delete fItJetInputArray;
  if(fItTrackInputArray) delete fItTrackInputArray;
  if(fItNeutralInputArray) delete fItNeutralInputArray;
  if(fItPVInputArray) delete fItPVInputArray;
  fCones.clear();
}

//------------------------------------------------------------------------------

void PileUpJetID::Process(){

  Candidate *jetCandidate, *constituent;
  TLorentzVector momentum, area;

  // loop over all input candidates
  fItJetInputArray->Reset();
  while((jetCandidate = static_cast<Candidate*>(fItJetInputArray->Next()))){

    momentum = jetCandidate->Momentum;
    area     = jetCandidate->Area;

    jetCandidate->FracPt.clear();
    jetCandidate->emFracPt.clear();
    jetCandidate->neutFracPt.clear();
    jetCandidate->chFracPt.clear();

    jetCandidate->FracPt.resize(fCones.size());
    jetCandidate->emFracPt.resize(fCones.size());
    jetCandidate->neutFracPt.resize(fCones.size());
    jetCandidate->chFracPt.resize(fCones.size());

    jetCandidate->dRMean  = 0.;
    jetCandidate->dR2Mean = 0.;
    jetCandidate->ptD     = 0.;
    jetCandidate->sumPt   = 0.;
    jetCandidate->sumPt2  = 0.;
    jetCandidate->dRMeanEm = 0.;
    jetCandidate->ptDNe    = 0.;
    jetCandidate->sumPtNe  = 0.;
    jetCandidate->nNeutral = 0.;
    jetCandidate->neuEMfrac  = 0.;
    jetCandidate->dRMeanNeut = 0.;
    jetCandidate->neuHadfrac = 0.;
    jetCandidate->dRMeanCh = 0.;
    jetCandidate->ptDCh = 0.;
    jetCandidate->sumPtCh = 0.;
    jetCandidate->nCharged = 0;
    jetCandidate->chgEMfrac = 0.;
    jetCandidate->chgHadfrac = 0;
    jetCandidate->betaClassic = 0;
    jetCandidate->betaClassicStar = 0;
    jetCandidate->beta = 0;
    jetCandidate->betaStar = 0;
    jetCandidate->constituents = 0;
    jetCandidate->dZ = 0;
    jetCandidate->d0 = 0;
    jetCandidate->etaW = 0;
    jetCandidate->phiW = 0;
    jetCandidate->jetW = 0;
    jetCandidate->majW = 0;
    jetCandidate->minW = 0;
    jetCandidate->dRLeadCent = 0;
    jetCandidate->dRLead2nd  = 0;
    jetCandidate->ptMean = 0;
    jetCandidate->ptRMS  = 0;
    jetCandidate->pt2A   = 0;
    jetCandidate->sumChPt = 0;
    jetCandidate->sumNePt = 0;
    jetCandidate->axis2   = 0;
    jetCandidate->leadFrac = 0;
    jetCandidate->secondFrac = 0;
    jetCandidate->thirdFrac  = 0;
    jetCandidate->fourthFrac = 0;
    jetCandidate->leadChFrac = 0;
    jetCandidate->secondChFrac = 0;
    jetCandidate->thirdChFrac = 0;
    jetCandidate->fourthChFrac = 0;
    jetCandidate->leadEmFrac = 0;
    jetCandidate->secondEmFrac = 0;
    jetCandidate->thirdEmFrac = 0;
    jetCandidate->fourthEmFrac = 0;
    jetCandidate->leadNeutFrac = 0;
    jetCandidate->secondNeutFrac = 0;
    jetCandidate->thirdNeutFrac = 0; 
    jetCandidate->fourthNeutFrac = 0;
    jetCandidate->pileupIDFlagCutBased = 0;
 
    TMatrixDSym covMatrix(2); covMatrix = 0.;
 
    // Loop on the constituent of each jet
    Candidate leadCand, trailCand, secondCand;
    Candidate leadCandEM;
    Candidate leadCandNeutral;
    Candidate leadCandCh;

    std::vector<float> frac, fracCh, fracEm, fracNeut;

    float sumTkPt = 0.;

    float sum_deta = 0.;
    float sum_dphi  = 0.;

    if (fUseConstituents) {
      TIter itConstituents(jetCandidate->GetCandidates());      
      while((constituent = static_cast<Candidate*>(itConstituents.Next()))) {
        float candPt     = constituent->Momentum.Pt();
        float candDr     = jetCandidate->Momentum.DeltaR(constituent->Momentum);
	float candPtFrac = candPt/jetCandidate->Momentum.Pt();
	float candDeta   = fabs(jetCandidate->Momentum.Eta()-constituent->Momentum.Eta());
	float candDphi   = jetCandidate->Momentum.DeltaPhi(constituent->Momentum);
	float candPtDr   = candPt * candDr;

	if(candPt > leadCand.Momentum.Pt()) {
	  secondCand = leadCand;
	  leadCand  = *constituent; 
	} 
        else if(candPt > secondCand.Momentum.Pt() && candPt < leadCand.Momentum.Pt()) {
	  secondCand = *constituent;
	}

        jetCandidate->dRMean  += candPtDr;
        jetCandidate->dR2Mean += candPtDr*candPtDr;
        covMatrix(0,0)     += candPt*candPt*candDeta*candDeta;
        covMatrix(0,1)     += candPt*candPt*candDeta*candDphi;
        covMatrix(1,1)     += candPt*candPt*candDphi*candDphi;
        jetCandidate->ptD     += candPt*candPt;
        jetCandidate->sumPt   += candPt;
        jetCandidate->sumPt2  += candPt*candPt;

	sum_deta   = candPt*candPt*candDeta;
	sum_dphi   = candPt*candPt*candDphi;
       
	size_t iCone = std::lower_bound(fCones.begin(),fCones.end(),candDr)-fCones.begin();        
	frac.push_back(candPtFrac);

	if( iCone < fCones.size()) jetCandidate->FracPt[iCone] += candPt;
	
        if( TMath::Abs(constituent->PID) == 22){ // look for gamma
	  if(candPt > leadCandEM.Momentum.Pt()) leadCandEM = *constituent; 
          if(iCone  < fCones.size() ) jetCandidate->emFracPt[iCone] += candPt; 
          fracEm.push_back(candPtFrac);
          jetCandidate->dRMeanEm  += candPtDr;
          jetCandidate->ptDNe     += candPt*candPt;
          jetCandidate->sumPtNe   += candPt;
          jetCandidate->neuEMfrac += constituent->Momentum.E();
          jetCandidate->nNeutral ++;
	}

        else if(constituent->Charge == 0 and TMath::Abs(constituent->PID) > 18 and TMath::Abs(constituent->PID) !=21 and TMath::Abs(constituent->PID) !=22 and 
                TMath::Abs(constituent->PID) !=25){ // look for neutral hadrons
 	          if(candPt > leadCandNeutral.Momentum.Pt()) leadCandNeutral = *constituent; 
                  if(iCone < fCones.size()) jetCandidate->neutFracPt[iCone] += candPt;
                  fracNeut.push_back(candPtFrac);
                  jetCandidate->dRMeanNeut += candPtDr;
                  jetCandidate->ptDNe      += candPt*candPt;
                  jetCandidate->sumPtNe    += candPt;
                  jetCandidate->neuHadfrac += constituent->Momentum.E();
                  jetCandidate->nNeutral ++;
	}

        else if(constituent->Charge != 0){ // look for charged particles
 	     if(candPt > leadCandCh.Momentum.Pt()) leadCandCh = *constituent; 
             if(iCone  < fCones.size()) jetCandidate->chFracPt[iCone] += candPt;
             fracCh.push_back(candPtFrac);

             jetCandidate->dRMeanCh  += candPtDr;
             jetCandidate->ptDCh     += candPt*candPt;
             jetCandidate->sumPtCh   += candPt;
             jetCandidate->nCharged ++;

             float tkpt = constituent->Momentum.Pt(); 
   	     sumTkPt += tkpt;

             if(TMath::Abs(constituent->PID) >= 11 and TMath::Abs(constituent->PID) < 18) jetCandidate->chgEMfrac += constituent->Momentum.E(); // EM charged objects
             else jetCandidate->chgHadfrac += constituent->Momentum.E();

             if(constituent->IsRecoPU == 0)   jetCandidate->betaClassic += tkpt;
             if(constituent->IsRecoPU != 0 )  jetCandidate->betaClassicStar += tkpt;
	     if(fabs(constituent->IsPU == 0)) jetCandidate->beta += tkpt;
	     if(fabs(constituent->IsPU != 0)) jetCandidate->betaStar += tkpt;

	}
      	
	// trailing candidate
	if(candPt < trailCand.Momentum.Pt()) {
	  trailCand = *constituent; 
	}
	
      }
    }
    else {
      // Not using constituents, using dr
      fItTrackInputArray->Reset();
      while ((constituent = static_cast<Candidate*>(fItTrackInputArray->Next()))) { // loop on charged particles 
        if (constituent->Momentum.DeltaR(jetCandidate->Momentum) < fParameterR) { // check for particles mathced to the jet

         float candPt     = constituent->Momentum.Pt();
         float candDr     = jetCandidate->Momentum.DeltaR(constituent->Momentum);
  	 float candPtFrac = candPt/jetCandidate->Momentum.Pt();
	 float candDeta   = fabs(jetCandidate->Momentum.Eta()-constituent->Momentum.Eta());
	 float candDphi   = jetCandidate->Momentum.DeltaPhi(constituent->Momentum);
	 float candPtDr   = candPt * candDr;

	 if(candPt > leadCand.Momentum.Pt()) {
	  secondCand = leadCand;
	  leadCand  = *constituent; 
	 } 
         else if(candPt > secondCand.Momentum.Pt() && candPt < leadCand.Momentum.Pt()) {
	  secondCand = *constituent;
	 }

         jetCandidate->dRMean  += candPtDr;
         jetCandidate->dR2Mean += candPtDr*candPtDr;
         covMatrix(0,0)     += candPt*candPt*candDeta*candDeta;
         covMatrix(0,1)     += candPt*candPt*candDeta*candDphi;
         covMatrix(1,1)     += candPt*candPt*candDphi*candDphi;
         jetCandidate->ptD     += candPt*candPt;
         jetCandidate->sumPt   += candPt;
         jetCandidate->sumPt2  += candPt*candPt;

	 sum_deta   = candPt*candPt*candDeta;
  	 sum_dphi   = candPt*candPt*candDphi;

	 size_t iCone = std::lower_bound(fCones.begin(),fCones.end(),candDr)-fCones.begin();        
	 frac.push_back(candPtFrac);
	 if( iCone < fCones.size()) jetCandidate->FracPt[iCone] += candPt;

         if( TMath::Abs(constituent->PID) == 22) continue; // should be a problem
         if(constituent->Charge == 0) continue;

 	 if(candPt > leadCandCh.Momentum.Pt()) leadCandCh = *constituent; 
         if(iCone  < fCones.size()) jetCandidate->chFracPt[iCone] += candPt;
         fracCh.push_back(candPtFrac);

         jetCandidate->dRMeanCh  += candPtDr;
         jetCandidate->ptDCh     += candPt*candPt;
         jetCandidate->sumPtCh   += candPt;
         jetCandidate->nCharged ++;

         float tkpt = constituent->Momentum.Pt(); 
   	 sumTkPt += tkpt;

         if(TMath::Abs(constituent->PID) >= 11 and TMath::Abs(constituent->PID) < 18) jetCandidate->chgEMfrac += constituent->Momentum.E(); // EM charged objects
         else jetCandidate->chgHadfrac += constituent->Momentum.E();

         if(constituent->IsRecoPU == 0)   jetCandidate->betaClassic += tkpt;
         if(constituent->IsRecoPU != 0 )  jetCandidate->betaClassicStar += tkpt;
	 if(fabs(constituent->IsPU == 0)) jetCandidate->beta += tkpt;
	 if(fabs(constituent->IsPU != 0)) jetCandidate->betaStar += tkpt;	 
	}
	// trailing candidate
	if(constituent->Momentum.Pt() < trailCand.Momentum.Pt()) {
	  trailCand = *constituent; 
	}      
      }
    
      fItNeutralInputArray->Reset();
      while ((constituent = static_cast<Candidate*>(fItNeutralInputArray->Next()))) { // Loop on neutral PF jetCandidates
	if (constituent->Momentum.DeltaR(jetCandidate->Momentum) < fParameterR) { // check within the jet

         float candPt     = constituent->Momentum.Pt();
         float candDr     = jetCandidate->Momentum.DeltaR(constituent->Momentum);
  	 float candPtFrac = candPt/jetCandidate->Momentum.Pt();
	 float candDeta   = fabs(jetCandidate->Momentum.Eta()-constituent->Momentum.Eta());
	 float candDphi   = jetCandidate->Momentum.DeltaPhi(constituent->Momentum);
	 float candPtDr   = candPt * candDr;

	 if(candPt > leadCand.Momentum.Pt()) {
	  secondCand = leadCand;
	  leadCand  = *constituent; 
	 } 
         else if(candPt > secondCand.Momentum.Pt() && candPt < leadCand.Momentum.Pt()) {
	  secondCand = *constituent;
	 }

         jetCandidate->dRMean  += candPtDr;
         jetCandidate->dR2Mean += candPtDr*candPtDr;
         covMatrix(0,0)     += candPt*candPt*candDeta*candDeta;
         covMatrix(0,1)     += candPt*candPt*candDeta*candDphi;
         covMatrix(1,1)     += candPt*candPt*candDphi*candDphi;
         jetCandidate->ptD     += candPt*candPt;
         jetCandidate->sumPt   += candPt;
         jetCandidate->sumPt2  += candPt*candPt;

 	 sum_deta   = candPt*candPt*candDeta;
	 sum_dphi   = candPt*candPt*candDphi;

	 size_t iCone = std::lower_bound(fCones.begin(),fCones.end(),candDr)-fCones.begin();        
	 frac.push_back(candPtFrac);
	 if( iCone < fCones.size()) jetCandidate->FracPt[iCone] += candPt;

         if(constituent->Charge != 0 ) continue; // in principle should be a problem

         if( TMath::Abs(constituent->PID) == 22){ // look for gamma
	  if(candPt > leadCandEM.Momentum.Pt()) leadCandEM = *constituent; 
          if(iCone  < fCones.size() ) jetCandidate->emFracPt[iCone] += candPt; 
          fracEm.push_back(candPtFrac);
          jetCandidate->dRMeanEm  += candPtDr;
          jetCandidate->ptDNe     += candPt*candPt;
          jetCandidate->sumPtNe   += candPt;
          jetCandidate->neuEMfrac += constituent->Momentum.E();
          jetCandidate->nNeutral ++;
	 }

         else if(constituent->Charge == 0 and TMath::Abs(constituent->PID) > 18 and TMath::Abs(constituent->PID) !=21 and TMath::Abs(constituent->PID) !=22 and 
                 TMath::Abs(constituent->PID) !=25){ // look for neutral hadrons
 	          if(candPt > leadCandNeutral.Momentum.Pt()) leadCandNeutral = *constituent; 
                  if(iCone < fCones.size()) jetCandidate->neutFracPt[iCone] += candPt;
                  fracNeut.push_back(candPtFrac);
                  jetCandidate->dRMeanNeut += candPtDr;
                  jetCandidate->ptDNe      += candPt*candPt;
                  jetCandidate->sumPtNe    += candPt;
                  jetCandidate->neuHadfrac += constituent->Momentum.E();
                  jetCandidate->nNeutral ++;
	 }
	 // trailing candidate
	 if(candPt < trailCand.Momentum.Pt()) {
	  trailCand = *constituent; 
	 }
	}
      }
    }

    // fix the final values
    assert(leadCand.Momentum.Pt() != 0);
    
    if (secondCand.Momentum.Pt() == 0 )     secondCand      = trailCand; 
    if (leadCandNeutral.Momentum.Pt() ==0 ) leadCandNeutral = trailCand; 
    if (leadCandEM.Momentum.Pt() ==0 )      leadCandEM = trailCand; 
    if (leadCandCh.Momentum.Pt() ==0 )      leadCandCh = trailCand; 

    jetCandidate->chgEMfrac  /= jetCandidate->Momentum.E();
    jetCandidate->neuEMfrac  /= jetCandidate->Momentum.E();
    jetCandidate->chgHadfrac /= jetCandidate->Momentum.E();
    jetCandidate->neuHadfrac /= jetCandidate->Momentum.E();
     
    jetCandidate->constituents = jetCandidate->GetCandidates()->GetSize();

     
    jetCandidate->dZ = fabs(jetCandidate->Position.Z()-dynamic_cast<Candidate*>(fPVInputArray->At(0))->Position.Z());
    jetCandidate->d0 = (-jetCandidate->Position.X()*dynamic_cast<Candidate*>(fPVInputArray->At(0))->Position.Y()+jetCandidate->Position.Y()*dynamic_cast<Candidate*>(fPVInputArray->At(0))->Position.X())/jetCandidate->Position.Pt();
    

    std::sort(frac.begin(),frac.end(),std::greater<float>());
    std::sort(fracCh.begin(),fracCh.end(),std::greater<float>());
    std::sort(fracEm.begin(),fracEm.end(),std::greater<float>());
    std::sort(fracNeut.begin(),fracNeut.end(),std::greater<float>());

    jetCandidate->leadFrac = 0; jetCandidate->secondFrac = 0; jetCandidate->thirdFrac = 0; jetCandidate->fourthFrac = 0; 
    jetCandidate->leadChFrac = 0; jetCandidate->secondChFrac = 0; jetCandidate->thirdChFrac = 0; jetCandidate->fourthChFrac = 0; 
    jetCandidate->leadEmFrac = 0; jetCandidate->secondEmFrac = 0; jetCandidate->thirdEmFrac = 0; jetCandidate->fourthEmFrac = 0; 
    jetCandidate->leadNeutFrac = 0; jetCandidate->secondNeutFrac = 0; jetCandidate->thirdNeutFrac = 0; jetCandidate->fourthNeutFrac = 0; 

    assign(frac,    jetCandidate->leadFrac,    jetCandidate->secondFrac,    jetCandidate->thirdFrac,    jetCandidate->fourthFrac);
    assign(fracCh,  jetCandidate->leadChFrac,  jetCandidate->secondChFrac,  jetCandidate->thirdChFrac,  jetCandidate->fourthChFrac);
    assign(fracEm,  jetCandidate->leadEmFrac,  jetCandidate->secondEmFrac,  jetCandidate->thirdEmFrac,  jetCandidate->fourthEmFrac);
    assign(fracNeut,jetCandidate->leadNeutFrac,jetCandidate->secondNeutFrac,jetCandidate->thirdNeutFrac,jetCandidate->fourthNeutFrac);

    covMatrix(0,0) /= jetCandidate->sumPt2;
    covMatrix(0,1) /= jetCandidate->sumPt2;
    covMatrix(1,1) /= jetCandidate->sumPt2;
    covMatrix(1,0)  = covMatrix(0,1);

    jetCandidate->etaW = 0.;
    jetCandidate->phiW = 0.;
    jetCandidate->jetW = 0.;

    jetCandidate->etaW = sqrt(covMatrix(0,0));
    jetCandidate->phiW = sqrt(covMatrix(1,1));
    jetCandidate->jetW = 0.5*(jetCandidate->etaW+jetCandidate->phiW);

    TVectorD eigVals(2); 
    eigVals = TMatrixDSymEigen(covMatrix).GetEigenValues();
    jetCandidate->majW = 0.; jetCandidate->minW = 0.; jetCandidate->dRLeadCent = 0.;
    jetCandidate->majW = sqrt(fabs(eigVals(0)));
    jetCandidate->minW = sqrt(fabs(eigVals(1)));
    if(jetCandidate->majW < jetCandidate->minW ) std::swap(jetCandidate->majW,jetCandidate->minW);
    jetCandidate->dRLeadCent = jetCandidate->Momentum.DeltaR(leadCand.Momentum);
    jetCandidate->dRLead2nd  = 0.; 
    if(secondCand.Momentum.Pt() > 0) jetCandidate->dRLead2nd = jetCandidate->Momentum.DeltaR(secondCand.Momentum);
    jetCandidate->dRMean     /= jetCandidate->Momentum.Pt();
    jetCandidate->dRMeanNeut /= jetCandidate->Momentum.Pt();
    jetCandidate->dRMeanEm   /= jetCandidate->Momentum.Pt();
    jetCandidate->dRMeanCh   /= jetCandidate->Momentum.Pt();
    jetCandidate->dR2Mean    /= jetCandidate->sumPt2;

    for(size_t iCone = 0; iCone < fCones.size(); ++iCone){
        jetCandidate->FracPt[iCone]     /= jetCandidate->Momentum.Pt();
	jetCandidate->emFracPt[iCone]   /= jetCandidate->Momentum.Pt();
	jetCandidate->neutFracPt[iCone] /= jetCandidate->Momentum.Pt();
	jetCandidate->chFracPt[iCone]   /= jetCandidate->Momentum.Pt();
    }
    
    double ptMean = jetCandidate->sumPt/jetCandidate->constituents;
    double ptRMS  = 0;
    for(unsigned int i0 = 0; i0 < frac.size(); i0++) ptRMS+=(frac[i0]-ptMean)*(frac[i0]-ptMean);
    ptRMS /= jetCandidate->constituents;
    ptRMS  = sqrt(ptRMS);

    jetCandidate->ptMean = 0.; jetCandidate->ptRMS = 0.; jetCandidate->pt2A = 0.; 
    jetCandidate->ptMean   = ptMean;
    jetCandidate->ptRMS    = ptRMS/jetCandidate->Momentum.Pt();
    jetCandidate->pt2A     = sqrt( jetCandidate->ptD/jetCandidate->constituents)/jetCandidate->Momentum.Pt();
    jetCandidate->ptD      = sqrt( jetCandidate->ptD) / jetCandidate->sumPt;
    jetCandidate->ptDCh    = sqrt( jetCandidate->ptDCh) / jetCandidate->sumPtCh;
    jetCandidate->ptDNe    = sqrt( jetCandidate->ptDNe)  / jetCandidate->sumPtNe;
    jetCandidate->sumPt    = jetCandidate->sumPt;
    jetCandidate->sumChPt  = jetCandidate->sumPtCh;
    jetCandidate->sumNePt  = jetCandidate->sumPtNe;

    if( sumTkPt != 0. ) {
        jetCandidate->beta        /= sumTkPt;
        jetCandidate->betaStar    /= sumTkPt;
        jetCandidate->betaClassic /= sumTkPt;
        jetCandidate->betaClassicStar /= sumTkPt;
    } 

    float a = 0., b = 0., c = 0.;
    float ave_deta = 0., ave_dphi = 0., ave_deta2 = 0., ave_dphi2 = 0.;
    if(jetCandidate->sumPt2 > 0){
      ave_deta = sum_deta/jetCandidate->sumPt2;
      ave_dphi = sum_dphi/jetCandidate->sumPt2;
      ave_deta2 = covMatrix(0,0)/jetCandidate->sumPt2;
      ave_dphi2 = covMatrix(1,1)/jetCandidate->sumPt2;
      a = ave_deta2 - ave_deta*ave_deta;                          
      b = ave_dphi2 - ave_dphi*ave_dphi;                          
      c = -(covMatrix(0,1)/jetCandidate->sumPt2-ave_deta*ave_dphi);                
    }
    float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
    if(a+b-delta > 0) jetCandidate->axis2 = sqrt(0.5*(a+b-delta));
    else jetCandidate->axis2 = 0;

    jetCandidate->pileupIDFlagCutBased = computeCutIDflag(jetCandidate->betaClassicStar,jetCandidate->dR2Mean,fPVInputArray->GetSize(),jetCandidate->Momentum.Pt(),jetCandidate->Momentum.Eta());
    
    // fill output jetCandidate without cutting
    fOutputArray->Add(jetCandidate);    
  }
}


//------------------------------------------------------------------------------
void PileUpJetID::assign(const std::vector<float> & vec, float & a, float & b, float & c, float & d ){
  size_t sz = vec.size();
  a = ( sz > 0 ? vec[0] : 0. );
  b = ( sz > 1 ? vec[1] : 0. );
  c = ( sz > 2 ? vec[2] : 0. );
  d = ( sz > 3 ? vec[3] : 0. );
}

std::pair<int,int> PileUpJetID::getJetIdKey(float jetPt, float jetEta){
  int ptId = 0;                                                                                                                                    
  if(jetPt > 10 && jetPt < 20) ptId = 1;                                                                                 
  if(jetPt > 20 && jetPt < 30) ptId = 2;                                                                                 
  if(jetPt > 30              ) ptId = 3;                                                                                          
  int etaId = 0;                                                                                                                                   
  if(fabs(jetEta) > 2.5  && fabs(jetEta) < 2.75) etaId = 1;                                                              
  if(fabs(jetEta) > 2.75 && fabs(jetEta) < 3.0 ) etaId = 2;                                                              
  if(fabs(jetEta) > 3.0  && fabs(jetEta) < 5.0 ) etaId = 3;                                 
  return std::pair<int,int>(ptId,etaId);
}

int PileUpJetID::computeCutIDflag(float betaClassicStar,float dR2Mean,float nvtx, float jetPt, float jetEta){
  std::pair<int,int> jetIdKey = getJetIdKey(jetPt,jetEta);
  float betaStarModified      = betaClassicStar/log(nvtx-0.64);
  int idFlag(0);
  if(betaStarModified < pileUpIDCut_betaStar[0][jetIdKey.first][jetIdKey.second] and 
     dR2Mean          < pileUpIDCut_RMS[0][jetIdKey.first][jetIdKey.second] ) idFlag += 1 <<  0;
  if(betaStarModified < pileUpIDCut_betaStar[1][jetIdKey.first][jetIdKey.second] and 
     dR2Mean          < pileUpIDCut_RMS[1][jetIdKey.first][jetIdKey.second] ) idFlag += 1 <<  1;
  if(betaStarModified < pileUpIDCut_betaStar[2][jetIdKey.first][jetIdKey.second] and 
     dR2Mean          < pileUpIDCut_RMS[2][jetIdKey.first][jetIdKey.second] ) idFlag += 1 <<  2;
  return idFlag;
}
