#include "puppiCleanContainer.hh"
#include "Math/SpecFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "TH2F.h"
#include "TMath.h"


using namespace std;

puppiCleanContainer::puppiCleanContainer(std::vector<RecoObj> inParticles,double iTrackerEta,bool iExperiment,bool iTuned){
    _isExperiment = iExperiment;
    _isTuned      = iTuned;
    fNeutralMinE  = 2.0;  //=> This can be tuned
    //Clear everything
    _recoParticles.resize(0);
    _pfParticles.resize(0);
    _genParticles.resize(0);
    _pfchsParticles.resize(0);
    _chargedPV.resize(0);
    _chargedNoPV.resize(0);
    _vals.resize(0);
    fTrackerEta   = iTrackerEta;
    //Link to the RecoObjects
    _recoParticles = inParticles;
    for (unsigned int i = 0; i < inParticles.size(); i++){
        fastjet::PseudoJet curPseudoJet;
        curPseudoJet.reset_PtYPhiM(inParticles[i].pt,inParticles[i].eta,inParticles[i].phi,inParticles[i].m);
        curPseudoJet.set_user_index(inParticles[i].id);
        // fill vector of pseudojets for internal references
        _pfParticles.push_back(curPseudoJet);
        if((inParticles[i].id == 0) && (inParticles[i].id == 2))  _genParticles.push_back( curPseudoJet);
        if(inParticles[i].id <= 2) _pfchsParticles.push_back(curPseudoJet);                                 //Remove Charged particles associated to other vertex
        if(inParticles[i].id == 2) _chargedPV.push_back(curPseudoJet);                                      //Take Charged particles associated to PV
        if(inParticles[i].id == 3) _chargedNoPV.push_back(curPseudoJet);
    }
    puppiWeights_chLV.resize(0);
    puppiWeights_all.resize(0);
    alphas_chLV.resize(0);
    alphas_all.resize(0);
    fNEta = 4; 
    fMed  = new double[fNEta];
    fRMS  = new double[fNEta];
    for(int i0 = 0; i0 < fNEta; i0++) fMed[i0] = 0; 
    for(int i0 = 0; i0 < fNEta; i0++) fRMS[i0] = 0; 
}
puppiCleanContainer::~puppiCleanContainer(){}

double puppiCleanContainer::goodVar(fastjet::PseudoJet &iPart,std::vector<fastjet::PseudoJet> &iParts, int iOpt) {
  double Rsub = 0.2;
  double lPup = 0;
  //Cryptic way of calling one of the Puppi Algos (there are indeed 10)
  if(iOpt == 14) return iPart.pt();
  if(iOpt == 10) return iPart.pt()+var_within_R(iOpt,iParts,iPart,Rsub);
  //if(iOpt == 10) return var_within_R(iOpt,iParts,iPart,Rsub);
  if(iOpt == 12) {
    lPup         = var_within_R(9,iParts,      iPart,Rsub);
    double lPup1 = var_within_R(9,_chargedNoPV,iPart,Rsub);
    if(lPup  != 0) lPup  = log(lPup);
    if(lPup1 != 0) lPup1 = log(lPup1);
    return lPup-lPup1;
  }
  if(iOpt == 13) { 
    lPup         = var_within_R(9,iParts,      iPart,Rsub);
    if(lPup != 0) lPup = log(lPup); 
    double Rsmall = 0.1;
    fastjet::Selector sel = fastjet::SelectorCircle(Rsmall);
    sel.set_reference(iPart);
    vector<fastjet::PseudoJet> near_particles = sel(_chargedNoPV);
    double lSum = 0; 
    for(int i0 = 0; i0 < near_particles.size(); i0++) { 
      double pPup = var_within_R(9,iParts,      near_particles[i0],Rsub);
      if(pPup != 0) pPup = log(pPup);
      lSum +=  pPup;
    }
    if(near_particles.size() > 0) lSum/=near_particles.size(); 
    lPup -= lSum;
    return lPup;
  }
  if(iOpt == 23) { 
    lPup         = var_within_R(11,iParts,      iPart,Rsub);
    double Rsmall = 0.1;
    fastjet::Selector sel = fastjet::SelectorCircle(Rsmall);
    sel.set_reference(iPart);
    vector<fastjet::PseudoJet> near_particles = sel(_chargedNoPV);
    double lSum = 0; 
    for(int i0 = 0; i0 < near_particles.size(); i0++) { 
      double pPup = var_within_R(11,iParts,      near_particles[i0],Rsub);
      lSum +=  pPup;
    }
    if(near_particles.size() > 0) lSum/=near_particles.size(); 
    lPup -= lSum;
    return lPup;
  }
  lPup = var_within_R(iOpt,iParts,iPart,Rsub);
  if(iOpt == 9 && lPup != 0) return log(lPup);
  if(iOpt == 6) lPup = lPup * iPart.pt()/pt_within_R(_pfParticles,iPart,Rsub);
  return lPup;
}
//In fact takes the median no the average
void puppiCleanContainer::getRMSAvg(int iOpt,std::vector<fastjet::PseudoJet> &iConstits,std::vector<fastjet::PseudoJet> &iParticles,double iQuant,double iPtRMS) { 
  std::vector<double>* lValsPU = new std::vector<double>[fNEta];
  int lNCount[fNEta];
  for(int i0 = 0; i0 < fNEta; i0++) lNCount[i0] = 0;
  for(unsigned int i0 = 0; i0 < iConstits.size(); i0++ ) { 
    double pVal = goodVar(iConstits[i0],iParticles,iOpt);
    _vals.push_back(pVal);
    if(std::isnan(pVal) || std::isinf(pVal)) cout << "====> Crap " << pVal << " == " << iConstits[i0].pt() << " -- " << iConstits[i0].eta() << endl;
    if(std::isnan(pVal) || std::isinf(pVal)) continue;
    lNCount[getEtaId(iConstits[i0].eta())]++;
    if(pVal == 0) continue;
    if(iConstits[i0].pt() < iPtRMS) continue;
    if( fabs(iConstits[i0].eta()) < 4.0  && iConstits[i0].user_index() % 4 != 3 ) continue;
    lValsPU[getEtaId(iConstits[i0].eta())].push_back(pVal);
  }
  for(int i0 = 0; i0 < fNEta; i0++) if(lValsPU[i0].size() > 0) std::sort (lValsPU[i0]    .begin(),lValsPU[i0]    .end());   
  for(int i0 = 0; i0 < fNEta; i0++) if(lValsPU[i0].size() > 0) fMed[i0] = lValsPU[i0]    [int(lValsPU[i0]    .size()*iQuant+0.5)]; 
  for(int i0 = 0; i0 < fNEta; i0++) fRMS[i0]     = 0; 

  for(int  i0 = 0; i0 < fNEta; i0++) { 
    for(unsigned int i1 = 0 ;i1 < lValsPU[i0]    .size(); i1++) {
      fRMS[i0]     += (lValsPU[i0]    [i1]-fMed[i0])*(lValsPU[i0][i1]-fMed[i0]);
    }
  }
  for(int i0 = 0; i0 < fNEta; i0++) fRMS[i0]/=lValsPU[i0].size();
  for(int i0 = 0; i0 < fNEta; i0++) (fRMS[i0] != 0 ) ? fRMS[i0] =sqrt(fRMS[i0]) : fRMS[i0] = 1;
}
int    puppiCleanContainer::getEtaId(float iEta) { 
  int lId = int((iEta+5.)/(10./fNEta));
  if(lId < 0      ) lId = 0; 
  if(lId > fNEta-1) lId = fNEta-1;
  return lId;
}
double puppiCleanContainer::compute(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp) {
    if(iOpt == 1 && iVal < iMed) return 0;
    if(iOpt == 1 && iVal > iMed) return 1;
    double lVal = (iVal-iMed)/iRMS;
    int lNDof = 1;
    if(iChi2Exp > 0) lVal += iChi2Exp;
    if(iChi2Exp > 0) lNDof++;
    return  TMath::Erf(lVal);//ROOT::Math::chisquared_cdf(lVal*fabs(lVal),lNDof);
}
double puppiCleanContainer::compute2(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp,double iVal1,double iMed1,double iRMS1) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal  = (iVal -iMed) /iRMS;
  double lVal1 = (iVal1-iMed1)/iRMS1;
  int lNDof = 2;
  if(iChi2Exp > 0) lVal += iChi2Exp;
  if(iChi2Exp > 0) lNDof++;
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal)+lVal1*fabs(lVal1),lNDof);
}
double puppiCleanContainer::compute3(int iOpt,double iVal,double iMed,double iRMS,double iChi2Exp,double iVal1,double iMed1,double iRMS1,double iVal2,double iMed2,double iRMS2) { 
  if(iOpt == 1 && iVal < iMed) return 0;
  if(iOpt == 1 && iVal > iMed) return 1;
  double lVal  = (iVal -iMed) /iRMS;
  double lVal1 = (iVal1-iMed1)/iRMS1;
  double lVal2 = (iVal2-iMed2)/iRMS2;
  int lNDof = 3;
  if(iChi2Exp > 0) lVal += iChi2Exp;
  if(iChi2Exp > 0) lNDof++;
  return  ROOT::Math::chisquared_cdf(lVal*fabs(lVal)+lVal1*fabs(lVal1)+lVal2*fabs(lVal2),lNDof);
}
double puppiCleanContainer::getChi2FromdZ(double iDZ) { 
  //We need to obtain prob of PU + (1-Prob of LV)
  // Prob(LV) = Gaus(dZ,sigma) where sigma = 1.5mm  (its really more like 1mm)
  double lProbLV = ROOT::Math::normal_cdf_c(fabs(iDZ),0.2)*2.; //*2 is to do it double sided
  double lProbPU = 1-lProbLV;
  if(lProbPU <= 0) lProbPU = 1e-16;   //Quick Trick to through out infs
  if(lProbPU >= 0) lProbPU = 1-1e-16; //Ditto
  double lChi2PU = TMath::ChisquareQuantile(lProbPU,1);
  return lChi2PU;
}
std::vector<fastjet::PseudoJet> puppiCleanContainer::puppiEvent     (int iOpt,double iQuant) {
    std::vector<fastjet::PseudoJet> particles;
    //Run through all compute mean and RMS
    _vals.resize(0);
    int lNEvents    = _recoParticles.size();
    
    //Log2 Puppi with Leading Vertex Particles Only
    getRMSAvg(9,_pfParticles,_chargedPV,iQuant,0.1);
    double *lMed0 = new double[fNEta];
    double *lRMS0 = new double[fNEta];
    for(int i0 = 0; i0 < fNEta; i0++) {lMed0[i0] = fMed[i0]; lRMS0[i0] = fRMS[i0];}
    
    getRMSAvg(9,_pfParticles,_pfParticles,iQuant,0.5);
    double *lMed1  = new double[fNEta];
    double *lRMS1  = new double[fNEta];
    for(int i0 = 0; i0 < fNEta; i0++) {lMed1[i0] = fMed[i0]; lRMS1[i0] = fRMS[i0];}
    
    //Sum pT in a cone 
    if(_isTuned) getRMSAvg(10,_pfParticles,_pfParticles,iQuant,0.5);
    double *lMed2  = new double[fNEta];
    double *lRMS2  = new double[fNEta];
    for(int i0 = 0; i0 < fNEta; i0++) {lMed2[i0] = fMed[i0]; lRMS2[i0] = fRMS[i0]*3.; } //3}
    
    for(int i0 = 0; i0 < lNEvents; i0++) {
	int pId = getEtaId(_pfParticles[i0].eta());
        // fill alpha values
        alphas_chLV.push_back(_vals[i0]);//-_vals[i0+3.*lNEvents]);
        alphas_all .push_back(_vals[i0+lNEvents]);
        // fill the p-values
	if(_vals[i0] == 0) _vals[i0] = lMed0[pId];
	puppiWeights_chLV.push_back( (_vals[i0]-lMed0[pId])         /lRMS0[pId]);
        puppiWeights_all .push_back( (_vals[i0+lNEvents]-lMed1[pId])/lRMS1[pId]);
        //Now compute the puppi Weight
        double pWeight = 1;
	double pChi2   = 0;
        if(_isExperiment) {
	  //Compute an Experimental Puppi Weight with delta Z info (very simple example)
	  pChi2 = getChi2FromdZ(_recoParticles[i0].dZ);
	  //Now make sure Neutrals are not set
	  if(_recoParticles[i0].pfType > 3) pChi2 = 0;
	}
	//Basic Puppi
	if(fabs(_pfParticles[i0].eta()) < fTrackerEta)      pWeight *= compute(0,_vals[i0]            ,lMed0[pId],lRMS0[pId],0.);
	if(fabs(_pfParticles[i0].eta()) > fTrackerEta)      pWeight *= compute(0,_vals[i0+lNEvents],   lMed1[pId],lRMS1[pId],0.);//,_vals[i0+2.*lNEvents],lMed2[pId],lRMS2[pId]);
     	if(_isTuned) if(fabs(_pfParticles[i0].eta()) > -1.) pWeight *= compute(0,_vals[i0+2.*lNEvents],lMed2[pId],lRMS2[pId],0.);
        //CHS
	if(_pfParticles[i0].user_index() == 2 ) pWeight = 1;
        if(_pfParticles[i0].user_index() == 3 ) pWeight = 0;
        
	//Basic Cuts
	if(std::isinf(pWeight)||std::isnan(pWeight)) continue;
	if(pWeight*_pfParticles[i0].pt()  < 0.05) continue;  //==> Elminate the low pT stuff 
	if(pWeight*_pfParticles[i0].pt()   < fNeutralMinE && _recoParticles[i0].pfType > 3 ) continue;  //threshold cut on the neutral E
        
	//Produce
	fastjet::PseudoJet curjet( pWeight*_pfParticles[i0].px(), pWeight*_pfParticles[i0].py(), pWeight*_pfParticles[i0].pz(), pWeight*_pfParticles[i0].e());
        curjet.set_user_index(_recoParticles[i0].id);
        particles.push_back(curjet);

    }
    return particles;
}


