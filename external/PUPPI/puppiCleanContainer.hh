#include "PUPPI/RecoObj.hh"
#include "PUPPI/puppiParticle.hh"
#include "PUPPI/puppiAlgoBin.hh"

#include "fastjet/internal/base.hh"
#include "fastjet/PseudoJet.hh"
#include <algorithm>

using namespace std;

//......................

class puppiCleanContainer{

 public:

  // basic constructor which takes input particles as RecoObj, the tracker eta extension and other two boolean info
  puppiCleanContainer(std::vector<RecoObj> inParticles,    // incoming particles of the event
                      std::vector<puppiAlgoBin> puppiAlgo, // vector with the definition of the puppi algorithm in different eta region (one for each eta) 
                      float minPuppiWeight  = 0.01,        // min puppi weight cut
                      bool  useExp = false                 // useDz vertex probability
  ); 

  ~puppiCleanContainer(); 

  // ----- get methods

  // get all the PF particles
  std::vector<fastjet::PseudoJet> pfParticles()  { return  fPFParticles_; }    
  // get all the PF charged from PV
  std::vector<fastjet::PseudoJet> pvParticles()  { return  fChargedPV_; }        
  // get all the PF charged from PU
  std::vector<fastjet::PseudoJet> puParticles()  { return  fChargedNoPV_; }    
  // get CHS particle collection
  std::vector<fastjet::PseudoJet> pfchsParticles(){ return fPFchsParticles_; }    
  // get puppi weight for all particles 
  std::vector<float> getPuppiWeights() { return fPuppiWeights_; };

  // process puppi
  std::vector<fastjet::PseudoJet> puppiEvent();
 
 protected:

   void    getRMSAvg(int iPuppiAlgo, std::vector<fastjet::PseudoJet> & particlesAll, std::vector<fastjet::PseudoJet> &chargedPV);        
   float   goodVar  (const fastjet::PseudoJet & particle, const std::vector<fastjet::PseudoJet> & particleAll, const int & pPupId, const float & coneSize);    
   void    computeMedRMS(const int & puppiAlgo);  
   float   compute(const float & val, const float & chi2, const puppiAlgoBin & puppiAlgo);

   // some get functions
   float getNeutralPtCut(const float&, const float&, const int&);
   int   getPuppiId(const float & pt, const float & eta, const std::vector<puppiAlgoBin> & puppiAlgos);
   float getChi2FromdZ(float iDZ);

   // other functions
   float  var_within_R(const int & pPupId, const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R);
   float  pt_within_R(const std::vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet & centre, const float & R);
   fastjet::PseudoJet flow_within_R(const vector<fastjet::PseudoJet> & particles, const fastjet::PseudoJet& centre, const float & R);
   
  
 private:    
 
  std::vector<RecoObj>            fRecoParticles_;
  std::vector<fastjet::PseudoJet> fPFParticles_;
  std::vector<fastjet::PseudoJet> fPFchsParticles_;    
  std::vector<fastjet::PseudoJet> fChargedPV_;
  std::vector<fastjet::PseudoJet> fChargedNoPV_;

  std::vector<puppiAlgoBin> puppiAlgo_;
  std::vector<float> fPuppiWeights_;

  int    fNPV_;             // NPV
  float  fMinPuppiWeight_;
  float  fPVFrac_;
  bool   fUseExp_ ;
    
};

