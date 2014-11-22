#ifndef JetPileUpSubtractor_h
#define JetPileUpSubtractor_h

/** \class JetPileUpSubtractor
 *
 *  Subtract pile-up contribution from jets using the fastjet area method
 *
 *  $Date: 2012-11-18 15:57:08 +0100 (Sun, 18 Nov 2012) $
 *  $Revision: 814 $
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <deque>

class TObjArray;
class DelphesFormula;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class Selector;
}


class JetPileUpSubtractor: public DelphesModule{

public:

  JetPileUpSubtractor();
  ~JetPileUpSubtractor();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fJetPTMin;
  Bool_t   fSafe4VAreaSubtraction;

  TIterator *fItJetInputArray; //!
  TIterator *fItRhoInputArray; //!

  const TObjArray *fJetInputArray; //!
  const TObjArray *fRhoInputArray; //!

  TObjArray *fOutputArray; //!

  // only for Safe4V subtraction
  void *fPluginRho; //!                                                                                                                                                                   
  void *fPlugin; //!                                                                                                                                                                   

  fastjet::JetDefinition *fDefinition; //!                                                                                                                                        
  fastjet::JetDefinition *fDefinitionRho; //!                                                                                                                                   
 
  Int_t    fJetAlgorithmRho;
  Int_t    fJetAlgorithm;

  Double_t fParameterRRho;
  Double_t fParameterR;

  Double_t fConeRadiusRho;
  Double_t fConeRadius;

  Double_t fSeedThreshold;
  Double_t fSeedThresholdRho;

  Double_t fConeAreaFraction;
  Double_t fConeAreaFractionRho;

  Double_t fAdjacencyCutRho;
  Double_t fAdjacencyCut;

  Double_t fOverlapThreshold;
  Double_t fOverlapThresholdRho;
  
  // --- FastJet Area method --------                                                                                                                                                    
  fastjet::AreaDefinition *fAreaDefinition;
  fastjet::AreaDefinition *fAreaDefinitionRho;
  Int_t  fAreaAlgorithm;
  Int_t  fAreaAlgorithmRho;

  // -- ghost based areas --                                                                                                                                                      
  Double_t fGhostEtaMax;
  Double_t fGhostEtaMaxRho;
  Int_t    fRepeat;
  Int_t    fRepeatRho;
  Int_t    fMaxIterations;
  Int_t    fMaxIterationsRho;
  Int_t    fMaxPairSize;
  Int_t    fMaxPairSizeRho;
  Int_t    fIratch;
  Int_t    fIratchRho;
  Double_t fGhostArea;
  Double_t fGhostAreaRho;
  Double_t fGridScatter;
  Double_t fGridScatterRho;
  Double_t fPtScatter;
  Double_t fPtScatterRho;
  Double_t fMeanGhostPt;
  Double_t fMeanGhostPtRho;

  // -- voronoi areas --                                                                                                                                                                   
  Double_t fEffectiveRfact;
  Double_t fEffectiveRfactRho;
  std::map< Double_t, Double_t > fEtaRangeMap; //!                                                                                                                                        

  TIterator *fItInputArray; //!                                                                                                                                                           
       
  const TObjArray *fInputArray; //!                                                                                                                                                        
        

  ClassDef(JetPileUpSubtractor, 1)
};

#endif
