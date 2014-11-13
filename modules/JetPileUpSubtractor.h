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
  void *fPlugin; //!                                                                                                                                                                   
  fastjet::JetDefinition *fDefinition; //!                                                                                                                                                  
  Int_t fJetAlgorithm;
  Double_t fParameterR;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Double_t fAdjacencyCut;
  Double_t fOverlapThreshold;
  
  // --- FastJet Area method --------                                                                                                                                                    
  fastjet::AreaDefinition *fAreaDefinition;
  Int_t  fAreaAlgorithm;

  // -- ghost based areas --                                                                                                                                                      
  Double_t fGhostEtaMax;
  Int_t    fRepeat;
  Int_t    fMaxIterations;
  Int_t    fMaxPairSize;
  Int_t    fIratch;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt;

  // -- voronoi areas --                                                                                                                                                                   
  Double_t fEffectiveRfact;
  std::map< Double_t, Double_t > fEtaRangeMap; //!                                                                                                                                        
  TIterator *fItInputArray; //!                                                                                                                                                           
  const TObjArray *fInputArray; //!                                                                                                                                                      


  

  ClassDef(JetPileUpSubtractor, 1)
};

#endif
