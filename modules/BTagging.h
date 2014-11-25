#ifndef BTagging_h
#define BTagging_h

/** \class BTagging
 *  Determines origin of jet,
 *  applies b-tagging efficiency (miss identification rate) formulas
 *  and sets b-tagging flags 
 *  $Date: 2013-04-26 12:39:14 +0200 (Fri, 26 Apr 2013) $
 *  $Revision: 1099 $
 *  \author P. Demin - UCL, Louvain-la-Neuve
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <map>

class TObjArray;
class DelphesFormula;

class ExRootFilter;
class BTaggingPartonClassifier;
class BTaggingLHEPartonClassifier;

class BTagging: public DelphesModule {

 public:

  BTagging();
  ~BTagging();

  void Init();
  void Process();
  void Finish();

  void GetAlgoFlavour(Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray);  
  void GetPhysicsFlavour(Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray);  

 private:

  Int_t fBitNumber;
  Double_t fDeltaR;
  std::map< Int_t, DelphesFormula * > fEfficiencyMap; //!
  
  BTaggingPartonClassifier *fClassifier; //!
  BTaggingLHEPartonClassifier *fClassifierLHE; //!
  
  ExRootFilter *fFilter;
  ExRootFilter *fFilterLHE;

  TIterator *fItPartonInputArray; //!  
  TIterator *fItLHEPartonInputArray; //!  
  TIterator *fItJetInputArray; //!
  TIterator *fItParticleInputArray; //!

  const TObjArray *fPartonInputArray; //! 
  const TObjArray *fLHEPartonInputArray; //! 
  const TObjArray *fJetInputArray; //!
  const TObjArray *fParticleInputArray; //!

  ClassDef(BTagging, 1)
};

#endif
