#ifndef JetFlavourAssociation_h
#define JetFlavourAssociation_h

/** \class JetFlavourAssociation
 *  Determines origin of jet,
 *  evaluate jet flavour
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
class PartonClassifier;
class LHEPartonClassifier;

class JetFlavourAssociation: public DelphesModule {

 public:

  JetFlavourAssociation();
  ~JetFlavourAssociation();

  void Init();
  void Process();
  void Finish();

  void GetAlgoFlavour   (Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray);  
  void GetPhysicsFlavour(Candidate* jet, TIter & itPartonArray, TIter & itLHEPartonArray);  

 private:

  Double_t  fDeltaR;
  
  PartonClassifier    *fClassifier; //!
  LHEPartonClassifier *fClassifierLHE; //!
  
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

  ClassDef(JetFlavourAssociation, 1)
};

#endif
