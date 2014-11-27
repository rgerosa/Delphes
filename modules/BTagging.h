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

class BTagging: public DelphesModule {

 public:

  BTagging();
  ~BTagging();

  void Init();
  void Process();
  void Finish();

 private:

  std::map< Int_t, DelphesFormula * > fEfficiencyMapLoose; //!
  std::map< Int_t, DelphesFormula * > fEfficiencyMapMedium; //!
  std::map< Int_t, DelphesFormula * > fEfficiencyMapTight; //!

  TIterator *fItJetInputArray; //!

  const TObjArray *fJetInputArray; //!

  ClassDef(BTagging, 1)
};

#endif
